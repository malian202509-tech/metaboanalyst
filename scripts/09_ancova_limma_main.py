"""09. ANCOVA / limma 主分析: 双轨 BMI_group 差异分析.

模型: log2(metab) ~ BMI_group + age + GA_decimal + GWG + GDM
强制协变量依据 (见 docs/数据处理与分析进展.md §5.6):
  - 孕妇年龄 SMD=0.21 (⚠) → 必入
  - 孕周 GA_decimal SMD=0.20 (⚠) → 必入 (兼吸收 GA_category 的不均衡)
  - 孕期增重 GWG SMD=0.15 (⚠) → 必入 (5 例缺失 → complete case)
  - GDM SMD=0.05 (√) 但机制相关 → 稳健性考虑保留
  - PIH SMD=0.56 (❗❗) 但 n=8 → 不入模型, 单独作临床发现报告
  - 出生体重 / 出生体重分类 → 中介, 禁止调整
  - 年份 → 已在 ComBat 阶段抹掉

统计:
  - OLS + HC3 稳健 SE → β_BMI (log2FC) + 95% CI + p_ols
  - limma EB 方差收缩 (Smyth 2004) → moderated t, p_limma
  - Mann-Whitney U → p_wilcoxon (非参数稳健性参照)
  - MAD 离群 (|x-median|/MAD>3.5) → 剔除后重跑 OLS → is_robust_to_outlier

三层入选 (互斥, 取最高, 基于 limma q 因为 limma 是主分析):
  - strict      : q_limma < 0.05
  - standard    : q_limma < 0.10 且 |log2FC| ≥ log2(1.2)
  - exploratory : p_limma < 0.05 且 |log2FC| ≥ log2(1.5)

high_confidence (SOP "两条路取交集"): tier ∈ {strict, standard}
  且 OLS+HC3 同时支持 (raw p<0.05 且 95% CI 不跨 0).
  用于剔除 limma 在小特征数下因方差池化产生的伪阳性.

输入:
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat.xlsx
  data/01_raw/分娩特征分析_2019-2021_主分析队列n165.xlsx
  data/02_preprocessed/sample_alignment_n165.csv

输出:
  results/tables/ancova_main_{80,50}.csv
  results/tables/ancova_audit.csv
"""
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
from utils.limma import moderated_t, bh_qvalues

COMBAT_DIR = ROOT / 'data' / '03_batch_corrected'
RAW_DIR = ROOT / 'data' / '01_raw'
ALIGN_DIR = ROOT / 'data' / '02_preprocessed'
OUT_DIR = ROOT / 'results' / 'tables'
OUT_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'

# 设计矩阵列: [截距, BMI_overweight (vs 正常), age, GA_decimal, GWG, GDM]
CONTRAST_COL = 'BMI_overweight'
COVARS = ['age', 'GA_decimal', 'GWG', 'GDM']

FC_STD_LOG2 = np.log2(1.2)      # 标准层 FC 门槛
FC_EXP_LOG2 = np.log2(1.5)      # 探索层 FC 门槛
MAD_K = 3.5                     # MAD 离群阈值


def load_covariates():
    """加载临床表 + sample_alignment, 返回以 omx_id 为索引的协变量表."""
    align = pd.read_csv(ALIGN_DIR / 'sample_alignment_n165.csv')
    clin = pd.read_excel(RAW_DIR / '分娩特征分析_2019-2021_主分析队列n165.xlsx')

    cov = (clin[['标本号', '孕妇年龄', 'GA_decimal', '孕期增重', 'GDM', 'BMI_group']]
           .rename(columns={'孕妇年龄': 'age', '孕期增重': 'GWG'})
           .merge(align[['标本号', 'omx_id']], on='标本号', how='left'))

    # 二分类编码: 正常=0, 超重肥胖=1
    cov[CONTRAST_COL] = (cov['BMI_group'] == '超重肥胖').astype(float)
    cov['GDM'] = cov['GDM'].astype(float)
    return cov.set_index('omx_id')


def detect_outliers_mad(y, k=MAD_K):
    """MAD 离群检测: |y-median|/MAD > k. 返回 bool 掩码 (True=离群)."""
    med = np.median(y)
    mad = stats.median_abs_deviation(y, scale='normal')  # scale='normal' 使 MAD 等价 SD
    if mad <= 0:
        return np.zeros_like(y, dtype=bool)
    return np.abs(y - med) / mad > k


def fit_per_feature(y, X, contrast_idx):
    """对单个代谢物拟合 OLS+HC3, 返回 (β, se_hc3, p_hc3, ci_lo, ci_hi, s2_ols, df_resid)."""
    model = sm.OLS(y, X).fit(cov_type='HC3')
    beta = float(model.params.iloc[contrast_idx])
    se_hc3 = float(model.bse.iloc[contrast_idx])
    p_hc3 = float(model.pvalues.iloc[contrast_idx])
    ci_lo, ci_hi = model.conf_int().iloc[contrast_idx].tolist()
    # 限定方差用 OLS 残差 (HC3 不适合 limma EB 假设的同方差)
    s2_ols = float(np.sum(model.resid ** 2) / model.df_resid)
    return beta, se_hc3, p_hc3, float(ci_lo), float(ci_hi), s2_ols, int(model.df_resid)


def refit_robust(y, X, mask_keep, contrast_idx):
    """剔除离群后重跑 OLS, 返回 p_robust + 剔除后 n. 若不能拟合返回 NaN."""
    n_keep = int(mask_keep.sum())
    if n_keep < X.shape[1] + 5:
        return np.nan, n_keep
    y_k = y[mask_keep]
    X_k = X.loc[mask_keep] if hasattr(X, 'loc') else X[mask_keep]
    # 同组至少各 3 例才有意义
    g_k = X_k[CONTRAST_COL].values
    if int((g_k == 1).sum()) < 3 or int((g_k == 0).sum()) < 3:
        return np.nan, n_keep
    try:
        m = sm.OLS(y_k, X_k).fit(cov_type='HC3')
        return float(m.pvalues.iloc[contrast_idx]), n_keep
    except Exception:
        return np.nan, n_keep


def assign_tier(q_limma, p_limma, log2fc):
    """三层互斥分级 (基于 limma q, limma 是主分析)."""
    afc = abs(log2fc)
    if pd.notna(q_limma) and q_limma < 0.05:
        return 'strict'
    if pd.notna(q_limma) and q_limma < 0.10 and afc >= FC_STD_LOG2:
        return 'standard'
    if pd.notna(p_limma) and p_limma < 0.05 and afc >= FC_EXP_LOG2:
        return 'exploratory'
    return '-'


def is_high_confidence(tier, p_ols, ci_lo, ci_hi):
    """高置信 = limma 主分析达 standard/strict 且 OLS+HC3 同时支持."""
    if tier not in {'strict', 'standard'}:
        return False
    if not (pd.notna(p_ols) and p_ols < 0.05):
        return False
    return (ci_lo > 0 and ci_hi > 0) or (ci_lo < 0 and ci_hi < 0)


def process_track(src_name, tag, cov_full):
    print(f'\n{"=" * 60}')
    print(f'轨道 {tag}: {src_name}')

    df = pd.read_excel(COMBAT_DIR / src_name)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    matrix = df[sample_cols].copy()
    n_feat, n_samp = matrix.shape
    print(f'  原始形状: {n_feat} 代谢物 × {n_samp} 样本')

    # 对齐协变量到组学样本顺序
    cov = cov_full.reindex(sample_cols)
    missing_cov = cov[[CONTRAST_COL] + COVARS].isna().any(axis=1)
    print(f'  协变量缺失样本: {int(missing_cov.sum())} '
          f'(主要 GWG: {int(cov["GWG"].isna().sum())})')

    keep_cols = [c for c, m in zip(sample_cols, missing_cov) if not m]
    cov_cc = cov.loc[keep_cols]
    matrix_cc = matrix[keep_cols]
    n_cc = len(keep_cols)
    n_normal = int((cov_cc[CONTRAST_COL] == 0).sum())
    n_over = int((cov_cc[CONTRAST_COL] == 1).sum())
    print(f'  Complete-case 后: n={n_cc} (正常 {n_normal} + 超重肥胖 {n_over})')

    # 设计矩阵
    X = sm.add_constant(cov_cc[[CONTRAST_COL] + COVARS].astype(float))
    contrast_idx = list(X.columns).index(CONTRAST_COL)
    XtX_inv = np.linalg.pinv(X.values.T @ X.values)
    se_ord_factor = float(XtX_inv[contrast_idx, contrast_idx])

    # 逐特征拟合
    rows = []
    s2_all, beta_all = [], []
    df_resid_all = None
    for feat_idx in range(n_feat):
        y = matrix_cc.iloc[feat_idx].values.astype(float)
        meta = df.loc[feat_idx, META_COLS].to_dict()

        beta, se_hc3, p_hc3, ci_lo, ci_hi, s2_ols, df_resid = fit_per_feature(y, X, contrast_idx)
        df_resid_all = df_resid  # 同 X → 所有特征一致

        # 分组均值 (log2 尺度)
        g = cov_cc[CONTRAST_COL].values
        mean_normal = float(np.mean(y[g == 0]))
        mean_over = float(np.mean(y[g == 1]))

        # Wilcoxon (Mann-Whitney U)
        try:
            w_stat, p_wil = stats.mannwhitneyu(y[g == 1], y[g == 0], alternative='two-sided')
            p_wil = float(p_wil)
        except ValueError:
            p_wil = np.nan

        # MAD 离群 + 稳健重拟合
        out_mask = detect_outliers_mad(y)
        n_out = int(out_mask.sum())
        p_robust, n_after = refit_robust(y, X, ~out_mask, contrast_idx)

        rows.append({
            **meta,
            'n': n_cc,
            'mean_log2_normal': round(mean_normal, 4),
            'mean_log2_overweight': round(mean_over, 4),
            'log2FC': round(beta, 4),
            'log2FC_CI_lo': round(ci_lo, 4),
            'log2FC_CI_hi': round(ci_hi, 4),
            'p_ols_hc3': p_hc3,
            'p_limma': np.nan,        # 占位, 跨特征算完再填
            'p_wilcoxon': p_wil,
            'n_outliers': n_out,
            'n_after_outlier_removal': n_after,
            'p_robust_ols': p_robust,
            'is_robust_to_outlier': bool(np.isfinite(p_robust) and p_robust < 0.05),
        })
        s2_all.append(s2_ols)
        beta_all.append(beta)

    # limma EB 跨特征
    s2_arr = np.array(s2_all)
    beta_arr = np.array(beta_all)
    mod = moderated_t(beta_arr, se_ord_factor, s2_arr, df=df_resid_all)
    print(f'  limma 先验: d_0={mod["df_prior"]:.2f}, s2_0={mod["s2_prior"]:.4f}, '
          f'df_total={float(np.atleast_1d(mod["df_total"])[0]):.2f}')

    # 回填 limma + BH
    for i, r in enumerate(rows):
        r['p_limma'] = float(mod['p_mod'][i])
        r['t_limma'] = round(float(mod['t_mod'][i]), 4)

    tbl = pd.DataFrame(rows)
    tbl['q_ols_hc3_bh'] = bh_qvalues(tbl['p_ols_hc3'].values)
    tbl['q_limma_bh'] = bh_qvalues(tbl['p_limma'].values)
    tbl['q_wilcoxon_bh'] = bh_qvalues(tbl['p_wilcoxon'].values)

    # 三层入选 (基于 limma q) + 高置信交集标志
    tbl['tier'] = [
        assign_tier(ql, pl, fc)
        for ql, pl, fc in zip(tbl['q_limma_bh'], tbl['p_limma'], tbl['log2FC'])
    ]
    tbl['high_confidence'] = [
        is_high_confidence(t, p, lo, hi)
        for t, p, lo, hi in zip(
            tbl['tier'], tbl['p_ols_hc3'],
            tbl['log2FC_CI_lo'], tbl['log2FC_CI_hi'],
        )
    ]
    tbl['direction'] = np.where(tbl['log2FC'] > 0, '↑(超重肥胖>正常)', '↓(超重肥胖<正常)')

    # P / q 格式化
    for c in ['p_ols_hc3', 'p_limma', 'p_wilcoxon', 'p_robust_ols',
              'q_ols_hc3_bh', 'q_limma_bh', 'q_wilcoxon_bh']:
        tbl[c] = tbl[c].apply(lambda x: round(float(x), 6) if pd.notna(x) else np.nan)

    # 排序: tier 先, 再 q_limma
    tier_rank = {'strict': 0, 'standard': 1, 'exploratory': 2, '-': 3}
    tbl = (tbl.assign(_tr=tbl['tier'].map(tier_rank))
              .sort_values(['_tr', 'q_limma_bh', 'p_limma'])
              .drop(columns='_tr')
              .reset_index(drop=True))

    # 列顺序
    front = ['Metabolite Name', 'Chinese Name', 'Class', 'Sub Class',
             'Super Class', 'Direct Parent', 'KEGG ID', 'HMDB ID',
             'Retention Time(min)', 'Cellular Locations',
             'n', 'tier', 'high_confidence', 'direction',
             'log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi',
             'mean_log2_normal', 'mean_log2_overweight',
             'p_ols_hc3', 'q_ols_hc3_bh',
             'p_limma', 't_limma', 'q_limma_bh',
             'p_wilcoxon', 'q_wilcoxon_bh',
             'n_outliers', 'n_after_outlier_removal',
             'p_robust_ols', 'is_robust_to_outlier']
    tbl = tbl[front]

    out_file = OUT_DIR / f'ancova_main_{tag}.csv'
    tbl.to_csv(out_file, index=False, encoding='utf-8-sig')
    print(f'  已写入: {out_file.relative_to(ROOT)}')

    # 方差诊断 (评估 limma EB 是否退化为池化方差)
    s2_min, s2_med, s2_max = float(s2_arr.min()), float(np.median(s2_arr)), float(s2_arr.max())
    s2_var_ratio = s2_max / s2_min if s2_min > 0 else np.inf
    log_s2_var = float(np.var(np.log(s2_arr), ddof=1))

    # 摘要
    n_strict = int((tbl['tier'] == 'strict').sum())
    n_standard = int((tbl['tier'] == 'standard').sum())
    n_exp = int((tbl['tier'] == 'exploratory').sum())
    n_robust = int(tbl['is_robust_to_outlier'].sum())
    n_high_conf = int(tbl['high_confidence'].sum())
    print(f'  入选: strict={n_strict}, standard={n_standard}, exploratory={n_exp}')
    print(f'  high_confidence (limma+OLS+CI 三认同): {n_high_conf}')
    print(f'  离群稳健 (剔 outlier 后仍 p<0.05): {n_robust} / {n_feat}')
    print(f'  特征方差: s²∈[{s2_min:.3f}, {s2_max:.3f}], median {s2_med:.3f}, '
          f'max/min ratio {s2_var_ratio:.1f}')

    d0_str = (f'{float(mod["df_prior"]):.2f}' if np.isfinite(mod['df_prior'])
              and float(mod['df_prior']) < 1e6 else '≈∞ (特征方差近似一致, EB 退化为池化方差)')
    return {
        'tier_track': tag,
        'src': src_name,
        'n_feat': n_feat,
        'n_sample_cc': n_cc,
        'n_normal': n_normal,
        'n_overweight': n_over,
        'df_resid': df_resid_all,
        'limma_d0': d0_str,
        'limma_s2_0': round(float(mod['s2_prior']), 4),
        's2_min': round(s2_min, 4),
        's2_median': round(s2_med, 4),
        's2_max': round(s2_max, 4),
        's2_max_min_ratio': round(s2_var_ratio, 2),
        'log_s2_var': round(log_s2_var, 4),
        'n_strict_qlimma05': n_strict,
        'n_standard_qlimma10_fc1.2': n_standard,
        'n_exploratory_plimma05_fc1.5': n_exp,
        'n_high_confidence': n_high_conf,
        'n_outlier_robust': n_robust,
        'out_file': out_file.name,
    }


def main():
    print('=== 09. ANCOVA / limma 主分析 (双轨) ===')
    print('模型: log2(metab) ~ BMI_group + age + GA_decimal + GWG + GDM')

    cov_full = load_covariates()
    print(f'  协变量表: {len(cov_full)} 样本, '
          f'BMI={cov_full["BMI_group"].value_counts().to_dict()}, '
          f'GWG 缺失 {int(cov_full["GWG"].isna().sum())}, '
          f'GDM 阳性 {int((cov_full["GDM"] == 1).sum())}')

    records = []
    records.append(process_track('ori_n165_filtered80_log2_combat.xlsx', '80', cov_full))
    records.append(process_track('ori_n165_filtered50_log2_combat.xlsx', '50', cov_full))

    audit = pd.DataFrame(records)
    audit_path = OUT_DIR / 'ancova_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
