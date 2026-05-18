"""09. ANCOVA / limma 主分析: 双轨 BMI_group 差异分析 (方案 D 矩阵 + 二轮协变量).

模型: log2(metab) ~ BMI_group + age + GA_decimal + year (B20 ref, 哑变量 A19/C21)

协变量决策 (2026-05-18 二轮锁定, 见 docs/数据处理与分析进展.md §10 & memory project-pipeline-state):
  - age          : Table 1 SMD=0.21 (⚠) → 必入
  - GA_decimal   : Table 1 SMD=0.20 (⚠) → 必入 (兼吸收 GA_category 不均衡)
  - year         : 方案 D 切换后新增 (年代/季节/储存复合效应, 与 batch 独立; PERMANOVA p=0.001)
  - 删除 GWG    : mediator-like, 写入敏感性分析
  - 删除 GDM    : Table 1 SMD=0.05 √, 远低 0.10 阈值
  - 删除 产次   : 队列未采集 → Limitations
  - 删除 分娩方式: Table 1 (病历补齐后) SMD=0.036 √, 不入模型
  - 不入 PIH    : SMD=0.557 但 n=8, 单独作临床发现报告
  - 不入 出生体重/分类: 中介 (暴露→中介→结局), 调整会洗效应
  - 不入 batch  : 已在 ComBat 阶段抹掉 (方案 D)

样本: n=165 全保留 (无 complete case 剔除, 因 4 个协变量均无缺失)
设计矩阵参数 = 6 (截距 + BMI + age + GA + year_A19 + year_C21), df_resid = 159

统计:
  - OLS + HC3 稳健 SE → β_BMI (log2FC) + 95% CI + p_ols
  - limma EB 方差收缩 (Smyth 2004) → moderated t, p_limma
  - Mann-Whitney U → p_wilcoxon (非参数稳健性参照)
  - OLS internal studentized residual → |r|>3 标记离群 (而非 y 的 MAD, 修正分布偏态误识)
  - 离群剔除后重拟合: 多列稳健性指标 (方向一致 + 效应衰减 + p<0.05)

三层入选 (互斥, 取最高, 基于 limma q):
  - strict      : q_limma < 0.05
  - standard    : q_limma < 0.10 且 |log2FC| ≥ log2(1.2)
  - exploratory : p_limma < 0.05 且 |log2FC| ≥ log2(1.5)

high_confidence: tier ∈ {strict, standard} 且 OLS+HC3 raw p<0.05 且 95% CI 不跨 0
is_robust_to_outlier (新版): 仅当 n_outliers>0 且 (方向一致 AND |效应衰减|≤50% AND p_robust<0.05) 才为 True;
                            n_outliers=0 时直接标 True (无离群无需检验)

输入:
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat.xlsx
  data/01_raw/分娩特征分析_2019-2021_主分析队列n165.xlsx
  data/02_preprocessed/sample_alignment_n165.csv (含 omx_id → year 派生)

输出:
  results/tables/ancova_main_{80,50}.csv  (50 轨多一列 also_in_80_main)
  results/tables/ancova_audit.csv
"""
import os
import sys
import tempfile
import numpy as np
# Windows PowerShell 默认 GBK, 重配 stdout 为 UTF-8 以支持中文与上下标
try:
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass
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

CONTRAST_COL = 'BMI_overweight'
# year 双哑变量 (B20 ref, 因 n=109 最大); 与 ComBat 内部编码顺序无关
YEAR_DUMMIES = ['year_A19', 'year_C21']
COVARS = ['age', 'GA_decimal'] + YEAR_DUMMIES

FC_STD_LOG2 = np.log2(1.2)
FC_EXP_LOG2 = np.log2(1.5)
STUD_RESID_THRESH = 3.0          # |internal studentized residual| > 3 → 离群
ATTENUATION_MAX = 0.50           # 效应衰减 >50% 视为不稳健


def load_covariates():
    """加载临床表 + sample_alignment, 返回以 omx_id 为索引的协变量表 (含 year 双哑变量)."""
    align = pd.read_csv(ALIGN_DIR / 'sample_alignment_n165.csv')
    clin = pd.read_excel(RAW_DIR / '分娩特征分析_2019-2021_主分析队列n165.xlsx')

    cov = (clin[['标本号', '孕妇年龄', 'GA_decimal', 'BMI_group']]
           .rename(columns={'孕妇年龄': 'age'})
           .merge(align[['标本号', 'omx_id']], on='标本号', how='left'))

    cov[CONTRAST_COL] = (cov['BMI_group'] == '超重肥胖').astype(float)

    # year 派生 (omx_id 前缀): A19=2019, B20=2020 (ref), C21=2021
    year_prefix = cov['omx_id'].str[:3]
    cov['year_A19'] = (year_prefix == 'A19').astype(float)
    cov['year_C21'] = (year_prefix == 'C21').astype(float)

    return cov.set_index('omx_id')


def detect_outliers_studentized(model, thresh=STUD_RESID_THRESH):
    """基于 OLS internal studentized residual 的离群检测.

    比 y 的 MAD 更标准: 残差去除了协变量与组间结构的影响,
    避免把分布偏态或组间效应误识为离群 (修正旧版 50 轨 n_outliers=31 的虚高).
    """
    influence = model.get_influence()
    r_stud = influence.resid_studentized_internal
    return np.abs(r_stud) > thresh


def fit_per_feature(y, X, contrast_idx):
    """OLS+HC3 主拟合; 返回 (β, se_hc3, p_hc3, ci_lo, ci_hi, s2_ols, df_resid, model_ols).

    model_ols 是无 HC3 的版本 (供 studentized residual 用; HC3 cov_type 不暴露 influence API)."""
    model_hc3 = sm.OLS(y, X).fit(cov_type='HC3')
    beta = float(model_hc3.params.iloc[contrast_idx])
    se_hc3 = float(model_hc3.bse.iloc[contrast_idx])
    p_hc3 = float(model_hc3.pvalues.iloc[contrast_idx])
    ci_lo, ci_hi = model_hc3.conf_int().iloc[contrast_idx].tolist()
    s2_ols = float(np.sum(model_hc3.resid ** 2) / model_hc3.df_resid)
    df_resid = int(model_hc3.df_resid)
    # 普通 OLS 仅供 influence (β 完全等价于 HC3 的 β, 只是 SE 不同)
    model_ord = sm.OLS(y, X).fit()
    return beta, se_hc3, p_hc3, float(ci_lo), float(ci_hi), s2_ols, df_resid, model_ord


def refit_after_outlier(y, X, mask_keep, contrast_idx, beta_full):
    """剔除离群后重拟合, 返回多列稳健性指标.

    Returns dict with keys:
      n_after, beta_after, p_robust_ols, direction_consistent, attenuation_pct
    """
    n_keep = int(mask_keep.sum())
    base = {
        'n_after': n_keep, 'beta_after': np.nan, 'p_robust_ols': np.nan,
        'direction_consistent': False, 'attenuation_pct': np.nan,
    }
    if n_keep < X.shape[1] + 5:
        return base
    y_k = y[mask_keep]
    X_k = X.loc[mask_keep] if hasattr(X, 'loc') else X[mask_keep]
    g_k = X_k[CONTRAST_COL].values
    if int((g_k == 1).sum()) < 3 or int((g_k == 0).sum()) < 3:
        return base
    try:
        m = sm.OLS(y_k, X_k).fit(cov_type='HC3')
        beta_after = float(m.params.iloc[contrast_idx])
        p_robust = float(m.pvalues.iloc[contrast_idx])
        same_sign = (np.sign(beta_after) == np.sign(beta_full)) and abs(beta_full) > 0
        if abs(beta_full) > 0:
            atten = (abs(beta_full) - abs(beta_after)) / abs(beta_full)  # >0 = 衰减
        else:
            atten = np.nan
        return {
            'n_after': n_keep, 'beta_after': beta_after, 'p_robust_ols': p_robust,
            'direction_consistent': bool(same_sign),
            'attenuation_pct': float(atten) if pd.notna(atten) else np.nan,
        }
    except Exception:
        return base


def is_robust_to_outlier(n_out, robust_info):
    """新版稳健性判定 (修正旧版仅看 p<0.05 的简化逻辑).

    - n_outliers=0: 无离群可剔, 自动稳健 (True)
    - n_outliers>0: 必须 方向一致 + |衰减|≤50% + p_robust<0.05
    """
    if n_out == 0:
        return True
    if not robust_info['direction_consistent']:
        return False
    atten = robust_info['attenuation_pct']
    if pd.isna(atten) or atten > ATTENUATION_MAX:
        return False
    p_r = robust_info['p_robust_ols']
    return bool(pd.notna(p_r) and p_r < 0.05)


def assign_tier(q_limma, p_limma, log2fc):
    afc = abs(log2fc)
    if pd.notna(q_limma) and q_limma < 0.05:
        return 'strict'
    if pd.notna(q_limma) and q_limma < 0.10 and afc >= FC_STD_LOG2:
        return 'standard'
    if pd.notna(p_limma) and p_limma < 0.05 and afc >= FC_EXP_LOG2:
        return 'exploratory'
    return '-'


def is_high_confidence(tier, p_ols, ci_lo, ci_hi):
    if tier not in {'strict', 'standard'}:
        return False
    if not (pd.notna(p_ols) and p_ols < 0.05):
        return False
    return (ci_lo > 0 and ci_hi > 0) or (ci_lo < 0 and ci_hi < 0)


def atomic_write_csv(df, target):
    """原子写: 先写临时文件, 再 os.replace. 避免 Excel/WPS 锁文件时直接 crash."""
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=target.stem + '_', suffix='.csv.tmp', dir=str(target.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, index=False, encoding='utf-8-sig')
        try:
            os.replace(tmp, target)
        except PermissionError as e:
            os.unlink(tmp)
            raise PermissionError(
                f'无法写入 {target.name}: 文件可能正被 Excel/WPS 打开. '
                f'请关闭后重跑. (原始错误: {e})'
            ) from e
    except Exception:
        if os.path.exists(tmp):
            try: os.unlink(tmp)
            except OSError: pass
        raise


def process_track(src_name, tag, cov_full, main80_features=None):
    print(f'\n{"=" * 60}')
    print(f'轨道 {tag}: {src_name}')

    df = pd.read_excel(COMBAT_DIR / src_name)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    matrix = df[sample_cols].copy()
    n_feat, n_samp = matrix.shape
    print(f'  原始形状: {n_feat} 代谢物 × {n_samp} 样本')

    cov = cov_full.reindex(sample_cols)
    missing = cov[[CONTRAST_COL] + COVARS].isna().any(axis=1)
    n_miss = int(missing.sum())
    if n_miss > 0:
        print(f'  ⚠ 协变量缺失样本: {n_miss} (理论应为 0, 请核对)')
        keep_cols = [c for c, m in zip(sample_cols, missing) if not m]
    else:
        keep_cols = sample_cols
    cov_cc = cov.loc[keep_cols]
    matrix_cc = matrix[keep_cols]
    n_cc = len(keep_cols)
    n_normal = int((cov_cc[CONTRAST_COL] == 0).sum())
    n_over = int((cov_cc[CONTRAST_COL] == 1).sum())
    print(f'  分析样本: n={n_cc} (正常 {n_normal} + 超重肥胖 {n_over})')
    year_dist = cov_cc[YEAR_DUMMIES].sum().to_dict()
    year_dist['B20_ref'] = n_cc - int(sum(year_dist.values()))
    print(f'  year 分布: {year_dist}')

    X = sm.add_constant(cov_cc[[CONTRAST_COL] + COVARS].astype(float))
    contrast_idx = list(X.columns).index(CONTRAST_COL)
    XtX_inv = np.linalg.pinv(X.values.T @ X.values)
    se_ord_factor = float(XtX_inv[contrast_idx, contrast_idx])

    rows = []
    s2_all, beta_all = [], []
    df_resid_all = None
    for feat_idx in range(n_feat):
        y = matrix_cc.iloc[feat_idx].values.astype(float)
        meta = df.loc[feat_idx, META_COLS].to_dict()

        beta, se_hc3, p_hc3, ci_lo, ci_hi, s2_ols, df_resid, model_ord = fit_per_feature(
            y, X, contrast_idx)
        df_resid_all = df_resid

        g = cov_cc[CONTRAST_COL].values
        mean_normal = float(np.mean(y[g == 0]))
        mean_over = float(np.mean(y[g == 1]))

        try:
            _, p_wil = stats.mannwhitneyu(y[g == 1], y[g == 0], alternative='two-sided')
            p_wil = float(p_wil)
        except ValueError:
            p_wil = np.nan

        out_mask = detect_outliers_studentized(model_ord)
        n_out = int(out_mask.sum())
        robust_info = refit_after_outlier(y, X, ~out_mask, contrast_idx, beta)
        robust_flag = is_robust_to_outlier(n_out, robust_info)

        rows.append({
            **meta,
            'n': n_cc,
            'mean_log2_normal': round(mean_normal, 4),
            'mean_log2_overweight': round(mean_over, 4),
            'log2FC': round(beta, 4),
            'log2FC_CI_lo': round(ci_lo, 4),
            'log2FC_CI_hi': round(ci_hi, 4),
            'p_ols_hc3': p_hc3,
            'p_limma': np.nan,
            'p_wilcoxon': p_wil,
            'n_outliers': n_out,
            'outlier_removed': bool(n_out > 0),
            'beta_after': round(robust_info['beta_after'], 4) if pd.notna(robust_info['beta_after']) else np.nan,
            'direction_consistent': robust_info['direction_consistent'],
            'attenuation_pct': round(robust_info['attenuation_pct'], 4) if pd.notna(robust_info['attenuation_pct']) else np.nan,
            'p_robust_ols': round(robust_info['p_robust_ols'], 6) if pd.notna(robust_info['p_robust_ols']) else np.nan,
            'is_robust_to_outlier': robust_flag,
        })
        s2_all.append(s2_ols)
        beta_all.append(beta)

    s2_arr = np.array(s2_all)
    beta_arr = np.array(beta_all)
    mod = moderated_t(beta_arr, se_ord_factor, s2_arr, df=df_resid_all)
    print(f'  limma 先验: d_0={mod["df_prior"]:.2f}, s2_0={mod["s2_prior"]:.4f}, '
          f'df_total={float(np.atleast_1d(mod["df_total"])[0]):.2f}')

    for i, r in enumerate(rows):
        r['p_limma'] = float(mod['p_mod'][i])
        r['t_limma'] = round(float(mod['t_mod'][i]), 4)

    tbl = pd.DataFrame(rows)
    tbl['q_ols_hc3_bh'] = bh_qvalues(tbl['p_ols_hc3'].values)
    tbl['q_limma_bh'] = bh_qvalues(tbl['p_limma'].values)
    tbl['q_wilcoxon_bh'] = bh_qvalues(tbl['p_wilcoxon'].values)

    tbl['tier'] = [
        assign_tier(ql, pl, fc)
        for ql, pl, fc in zip(tbl['q_limma_bh'], tbl['p_limma'], tbl['log2FC'])
    ]
    tbl['high_confidence'] = [
        is_high_confidence(t, p, lo, hi)
        for t, p, lo, hi in zip(tbl['tier'], tbl['p_ols_hc3'],
                                tbl['log2FC_CI_lo'], tbl['log2FC_CI_hi'])
    ]
    tbl['direction'] = np.where(tbl['log2FC'] > 0, '↑(超重肥胖>正常)', '↓(超重肥胖<正常)')

    if tag == '50' and main80_features is not None:
        tbl['also_in_80_main'] = tbl['Metabolite Name'].isin(main80_features)

    for c in ['p_ols_hc3', 'p_limma', 'p_wilcoxon',
              'q_ols_hc3_bh', 'q_limma_bh', 'q_wilcoxon_bh']:
        tbl[c] = tbl[c].apply(lambda x: round(float(x), 6) if pd.notna(x) else np.nan)

    tier_rank = {'strict': 0, 'standard': 1, 'exploratory': 2, '-': 3}
    tbl = (tbl.assign(_tr=tbl['tier'].map(tier_rank))
              .sort_values(['_tr', 'q_limma_bh', 'p_limma'])
              .drop(columns='_tr')
              .reset_index(drop=True))

    front = ['Metabolite Name', 'Chinese Name', 'Class', 'Sub Class',
             'Super Class', 'Direct Parent', 'KEGG ID', 'HMDB ID',
             'Retention Time(min)', 'Cellular Locations',
             'n', 'tier', 'high_confidence', 'direction',
             'log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi',
             'mean_log2_normal', 'mean_log2_overweight',
             'p_ols_hc3', 'q_ols_hc3_bh',
             'p_limma', 't_limma', 'q_limma_bh',
             'p_wilcoxon', 'q_wilcoxon_bh',
             'n_outliers', 'outlier_removed', 'beta_after',
             'direction_consistent', 'attenuation_pct',
             'p_robust_ols', 'is_robust_to_outlier']
    if 'also_in_80_main' in tbl.columns:
        front.append('also_in_80_main')
    tbl = tbl[front]

    out_file = OUT_DIR / f'ancova_main_{tag}.csv'
    atomic_write_csv(tbl, out_file)
    print(f'  已写入: {out_file.relative_to(ROOT)}')

    s2_min, s2_med, s2_max = float(s2_arr.min()), float(np.median(s2_arr)), float(s2_arr.max())
    s2_var_ratio = s2_max / s2_min if s2_min > 0 else np.inf

    n_strict = int((tbl['tier'] == 'strict').sum())
    n_standard = int((tbl['tier'] == 'standard').sum())
    n_exp = int((tbl['tier'] == 'exploratory').sum())
    n_robust = int(tbl['is_robust_to_outlier'].sum())
    n_high_conf = int(tbl['high_confidence'].sum())
    n_with_outlier = int(tbl['outlier_removed'].sum())
    print(f'  入选: strict={n_strict}, standard={n_standard}, exploratory={n_exp}')
    print(f'  high_confidence (limma+OLS+CI 三认同): {n_high_conf}')
    print(f'  含离群的特征: {n_with_outlier} / {n_feat} (中位 n_out={int(tbl["n_outliers"].median())})')
    print(f'  离群稳健 (方向一致+衰减≤50%+p<0.05): {n_robust} / {n_feat}')
    print(f'  特征方差: s²∈[{s2_min:.3f}, {s2_max:.3f}], median {s2_med:.3f}, max/min {s2_var_ratio:.1f}')

    d0_str = (f'{float(mod["df_prior"]):.2f}' if np.isfinite(mod['df_prior'])
              and float(mod['df_prior']) < 1e6 else '≈∞ (特征方差近似一致, EB 退化为池化方差)')
    return {
        'tier_track': tag,
        'src': src_name,
        'n_feat': n_feat,
        'n_sample': n_cc,
        'n_normal': n_normal,
        'n_overweight': n_over,
        'df_resid': df_resid_all,
        'limma_d0': d0_str,
        'limma_s2_0': round(float(mod['s2_prior']), 4),
        's2_min': round(s2_min, 4),
        's2_median': round(s2_med, 4),
        's2_max': round(s2_max, 4),
        's2_max_min_ratio': round(s2_var_ratio, 2),
        'n_strict_qlimma05': n_strict,
        'n_standard_qlimma10_fc1.2': n_standard,
        'n_exploratory_plimma05_fc1.5': n_exp,
        'n_high_confidence': n_high_conf,
        'n_features_with_outlier': n_with_outlier,
        'n_outlier_robust': n_robust,
        'out_file': out_file.name,
    }, set(tbl['Metabolite Name'].tolist())


def main():
    print('=== 09. ANCOVA / limma 主分析 (双轨, 方案 D + 协变量二轮锁定版) ===')
    print('模型: log2(metab) ~ BMI_group + age + GA_decimal + year (B20 ref)')

    cov_full = load_covariates()
    print(f'  协变量表: {len(cov_full)} 样本, '
          f'BMI={cov_full["BMI_group"].value_counts().to_dict()}')

    records = []
    rec80, feat80 = process_track('ori_n165_filtered80_log2_combat.xlsx', '80', cov_full)
    records.append(rec80)
    rec50, _ = process_track('ori_n165_filtered50_log2_combat.xlsx', '50', cov_full,
                              main80_features=feat80)
    records.append(rec50)

    audit = pd.DataFrame(records)
    audit_path = OUT_DIR / 'ancova_audit.csv'
    atomic_write_csv(audit, audit_path)
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
