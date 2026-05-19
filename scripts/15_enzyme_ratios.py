"""15. 酶活性代理比值 ANCOVA — 按 COX/LOX/CYP/sEH 四酶严格组织 (v2 重写, 2026-05-19).

11 个比值:
  COX  (2):  total PG_AA / AA                      — AA 系 COX 总通量
             PGE1 / PGE2                            — DGLA vs AA COX 偏移
  LOX  (3):  12-HETE / AA                          — 纯 12-LOX 通量 (AA)
             15-HETE / AA                          — 纯 15-LOX 通量 (AA)
             12-HEPE / EPA                          — 12-LOX 通量 (EPA)
  CYP  (1):  20-HDoHE / DHA                        — CYP4 ω-3 羟化通量 (DHA, 严格版)
                                                     (2026-05-19 修订: 16-HDoHE 因非酶主导剔除, 见 12 号脚本)
  sEH  (5):  8,9-DiHETrE / 8,9-EpETrE              — sEH 8,9 位 (AA)
             11,12-DiHETrE / 11,12-EpETrE          — sEH 11,12 位 (AA)
             9,10-DiHOME / 9,10-EpOME              — sEH 9,10 位 (LA)
             12,13-DiHOME / 12,13-EpOME            — sEH 12,13 位 (LA)
             19,20-DiHDPA / 19,20-EpDPA            — sEH 19,20 位 (DPA)

关键设计:
  - 用 raw 浓度 (QRILC 插补后, pre-log2, pre-ComBat) 算比值; 比值层面 batch effect 自然消减
    (同一样本同一批次的相关分子受相同乘性偏移, 比值代数相消)
  - ANCOVA 加 batch 哑变量作二重保护 (替代 ComBat 校正)
  - 通路严格: 剔除 11/8-HETE (Auto-ox 混杂), 7/10/11/13-HDoHE (非酶氧化主导)
  - sEH 按位置异构体配对 (不做"全 DHET / 部分 EET"的错配比值)

数据来源:
  67 主轨 metabolite : data/02_preprocessed/ori_n165_filtered80_imputed.xlsx (QRILC 插补 raw)
  EPA + DHA 母体     : data/01_raw/ori_allmet_QC.xlsx (panel 缺, 从 673-met 全表抽取;
                                                       100% 检出, QC CV%=5.5%/17.4%)
  上机批次          : data/01_raw/上机顺序和批次表.xlsx (C18 sheet, batch=1/2)
  临床协变量        : data/01_raw/分娩特征分析_2019-2021_主分析队列n165.xlsx

模型: log2(ratio) ~ BMI_overweight + age + GA_decimal + year_A19 + year_C21 + batch_2
  n=165, 7 参数, df_resid = 158
  OLS + HC3 稳健 SE; Mann-Whitney U 互校; BH (11 比值) 校正

输出:
  results/tables/enzyme_ratios_ancova.csv       (11 比值统计 + enzyme 列)
  results/tables/enzyme_ratios_matrix.csv       (11 × 165, log2 比值)
  results/figures/04_Metabolite_boxplots/enzyme_ratios_grid.png  (3×4 网格, 11 比值 + 1 空)
"""
import os
import sys
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
PREP_DIR = ROOT / 'data' / '02_preprocessed'
RAW_DIR = ROOT / 'data' / '01_raw'
TABLES = ROOT / 'results' / 'tables'
FIG_DIR = ROOT / 'results' / 'figures' / '04_Metabolite_boxplots'
FIG_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'
BMI_COLORS = {'正常': '#67A9CF', '超重肥胖': '#EF8A62'}
ENZYME_COLORS = {'COX': '#D62728', 'LOX': '#2CA02C', 'CYP': '#1F77B4', 'sEH': '#9467BD'}

# Exact metabolite names per panel (ori_n165_filtered80_imputed.xlsx + ori_allmet_QC.xlsx)
EPA = '5Z,8Z,11Z,14Z,17Z-Eicosapentaenoic acid'        # only in ori_allmet_QC
DHA = '4Z,7Z,10Z,13Z,16Z,19Z-Docosahexaenoic acid'     # only in ori_allmet_QC

# (enzyme_group, ratio_name, rationale, numerator metab list, denominator metab list)
RATIOS = [
    # ===== COX =====
    ('COX', 'total_PG_AA / AA', 'AA 系 COX 总通量',
     ['Prostaglandin D2', 'Prostaglandin E2', 'Prostaglandin F2α',
      'Prostaglandin J2', 'Prostaglandin A2', '11β-Prostaglandin E2',
      'Thromboxane B2', '12S-Hydroxy-5Z,8E,10E-heptadecatrienoic acid'],
     ['Arachidonic acid']),

    ('COX', 'PGE1 / PGE2', 'DGLA vs AA COX 偏移 (PGE1 抗炎 vs PGE2 促炎)',
     ['Prostaglandin E1'], ['Prostaglandin E2']),

    # ===== LOX =====
    ('LOX', '12-HETE / AA', '纯 12-LOX 通量 (AA 系)',
     ['12-Hydroxy-5Z,8Z,10E,14Z-eicosatetraenoic acid'],
     ['Arachidonic acid']),

    ('LOX', '15-HETE / AA', '纯 15-LOX 通量 (AA 系)',
     ['15-Hydroxy-5Z,8Z,11Z,13E-eicosatetraenoic acid'],
     ['Arachidonic acid']),

    ('LOX', '12-HEPE / EPA', '12-LOX 通量 (EPA 系)',
     ['12-Hydroxy-5,8,10,14,17-eicosapentaenoic acid'],
     [EPA]),

    # ===== CYP =====
    # 2026-05-19 修订: 原 (16+20-HDoHE)/DHA 混合 CYP+非酶, 改为纯 20-HDoHE/DHA.
    # 依据: 16-HDoHE 实证 log2FC+0.173 与非酶组同质 (VanRollins 2008: CYP4 偏好 ω-3 即 C20,
    # ω-7 即 C16 不是 CYP4 偏好位; Yin & Porter 2005: 16-HDHA 17E 命名契合自由基机制).
    # 剔除 16-HDoHE 后 CYP 信号反而增强 (FC 1.31× → 1.35×, p 0.009 → 0.008).
    ('CYP', '20-HDoHE / DHA', 'CYP4 ω-3 羟化通量 (DHA 系, 严格版)',
     ['20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid'],
     [DHA]),

    # ===== sEH (按位置异构体配对) =====
    ('sEH', '8,9-DiHETrE / 8,9-EpETrE', 'sEH 8,9 位 (AA 通路)',
     ['8,9-DiHydroxy-5Z,11Z,14Z-eicosatrienoic acid'],
     ['8,9-Dpoxy-5Z,11Z,14Z-eicosatrienoic acid']),

    ('sEH', '11,12-DiHETrE / 11,12-EpETrE', 'sEH 11,12 位 (AA 通路)',
     ['11,12-DiHydroxy-5Z,8Z,14Z-eicosatrienoic acid'],
     ['11,12-Epoxy-5Z,8Z,14Z-eicosatrienoic acid']),

    ('sEH', '9,10-DiHOME / 9,10-EpOME', 'sEH 9,10 位 (LA 通路)',
     ['9,10-DiHydroxy-12Z-octadecenoic acid'],
     ['9,10-Epoxy-12Z-octadecenoic acid']),

    ('sEH', '12,13-DiHOME / 12,13-EpOME', 'sEH 12,13 位 (LA 通路)',
     ['12,13 -DiHydroxy-9Z-octadecenoic acid'],
     ['12,13-Epoxy-9Z-octadecenoic acid']),

    ('sEH', '19,20-DiHDPA / 19,20-EpDPA', 'sEH 19,20 位 (DPA 通路)',
     ['19,20-DiHydroxy-4Z,7Z,10Z,13Z,16Z-docosapentaenoic acid'],
     ['19(20)-Epoxy-4Z,7Z,10Z,13Z,16Z-docosapentaenoic acid']),
]


def load_raw_matrix():
    """Load 67 主轨 (QRILC 插补) + EPA/DHA (ori_allmet_QC raw), align to 165 omx_ids."""
    df80 = pd.read_excel(PREP_DIR / 'ori_n165_filtered80_imputed.xlsx')
    sample_cols = [c for c in df80.columns if c not in META_COLS and c != UNIT_COL]
    mat80 = df80.set_index('Metabolite Name')[sample_cols]
    print(f'  80 主轨 raw 矩阵: {mat80.shape[0]} 特征 × {mat80.shape[1]} 样本')

    df_all = pd.read_excel(RAW_DIR / 'ori_allmet_QC.xlsx')
    all_sample_cols = [c for c in df_all.columns if c not in META_COLS and c != UNIT_COL]
    epa_dha = df_all[df_all['Metabolite Name'].isin([EPA, DHA])].set_index(
        'Metabolite Name')[all_sample_cols]
    # Align EPA/DHA columns to the 165 main-cohort omx_ids
    missing_in_qc = [s for s in sample_cols if s not in all_sample_cols]
    if missing_in_qc:
        raise RuntimeError(f'  {len(missing_in_qc)} 主队列样本不在 QC 表内: {missing_in_qc[:5]}')
    epa_dha_165 = epa_dha[sample_cols]
    print(f'  EPA/DHA (from 673-met 全表): {epa_dha_165.shape[0]} 特征 (检出 EPA/DHA 应 100%)')

    # Detection diagnostics for EPA/DHA on 165
    for nm in [EPA, DHA]:
        vals = pd.to_numeric(epa_dha_165.loc[nm], errors='coerce')
        nz = vals[(vals > 0) & vals.notna()]
        print(f'    [{nm}] 检出 {len(nz)}/{len(sample_cols)} '
              f'({len(nz)/len(sample_cols)*100:.1f}%) 中位 {nz.median():.3g} ng/g')

    full = pd.concat([mat80, epa_dha_165])
    print(f'  合并矩阵: {full.shape[0]} 特征 × {full.shape[1]} 样本')
    return full, sample_cols


def load_batch_table():
    """Read C18 sheet, return omx_id -> batch (1/2)."""
    df = pd.read_excel(RAW_DIR / '上机顺序和批次表.xlsx',
                       sheet_name='C18柱-上机顺序和批次表', header=1, usecols=range(4))
    df.columns = ['order', 'class', 'omx_id', 'batch']
    df['batch'] = pd.to_numeric(df['batch'], errors='coerce').astype('Int64')
    return df[df['class'] == 'Subject'][['omx_id', 'batch']].reset_index(drop=True)


def build_covariates(sample_cols):
    align = pd.read_csv(PREP_DIR / 'sample_alignment_n165.csv')
    clin = pd.read_excel(RAW_DIR / '分娩特征分析_2019-2021_主分析队列n165.xlsx')
    batch_tbl = load_batch_table()
    cov = (clin[['标本号', '孕妇年龄', 'GA_decimal', 'BMI_group']]
           .rename(columns={'孕妇年龄': 'age'})
           .merge(align[['标本号', 'omx_id']], on='标本号', how='left')
           .merge(batch_tbl, on='omx_id', how='left'))
    cov['BMI_overweight'] = (cov['BMI_group'] == '超重肥胖').astype(float)
    yp = cov['omx_id'].str[:3]
    cov['year_A19'] = (yp == 'A19').astype(float)
    cov['year_C21'] = (yp == 'C21').astype(float)
    cov['batch_2'] = (cov['batch'].astype(float) == 2.0).astype(float)
    miss = cov[['batch', 'age', 'GA_decimal']].isna().sum()
    if miss.any():
        raise RuntimeError(f'协变量有缺失: {miss.to_dict()}')
    return cov.set_index('omx_id').loc[sample_cols]


def fit_ancova(y, X, contrast_idx):
    """OLS+HC3 → (beta, se, p, ci_lo, ci_hi)."""
    m = sm.OLS(y, X).fit(cov_type='HC3')
    beta = float(m.params.iloc[contrast_idx])
    se = float(m.bse.iloc[contrast_idx])
    p = float(m.pvalues.iloc[contrast_idx])
    ci_lo, ci_hi = m.conf_int().iloc[contrast_idx].tolist()
    return beta, se, p, float(ci_lo), float(ci_hi)


def compute_ratio_log2(matrix_raw, num_list, denom_list, label):
    """Sum numerator metabs on raw scale / sum denominator metabs on raw scale, then log2."""
    num_missing = [m for m in num_list if m not in matrix_raw.index]
    den_missing = [m for m in denom_list if m not in matrix_raw.index]
    if num_missing or den_missing:
        raise RuntimeError(f'  [{label}] 缺分子 {num_missing} | 缺分母 {den_missing}')
    num_sum = matrix_raw.loc[num_list].sum(axis=0).astype(float).values
    den_sum = matrix_raw.loc[denom_list].sum(axis=0).astype(float).values
    if (num_sum <= 0).any() or (den_sum <= 0).any():
        raise RuntimeError(f'  [{label}] raw 求和出现 <=0 (插补失败?)')
    return np.log2(num_sum) - np.log2(den_sum)


def atomic_to_csv(df, path):
    fd, tmp = tempfile.mkstemp(suffix='.csv', dir=str(path.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, index=False, encoding='utf-8-sig')
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def atomic_to_csv_with_index(df, path):
    fd, tmp = tempfile.mkstemp(suffix='.csv', dir=str(path.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, encoding='utf-8-sig')
        os.replace(tmp, path)
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def draw_ratio_boxplot(ax, log2_ratio, groups, palette, title_lines, ylabel,
                        enzyme=None):
    grp_order = ['正常', '超重肥胖']
    data = [log2_ratio[groups == g] for g in grp_order]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.55, patch_artist=True,
                    boxprops=dict(linewidth=1.0),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=1.0),
                    flierprops=dict(marker='', markersize=0))
    for b, g in zip(bp['boxes'], grp_order):
        b.set_facecolor(palette[g])
        b.set_alpha(0.5)
        b.set_edgecolor(palette[g])
    rng = np.random.default_rng(42)
    for x_pos, g in enumerate(grp_order):
        yvals = log2_ratio[groups == g]
        x_jit = rng.normal(x_pos, 0.07, size=len(yvals))
        ax.scatter(x_jit, yvals, s=14, c=palette[g], edgecolors='black',
                   linewidths=0.3, alpha=0.55, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'{g}\n(n={int((groups == g).sum())})' for g in grp_order],
                       fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    title_color = ENZYME_COLORS.get(enzyme, 'black')
    ax.set_title('\n'.join(title_lines), fontsize=9.5, pad=6, color=title_color)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    # Color the spine top to indicate enzyme group
    if enzyme:
        ax.spines['top'].set_color(title_color)
        ax.spines['top'].set_linewidth(2.5)


def main():
    print('=== 15. 酶活性比值 ANCOVA v2 (按 COX/LOX/CYP/sEH 严格组织) ===\n')

    matrix_raw, sample_cols = load_raw_matrix()
    cov_idx = build_covariates(sample_cols)
    groups = cov_idx['BMI_group'].values
    print(f'\n  协变量队列: n={len(cov_idx)}; '
          f'BMI {pd.Series(groups).value_counts().to_dict()}; '
          f"batch {cov_idx['batch'].value_counts().sort_index().to_dict()}; "
          f"year_A19={int(cov_idx['year_A19'].sum())}, "
          f"year_C21={int(cov_idx['year_C21'].sum())}")

    X = sm.add_constant(cov_idx[['BMI_overweight', 'age', 'GA_decimal',
                                  'year_A19', 'year_C21', 'batch_2']].astype(float))
    contrast_idx = list(X.columns).index('BMI_overweight')
    print(f'  设计矩阵: n={X.shape[0]}, 参数 p={X.shape[1]}, df_resid={X.shape[0] - X.shape[1]}\n')

    ratio_matrix = pd.DataFrame(index=sample_cols)
    results = []
    print(f'  {"#":<3} {"enzyme":<5} {"ratio":<35} {"log2FC":>8} {"p_OLS":>7} '
          f'{"p_Wilc":>7} {"dir":<3}')
    print('  ' + '-' * 80)
    for i, (enzyme, name, rationale, num_list, denom_list) in enumerate(RATIOS, 1):
        log2_ratio = compute_ratio_log2(matrix_raw, num_list, denom_list, name)
        ratio_matrix[name] = log2_ratio

        beta, se, p, ci_lo, ci_hi = fit_ancova(log2_ratio, X, contrast_idx)
        try:
            _, p_wil = stats.mannwhitneyu(
                log2_ratio[groups == '超重肥胖'],
                log2_ratio[groups == '正常'], alternative='two-sided')
            p_wil = float(p_wil)
        except ValueError:
            p_wil = np.nan
        mean_norm = float(np.mean(log2_ratio[groups == '正常']))
        mean_over = float(np.mean(log2_ratio[groups == '超重肥胖']))
        direction = '↑' if beta > 0 else '↓'

        results.append({
            'enzyme': enzyme,
            'ratio': name,
            'rationale': rationale,
            'n_numerator': len(num_list),
            'n_denominator': len(denom_list),
            'mean_log2_normal': round(mean_norm, 4),
            'mean_log2_overweight': round(mean_over, 4),
            'log2FC': round(beta, 4),
            'log2FC_CI_lo': round(ci_lo, 4),
            'log2FC_CI_hi': round(ci_hi, 4),
            'fold_change_linear': round(2.0 ** beta, 4),
            'direction': direction,
            'p_ols_hc3': round(p, 6),
            'p_wilcoxon': round(p_wil, 6) if pd.notna(p_wil) else np.nan,
        })
        print(f'  {i:<3} {enzyme:<5} {name:<35} {beta:+.3f}   {p:.4f}  '
              f'{p_wil:.4f}   {direction}')

    res = pd.DataFrame(results)
    res['q_ols_hc3_bh'] = multipletests(res['p_ols_hc3'], method='fdr_bh')[1].round(6)
    res['q_wilcox_bh'] = multipletests(res['p_wilcoxon'], method='fdr_bh')[1].round(6)

    # Column order
    res = res[['enzyme', 'ratio', 'rationale',
               'n_numerator', 'n_denominator',
               'mean_log2_normal', 'mean_log2_overweight',
               'log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi', 'fold_change_linear',
               'direction',
               'p_ols_hc3', 'q_ols_hc3_bh',
               'p_wilcoxon', 'q_wilcox_bh']]

    out_csv = TABLES / 'enzyme_ratios_ancova.csv'
    atomic_to_csv(res, out_csv)
    print(f'\n  ✓ {out_csv.relative_to(ROOT)}')

    mat_csv = TABLES / 'enzyme_ratios_matrix.csv'
    atomic_to_csv_with_index(ratio_matrix, mat_csv)
    print(f'  ✓ {mat_csv.relative_to(ROOT)}')

    # Plot: 3 rows × 4 cols, 11 ratios + 1 blank
    fig, axes = plt.subplots(3, 4, figsize=(16, 11))
    axes_flat = axes.flatten()
    for i, r in res.iterrows():
        ax = axes_flat[i]
        log2_ratio = ratio_matrix[r['ratio']].values
        ann_q = f"q_BH={r['q_ols_hc3_bh']:.3f}" if pd.notna(r['q_ols_hc3_bh']) else ''
        title = [
            f"[{r['enzyme']}] {r['ratio']}",
            f"log2FC={r['log2FC']:+.3f}  FC={r['fold_change_linear']:.2f}×  "
            f"p={r['p_ols_hc3']:.4f}  {ann_q}",
            r['rationale'],
        ]
        draw_ratio_boxplot(ax, log2_ratio, groups, BMI_COLORS, title,
                            ylabel='log2 (ratio, raw)', enzyme=r['enzyme'])
    for j in range(len(res), len(axes_flat)):
        axes_flat[j].axis('off')

    # Top legend strip for enzyme colors
    handles = [matplotlib.patches.Patch(facecolor=c, edgecolor='black', label=e)
               for e, c in ENZYME_COLORS.items()]
    fig.legend(handles=handles, loc='upper center', ncol=4, fontsize=11,
               frameon=False, bbox_to_anchor=(0.5, 0.985))
    fig.suptitle('Enzyme activity proxy ratios — COX / LOX / CYP / sEH (n=165, '
                 'ANCOVA-adjusted for age, GA, year, batch)',
                 fontsize=13, fontweight='bold', y=0.965)
    plt.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.05,
                        hspace=0.65, wspace=0.32)
    fig_path = FIG_DIR / 'enzyme_ratios_grid.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {fig_path.relative_to(ROOT)}')

    # Print summary by enzyme group
    print('\n  === Summary by enzyme ===')
    for enz in ['COX', 'LOX', 'CYP', 'sEH']:
        sub = res[res['enzyme'] == enz]
        n_up = (sub['direction'] == '↑').sum()
        n_dn = (sub['direction'] == '↓').sum()
        n_sig_raw = (sub['p_ols_hc3'] < 0.05).sum()
        n_sig_q = (sub['q_ols_hc3_bh'] < 0.05).sum()
        print(f'    {enz:<4} n={len(sub)} | ↑{n_up} ↓{n_dn} | '
              f'raw p<0.05: {n_sig_raw} | BH q<0.05: {n_sig_q}')


if __name__ == '__main__':
    main()
