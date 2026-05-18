"""15. 酶活性代理比值 ANCOVA (SOP §5.2).

设计 6 个比值, 反映核心氧化脂质代谢酶通量:

  1. HETE_total / AA              → LOX 通量 (AA → HETE; 包含 8/11/12/15-HETE)
  2. DiHETrE_total / EpETrE_total → sEH 活性 (AA 通路, **核心指标**)
  3. DiHOME_total / EpOME_total   → sEH 活性 (LA 通路)
  4. PGE2 / PGD2                  → COX 下游促炎 vs 消退倾斜
  5. HODE_total / LA              → LA-LOX 活性
  6. HDoHE_omega / HDoHE_LOX      → DHA 通路 CYP-ω vs LOX 偏移

数据处理:
  - 输入: ComBat 校正后 log2 矩阵 (80 主轨 67 特征; 50 轨可补)
  - 反 log2 → 同家族求和 (raw scale) → log2 → 减法 (= log2 比值)
  - ANCOVA: log2(ratio) ~ BMI + age + GA_decimal + year (B20 ref, year_A19 + year_C21)
  - 输出: 比值矩阵 + ANCOVA 统计 + 6 个 boxplot 合一张图

输入:
  data/03_batch_corrected/ori_n165_filtered80_log2_combat.xlsx (主)
  data/02_preprocessed/sample_alignment_n165.csv
  data/01_raw/分娩特征分析_2019-2021_主分析队列n165.xlsx

输出:
  results/tables/enzyme_ratios_ancova.csv
  results/tables/enzyme_ratios_matrix.csv (6 比值 × 165 样本, log2)
  results/figures/04_Metabolite_boxplots/enzyme_ratios_grid.png
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
COMBAT_DIR = ROOT / 'data' / '03_batch_corrected'
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

# (比值名, 分子代谢物列表, 分母代谢物列表, 反映)
RATIOS = [
    ('HETE/AA',
     ['11-Hydroxy-5Z,8Z,11E,14Z-eicosatetraenoic acid',
      '15-Hydroxy-5Z,8Z,11Z,13E-eicosatetraenoic acid',
      '8-Hydroxy-5Z,9E,11Z,14Z-eicosatetraenoic acid',
      '12-Hydroxy-5Z,8Z,10E,14Z-eicosatetraenoic acid'],
     ['Arachidonic acid'],
     'LOX 通量 (AA → HETE)'),

    ('DiHETrE/EpETrE',
     ['5,6-DiHydroxy-8Z,11Z,14Z-eicosatrienoic acid',
      '8,9-DiHydroxy-5Z,11Z,14Z-eicosatrienoic acid',
      '11,12-DiHydroxy-5Z,8Z,14Z-eicosatrienoic acid',
      '14,15-DiHydroxy-5Z,8Z,11Z-eicosatrienoic acid'],
     ['8,9-Dpoxy-5Z,11Z,14Z-eicosatrienoic acid',
      '11,12-Epoxy-5Z,8Z,14Z-eicosatrienoic acid'],
     'sEH 活性 (AA 通路, 核心)'),

    ('DiHOME/EpOME',
     ['9,10-DiHydroxy-12Z-octadecenoic acid',
      '12,13 -DiHydroxy-9Z-octadecenoic acid'],
     ['9,10-Epoxy-12Z-octadecenoic acid',
      '12,13-Epoxy-9Z-octadecenoic acid'],
     'sEH 活性 (LA 通路)'),

    ('PGE2/PGD2',
     ['Prostaglandin E2'],
     ['Prostaglandin D2'],
     'COX 下游促炎 vs 消退'),

    ('HODE/LA',
     ['9-Hydroxy-10E,12Z-octadecadienoic acid',
      '13-Hydroxy-9Z,11E-octadecadienoic acid'],
     ['Linoleic acid'],
     'LA-LOX 活性'),

    ('HDoHE_omega/HDoHE_LOX',
     ['16-Hydroxy-4Z,7Z,10Z,13Z,17E,19Z-docosahexaenoic acid',
      '20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid'],
     ['14-Hydroxy-4Z,7Z,10Z,12E,16Z,19Z-docosahexaenoic acid',
      '10-Hydroxy-4Z,7Z,11E,13Z,16Z,19Z-docosahexaenoic acid',
      '11-Hydroxy-4Z,7Z,9E,13Z,16Z,19Z-docosahexaenoic acid',
      '13-Hydroxy-4Z,7Z,10Z,14E,16Z,19Z-docosahexaenoic acid',
      '7-Hydroxy-4Z,8E,10Z,13Z,16Z,19Z-docosahexaenoic acid'],
     'DHA 通路 CYP-ω vs LOX'),
]


def sum_in_raw(matrix_log2, metab_names):
    """log2 矩阵中给定 metabolite 子集 → 反 log2 → 按样本求和 → 返回 raw 浓度向量."""
    avail = [m for m in metab_names if m in matrix_log2.index]
    if not avail:
        return None, []
    sub = matrix_log2.loc[avail]
    raw = np.power(2.0, sub.values)  # log2 → raw
    summed = raw.sum(axis=0)
    return summed, avail


def fit_ancova(y, X, contrast_idx):
    """OLS+HC3 → (β, se, p, ci_lo, ci_hi)"""
    m = sm.OLS(y, X).fit(cov_type='HC3')
    beta = float(m.params.iloc[contrast_idx])
    se = float(m.bse.iloc[contrast_idx])
    p = float(m.pvalues.iloc[contrast_idx])
    ci_lo, ci_hi = m.conf_int().iloc[contrast_idx].tolist()
    return beta, se, p, float(ci_lo), float(ci_hi)


def draw_ratio_boxplot(ax, log2_ratio, groups, palette, title_lines, ylabel):
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
    ax.set_xticklabels([f'{g}\n(n={int((groups == g).sum())})' for g in grp_order], fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title('\n'.join(title_lines), fontsize=10, pad=6)
    ax.grid(axis='y', alpha=0.3, linestyle='--')


def main():
    print('=== 15. 酶活性比值 ANCOVA ===\n')

    df = pd.read_excel(COMBAT_DIR / 'ori_n165_filtered80_log2_combat.xlsx')
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    matrix = df.set_index('Metabolite Name')[sample_cols]
    n_samp = len(sample_cols)
    print(f'  输入: 80 主轨 ({len(matrix)} 特征 × {n_samp} 样本)')

    # 协变量
    align = pd.read_csv(PREP_DIR / 'sample_alignment_n165.csv')
    clin = pd.read_excel(RAW_DIR / '分娩特征分析_2019-2021_主分析队列n165.xlsx')
    cov = (clin[['标本号', '孕妇年龄', 'GA_decimal', 'BMI_group']]
           .rename(columns={'孕妇年龄': 'age'})
           .merge(align[['标本号', 'omx_id']], on='标本号', how='left'))
    cov['BMI_overweight'] = (cov['BMI_group'] == '超重肥胖').astype(float)
    yp = cov['omx_id'].str[:3]
    cov['year_A19'] = (yp == 'A19').astype(float)
    cov['year_C21'] = (yp == 'C21').astype(float)
    cov_idx = cov.set_index('omx_id').loc[sample_cols]
    groups = cov_idx['BMI_group'].values

    X = sm.add_constant(cov_idx[['BMI_overweight', 'age', 'GA_decimal',
                                  'year_A19', 'year_C21']].astype(float))
    contrast_idx = list(X.columns).index('BMI_overweight')

    # 计算比值矩阵 + ANCOVA
    ratio_matrix = pd.DataFrame(index=sample_cols)
    results = []
    for ratio_name, num_list, denom_list, rationale in RATIOS:
        num_sum, num_avail = sum_in_raw(matrix, num_list)
        denom_sum, denom_avail = sum_in_raw(matrix, denom_list)
        if num_sum is None or denom_sum is None:
            print(f'  ⚠ {ratio_name}: 跳过 (分子或分母在 80 轨内无可用代谢物)')
            continue
        log2_ratio = np.log2(num_sum) - np.log2(denom_sum)
        ratio_matrix[ratio_name] = log2_ratio

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
        direction = '↑(超重肥胖>正常)' if beta > 0 else '↓(超重肥胖<正常)'

        results.append({
            'ratio': ratio_name,
            'rationale': rationale,
            'n_numerator_metabs': len(num_avail),
            'n_denominator_metabs': len(denom_avail),
            'mean_log2_normal': round(mean_norm, 4),
            'mean_log2_overweight': round(mean_over, 4),
            'log2FC': round(beta, 4),
            'log2FC_CI_lo': round(ci_lo, 4),
            'log2FC_CI_hi': round(ci_hi, 4),
            'direction': direction,
            'p_ols_hc3': round(p, 6),
            'p_wilcoxon': round(p_wil, 6) if pd.notna(p_wil) else np.nan,
        })
        print(f'  {ratio_name:30s}  log2FC={beta:+.3f}  p_OLS={p:.4f}  p_Wilcox={p_wil:.4f}  {direction}')

    res = pd.DataFrame(results)
    # BH 校正
    from statsmodels.stats.multitest import multipletests
    res['q_ols_hc3_bh'] = multipletests(res['p_ols_hc3'], method='fdr_bh')[1]
    res['q_ols_hc3_bh'] = res['q_ols_hc3_bh'].round(6)

    out_csv = TABLES / 'enzyme_ratios_ancova.csv'
    res.to_csv(out_csv, index=False, encoding='utf-8-sig')
    print(f'\n  ✓ {out_csv.relative_to(ROOT)}')

    mat_csv = TABLES / 'enzyme_ratios_matrix.csv'
    ratio_matrix.to_csv(mat_csv, encoding='utf-8-sig')
    print(f'  ✓ {mat_csv.relative_to(ROOT)}')

    # 6 比值合并 boxplot (2×3)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for i, r in res.iterrows():
        ax = axes[i // 3, i % 3]
        log2_ratio = ratio_matrix[r['ratio']].values
        title = [
            f'{r["ratio"]}',
            f'log2FC={r["log2FC"]:+.3f}  p={r["p_ols_hc3"]:.4f}  q={r["q_ols_hc3_bh"]:.3f}',
            r['rationale'],
        ]
        draw_ratio_boxplot(ax, log2_ratio, groups, BMI_COLORS, title,
                            ylabel='log2 (ratio)')
    fig.suptitle('Enzyme activity proxy ratios (n=165, ANCOVA-adjusted for age/GA/year)',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.06, hspace=0.55, wspace=0.30)
    fig_path = FIG_DIR / 'enzyme_ratios_grid.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {fig_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
