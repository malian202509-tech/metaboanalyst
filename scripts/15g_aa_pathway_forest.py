"""15g. AA 池通路比值森林图 — 仿 dha_pathway_revised_forest 极简风格.

6 条 AA 池通路级聚合比值 (raw 浓度求和 → log2 → ANCOVA):
  1. COX 总通量      : total_PG_AA / AA                              (已在 15)
  2. 12-LOX          : 12-HETE / AA                                  (已在 15)
  3. 15-LOX          : 15-HETE / AA                                  (已在 15)
  4. CYP-Epo         : (8,9 + 11,12-EpETrE) / AA                     (新算, panel 仅 2 个位置)
  5. sEH             : (8,9 + 11,12-DiHETrE) / (8,9 + 11,12-EpETrE)  (新算, 位置配对)
  6. CYP4-ω          : (16 + 18-HETE) / AA                           (新算)
  7. 非酶氧化        : (11 + 8-HETE) / AA                            (新算)

模型: log2(ratio) ~ BMI + age + GA + year_A19 + year_C21 + batch_2
图: 横向 forest plot — 圆点 + 95% CI + 右侧 p 标注 + 0 处虚线; 通路色

输出:
  results/tables/aa_pathway_ratios_forest.csv
  results/figures/04_Metabolite_boxplots/aa_pathway_forest.png
  results/figures/04_Metabolite_boxplots/aa_pathway_forest.svg (矢量, 论文用)
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'

NAME = {
    'AA':           'Arachidonic acid',
    '12-HETE':      '12-Hydroxy-5Z,8Z,10E,14Z-eicosatetraenoic acid',
    '15-HETE':      '15-Hydroxy-5Z,8Z,11Z,13E-eicosatetraenoic acid',
    '11-HETE':      '11-Hydroxy-5Z,8Z,11E,14Z-eicosatetraenoic acid',
    '8-HETE':       '8-Hydroxy-5Z,9E,11Z,14Z-eicosatetraenoic acid',
    '16-HETE':      '16-Hydroxy-5Z,8Z,11Z,14Z-eicosatetraenoic acid',
    '18-HETE':      '18-Hydroxy-5Z,8Z,11Z,14Z-eicosatetraenoic acid',
    '8,9-EpETrE':   '8,9-Dpoxy-5Z,11Z,14Z-eicosatrienoic acid',
    '11,12-EpETrE': '11,12-Epoxy-5Z,8Z,14Z-eicosatrienoic acid',
    '8,9-DiHETrE':  '8,9-DiHydroxy-5Z,11Z,14Z-eicosatrienoic acid',
    '11,12-DiHETrE':'11,12-DiHydroxy-5Z,8Z,14Z-eicosatrienoic acid',
    # COX PG/TX (8 个 AA 系 PG total)
    'PGD2':         'Prostaglandin D2',
    'PGE2':         'Prostaglandin E2',
    'PGF2α':        'Prostaglandin F2α',
    'PGJ2':         'Prostaglandin J2',
    'PGA2':         'Prostaglandin A2',
    '11β-PGE2':     '11β-Prostaglandin E2',
    'TXB2':         'Thromboxane B2',
    '12-HHT':       '12S-Hydroxy-5Z,8E,10E-heptadecatrienoic acid',
}

# 通路定义 (顺序按 DHA 图风格: 从最强到最弱)
# 显示顺序从上到下: CYP4-ω -> 12-LOX -> COX -> CYP-Epo -> sEH -> 15-LOX -> 非酶
PATHWAYS = [
    # display_name, color, numerator_list, denominator_list
    ('COX 总通量',  '#D62728',
     ['PGD2','PGE2','PGF2α','PGJ2','PGA2','11β-PGE2','TXB2','12-HHT'],
     ['AA']),
    ('12-LOX',     '#1B7E3A',
     ['12-HETE'], ['AA']),
    ('15-LOX',     '#5BA85B',
     ['15-HETE'], ['AA']),
    ('CYP-Epo',    '#1F77B4',
     ['8,9-EpETrE','11,12-EpETrE'], ['AA']),
    ('sEH',        '#9467BD',
     ['8,9-DiHETrE','11,12-DiHETrE'], ['8,9-EpETrE','11,12-EpETrE']),
    ('CYP4-ω',     '#17BECF',
     ['16-HETE','18-HETE'], ['AA']),
    ('非酶氧化',    '#FF7F0E',
     ['11-HETE','8-HETE'], ['AA']),
]


def build_covariates(sample_cols):
    align = pd.read_csv(ROOT / 'data/02_preprocessed/sample_alignment_n165.csv')
    clin = pd.read_excel(ROOT / 'data/01_raw/分娩特征分析_2019-2021_主分析队列n165.xlsx')
    bxl = pd.read_excel(ROOT / 'data/01_raw/上机顺序和批次表.xlsx',
                        sheet_name='C18柱-上机顺序和批次表', header=1, usecols=range(4))
    bxl.columns = ['order', 'class', 'omx_id', 'batch']
    bxl['batch'] = pd.to_numeric(bxl['batch'], errors='coerce').astype('Int64')
    btbl = bxl[bxl['class'] == 'Subject'][['omx_id', 'batch']]
    cov = (clin[['标本号', '孕妇年龄', 'GA_decimal', 'BMI_group']]
           .rename(columns={'孕妇年龄': 'age'})
           .merge(align[['标本号', 'omx_id']], on='标本号', how='left')
           .merge(btbl, on='omx_id', how='left'))
    cov['BMI_overweight'] = (cov['BMI_group'] == '超重肥胖').astype(float)
    yp = cov['omx_id'].str[:3]
    cov['year_A19'] = (yp == 'A19').astype(float)
    cov['year_C21'] = (yp == 'C21').astype(float)
    cov['batch_2'] = (cov['batch'].astype(float) == 2.0).astype(float)
    return cov.set_index('omx_id').loc[sample_cols]


def main():
    print('=== 15g. AA 池通路森林图 ===\n')

    df80 = pd.read_excel(ROOT / 'data/02_preprocessed/ori_n165_filtered80_imputed.xlsx')
    sample_cols = [c for c in df80.columns if c not in META_COLS and c != UNIT_COL]
    mat = df80.set_index('Metabolite Name')[sample_cols]

    cov_idx = build_covariates(sample_cols)
    groups = cov_idx['BMI_group'].values
    X = sm.add_constant(cov_idx[['BMI_overweight', 'age', 'GA_decimal',
                                  'year_A19', 'year_C21', 'batch_2']].astype(float))

    rows = []
    print(f'  {"通路":<14} {"log2FC":>7} {"FC":>5}  {"95% CI":>22}  {"p_OLS":>8}  {"p_Wilc":>8}')
    print('  ' + '-' * 78)
    for dispname, color, num_keys, den_keys in PATHWAYS:
        nums = [NAME[k] for k in num_keys]
        dens = [NAME[k] for k in den_keys]
        miss = [n for n in nums + dens if n not in mat.index]
        if miss:
            raise RuntimeError(f'  {dispname}: 缺 {miss}')
        n_sum = mat.loc[nums].astype(float).sum(axis=0).values
        d_sum = mat.loc[dens].astype(float).sum(axis=0).values
        y = np.log2(n_sum) - np.log2(d_sum)
        m = sm.OLS(y, X).fit(cov_type='HC3')
        beta = float(m.params['BMI_overweight'])
        p = float(m.pvalues['BMI_overweight'])
        ci_lo, ci_hi = [float(x) for x in m.conf_int().loc['BMI_overweight'].tolist()]
        _, p_wil = stats.mannwhitneyu(y[groups == '超重肥胖'], y[groups == '正常'],
                                       alternative='two-sided')
        rows.append({
            'pathway': dispname, 'color': color,
            'numerator': '+'.join(num_keys), 'denominator': '+'.join(den_keys),
            'log2FC': round(beta, 4), 'CI_lo': round(ci_lo, 4),
            'CI_hi': round(ci_hi, 4), 'fold_change_linear': round(2 ** beta, 4),
            'p_ols_hc3': round(p, 6), 'p_wilcoxon': round(float(p_wil), 6),
        })
        ci_str = f'({ci_lo:+.3f},{ci_hi:+.3f})'
        print(f'  {dispname:<14} {beta:+.3f}  {2**beta:.2f}× {ci_str:>22}  '
              f'{p:>8.4f}  {float(p_wil):>8.4f}')

    res = pd.DataFrame(rows)
    res.to_csv(ROOT / 'results/tables/aa_pathway_ratios_forest.csv',
               index=False, encoding='utf-8-sig')
    print(f'\n  ✓ results/tables/aa_pathway_ratios_forest.csv')

    # ===== Forest plot (仿 dha_pathway_revised_forest 风格) =====
    fig, ax = plt.subplots(figsize=(13, 6.5))
    n = len(res)
    ys = np.arange(n)[::-1]  # 从上到下: 第 1 条最上

    for i, r in res.iterrows():
        y = ys[i]
        c = r['color']
        # 误差线 (两端短竖杠) + 圆点
        ax.errorbar(r['log2FC'], y,
                    xerr=[[r['log2FC'] - r['CI_lo']], [r['CI_hi'] - r['log2FC']]],
                    fmt='o', color=c, markersize=14, capsize=8,
                    linewidth=2.5, capthick=2.0, markeredgewidth=0)
        # 右侧 p 值标注 (与点同色, 粗体)
        ax.text(r['CI_hi'] + 0.04, y, f"p = {r['p_ols_hc3']:.4f}",
                va='center', fontsize=11.5, color=c, fontweight='bold')

    # 0 处垂直虚线
    ax.axvline(0, color='#444444', linewidth=1.0, linestyle='--', alpha=0.6)

    # y 轴: 通路名 (彩色 + 粗体)
    ax.set_yticks(ys)
    ax.set_yticklabels([])  # 自定义
    for i, r in res.iterrows():
        ax.text(-0.02, ys[i], r['pathway'],
                transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=14, fontweight='bold',
                color=r['color'])

    # x 轴
    ax.set_xlabel('log2FC', fontsize=13)
    # 自动扩展 x 范围, 给右侧 p 标注留空间
    xlo = min(res['CI_lo'].min(), -0.05) - 0.1
    xhi = max(res['CI_hi'].max(), 0.05) + 0.4
    ax.set_xlim(xlo, xhi)

    # 网格 (浅灰虚线, 仅 x 方向; y 方向不画)
    ax.grid(axis='x', alpha=0.35, linestyle='--', color='#bbbbbb', linewidth=0.8)
    ax.set_axisbelow(True)

    # 去掉上 + 右 spine, 留左下
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_color('#666666')
    ax.spines['bottom'].set_color('#666666')

    # 顶部小标注 (可选 — 题外信息)
    ax.set_title('AA 池 7 条通路比值 ANCOVA — 通路级聚合效应量森林图\n'
                 '(超重肥胖 vs 正常, n=165; 调整 age/GA/year/batch; raw 浓度比值)',
                 fontsize=12, pad=14)

    plt.tight_layout()
    out_png = ROOT / 'results/figures/04_Metabolite_boxplots/aa_pathway_forest.png'
    out_svg = ROOT / 'results/figures/04_Metabolite_boxplots/aa_pathway_forest.svg'
    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(out_svg, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  ✓ {out_png.relative_to(ROOT)}')
    print(f'  ✓ {out_svg.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
