"""15h. LA 池通路比值森林图 — 仿 dha_pathway_revised_forest / 15g 风格.

LA (亚油酸 18:2 ω-6) 池在 panel 内完整覆盖, 不需从全表补救.
没有 COX (LA 不是 COX 底物), 没有 CYP4-ω (panel 不收 LA ω 产物).

6 条 LA 池通路级聚合比值 (raw 浓度求和 → log2 → ANCOVA):
  1. CYP-Epo          : (9,10 + 12,13-EpOME) / LA                   位置配对完整
  2. sEH              : (9,10 + 12,13-DiHOME) / (9,10 + 12,13-EpOME) 完整配对!
  3. 13-LOX 酶为主     : 13-HODE / LA                                evidence=4 substantial
  4. 非酶 (9-HODE)    : 9-HODE / LA                                  evidence=2 dominant
  5. HpODE 中间体     : 13-HpODE / LA                                evidence=4 substantial
  6. oxoODE 失活/进一步: (9 + 13-oxoODE) / LA                         evidence=3 substantial

模型: log2(ratio) ~ BMI + age + GA + year_A19 + year_C21 + batch_2
图: 横向 forest plot — 圆点 + 95% CI + 右侧 p; 通路色

输出:
  results/tables/la_pathway_ratios_forest.csv
  results/figures/04_Metabolite_boxplots/la_pathway_forest.png + svg
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
    'LA':            'Linoleic acid',
    '9-HODE':        '9-Hydroxy-10E,12Z-octadecadienoic acid',
    '13-HODE':       '13-Hydroxy-9Z,11E-octadecadienoic acid',
    '13-HpODE':      '13-Hydroperoxy-9Z,11E-octadecadienoic acid',
    '9-oxoODE':      '9-Oxo-10E,12Z-octadecadienoic acid',
    '13-oxoODE':     '13-Oxo-9Z,11E-octadecadienoicacid',
    '9,10-EpOME':    '9,10-Epoxy-12Z-octadecenoic acid',
    '12,13-EpOME':   '12,13-Epoxy-9Z-octadecenoic acid',
    '9,10-DiHOME':   '9,10-DiHydroxy-12Z-octadecenoic acid',
    '12,13-DiHOME':  '12,13 -DiHydroxy-9Z-octadecenoic acid',
}

PATHWAYS = [
    # display_name, color, numerator_keys, denominator_keys
    ('CYP-Epo',         '#1F77B4',
     ['9,10-EpOME','12,13-EpOME'], ['LA']),
    ('sEH',             '#9467BD',
     ['9,10-DiHOME','12,13-DiHOME'], ['9,10-EpOME','12,13-EpOME']),
    ('13-LOX (酶为主)',  '#1B7E3A',
     ['13-HODE'], ['LA']),
    ('非酶氧化 (9-HODE)', '#FF7F0E',
     ['9-HODE'], ['LA']),
    ('HpODE 中间体',     '#9C27B0',
     ['13-HpODE'], ['LA']),
    ('oxoODE 失活',      '#795548',
     ['9-oxoODE','13-oxoODE'], ['LA']),
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
    print('=== 15h. LA 池通路森林图 ===\n')

    df80 = pd.read_excel(ROOT / 'data/02_preprocessed/ori_n165_filtered80_imputed.xlsx')
    sample_cols = [c for c in df80.columns if c not in META_COLS and c != UNIT_COL]
    mat = df80.set_index('Metabolite Name')[sample_cols]
    # 验证全部在 80 主轨
    for k, v in NAME.items():
        if v not in mat.index:
            raise RuntimeError(f'  {k} 不在 80 主轨: {v}')

    cov_idx = build_covariates(sample_cols)
    groups = cov_idx['BMI_group'].values
    X = sm.add_constant(cov_idx[['BMI_overweight', 'age', 'GA_decimal',
                                  'year_A19', 'year_C21', 'batch_2']].astype(float))

    rows = []
    print(f'  {"通路":<18} {"log2FC":>7} {"FC":>5}  {"95% CI":>22}  {"p_OLS":>8}  {"p_Wilc":>8}')
    print('  ' + '-' * 84)
    for dispname, color, num_keys, den_keys in PATHWAYS:
        nums = [NAME[k] for k in num_keys]
        dens = [NAME[k] for k in den_keys]
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
        print(f'  {dispname:<18} {beta:+.3f}  {2**beta:.2f}× {ci_str:>22}  '
              f'{p:>8.4f}  {float(p_wil):>8.4f}')

    res = pd.DataFrame(rows)
    res.to_csv(ROOT / 'results/tables/la_pathway_ratios_forest.csv',
               index=False, encoding='utf-8-sig')
    print(f'\n  ✓ results/tables/la_pathway_ratios_forest.csv')

    # ===== Forest plot =====
    fig, ax = plt.subplots(figsize=(13, 6.5))
    n = len(res)
    ys = np.arange(n)[::-1]

    for i, r in res.iterrows():
        y = ys[i]
        c = r['color']
        ax.errorbar(r['log2FC'], y,
                    xerr=[[r['log2FC'] - r['CI_lo']], [r['CI_hi'] - r['log2FC']]],
                    fmt='o', color=c, markersize=14, capsize=8,
                    linewidth=2.5, capthick=2.0, markeredgewidth=0)
        ax.text(r['CI_hi'] + 0.04, y, f"p = {r['p_ols_hc3']:.4f}",
                va='center', fontsize=11.5, color=c, fontweight='bold')

    ax.axvline(0, color='#444444', linewidth=1.0, linestyle='--', alpha=0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels([])
    for i, r in res.iterrows():
        ax.text(-0.02, ys[i], r['pathway'],
                transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=14, fontweight='bold',
                color=r['color'])

    ax.set_xlabel('log2FC', fontsize=13)
    xlo = min(res['CI_lo'].min(), -0.05) - 0.1
    xhi = max(res['CI_hi'].max(), 0.05) + 0.4
    ax.set_xlim(xlo, xhi)

    ax.grid(axis='x', alpha=0.35, linestyle='--', color='#bbbbbb', linewidth=0.8)
    ax.set_axisbelow(True)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_color('#666666')
    ax.spines['bottom'].set_color('#666666')

    ax.set_title('LA 池 6 条通路比值 ANCOVA — 通路级聚合效应量森林图\n'
                 '(超重肥胖 vs 正常, n=165; 调整 age/GA/year/batch; raw 浓度比值; '
                 'panel 内 LA 系完整覆盖)',
                 fontsize=12, pad=14)

    plt.tight_layout()
    out_png = ROOT / 'results/figures/04_Metabolite_boxplots/la_pathway_forest.png'
    out_svg = ROOT / 'results/figures/04_Metabolite_boxplots/la_pathway_forest.svg'
    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(out_svg, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  ✓ {out_png.relative_to(ROOT)}')
    print(f'  ✓ {out_svg.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
