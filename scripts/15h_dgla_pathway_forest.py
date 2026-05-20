"""15h. DGLA 池通路比值森林图 — 与 15g (AA 池) 对位.

4 条 DGLA 池通路级聚合比值 (raw 浓度求和 → log2 → ANCOVA), 外加 1 条 FADS1 底物比:
  1. DGLA-COX 总通量 : (PGE1+PGD1+PGF1α + 15-keto-PGE1 + 13,14-dh-15k-PGF1α) / DGLA
  2. DGLA-12-LOX      : 12-HETrE / DGLA
  3. DGLA-15-LOX      : 15-HETrE / DGLA
  4. DGLA-8-LOX/Auto  : 8-HETrE / DGLA
  5. FADS1 底物比      : DGLA / AA   (Δ5-去饱和酶通量代理: ↓ = FADS1 活跃, DGLA 被消耗)

数据源:
  - 9 个 DGLA-oxylipin    : data/02_preprocessed/ori_n165_filtered80_imputed.xlsx
  - DGLA (Dihomo-γ-linolenic acid) 母体: data/01_raw/ori_allmet_QC.xlsx
                                          (panel 缺, 从 673-met 全表抽取; 100% 检出, CV=37%)
  - AA 母体              : data/02_preprocessed/ori_n165_filtered80_imputed.xlsx (panel 内)

注意:
  - 与 15g 同酶同色 (12-LOX 深绿, 15-LOX 浅绿, COX 红, Auto 橙, FADS1 灰).
  - 旧版本曾用 AA 作 DGLA 通路的分母代理 (2026-05-20 修正; 详见 commit log).

模型: log2(ratio) ~ BMI + age + GA + year_A19 + year_C21 + batch_2
图: 横向 forest plot — 圆点 + 95% CI + 右侧 p 标注 + 0 处虚线; 通路色

输出:
  results/tables/dgla_pathway_ratios_forest.csv
  results/figures/04_Metabolite_boxplots/dgla_pathway_forest.png
  results/figures/04_Metabolite_boxplots/dgla_pathway_forest.svg (矢量, 论文用)
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
    'AA':                       'Arachidonic acid',
    'DGLA':                     'Dihomo-γ-linolenic acid',   # from ori_allmet_QC.xlsx (panel 缺)
    # DGLA-COX (PG-1 系)
    'PGE1':                     'Prostaglandin E1',
    'PGD1':                     'Prostaglandin D1',
    'PGF1α':                    'Prostaglandin F1α',
    '15-keto-PGE1':             '15-Keto prostaglandin E1',
    '13,14-dh-15k-PGF1α':       '13,14-Dihydro-15-keto Prostaglandin F1α',
    # DGLA-LOX (HETrE 系)
    '12-HETrE':                 '12-Hydroxy-8Z,10E,14Z-eicosatrienoic acid',
    '15-HETrE':                 '15-Hydroxy-8Z,11Z,13E-eicosatrienoic acid',
    '8-HETrE':                  '8-Hydroxy-9E,11Z,14Z-eicosatrienoic acid',
}

# 通路定义 (与 15g 同酶同色: COX 红, 12-LOX 深绿, 15-LOX 浅绿, Auto 橙; FADS1 灰)
PATHWAYS = [
    ('COX 总通量',  '#D62728',
     ['PGE1', 'PGD1', 'PGF1α', '15-keto-PGE1', '13,14-dh-15k-PGF1α'],
     ['DGLA']),
    ('12-LOX',     '#1B7E3A',
     ['12-HETrE'], ['DGLA']),
    ('15-LOX',     '#5BA85B',
     ['15-HETrE'], ['DGLA']),
    ('8-LOX/Auto', '#FF7F0E',
     ['8-HETrE'], ['DGLA']),
    ('DGLA/AA',    '#666666',
     ['DGLA'], ['AA']),
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
    print('=== 15h. DGLA 池通路森林图 ===\n')

    df80 = pd.read_excel(ROOT / 'data/02_preprocessed/ori_n165_filtered80_imputed.xlsx')
    sample_cols = [c for c in df80.columns if c not in META_COLS and c != UNIT_COL]
    mat80 = df80.set_index('Metabolite Name')[sample_cols]

    # DGLA 母体: panel 缺, 从 ori_allmet_QC.xlsx 全表抽取 (与 15_enzyme_ratios.py EPA/DHA 同模式)
    df_all = pd.read_excel(ROOT / 'data/01_raw/ori_allmet_QC.xlsx')
    all_sample_cols = [c for c in df_all.columns if c not in META_COLS and c != UNIT_COL]
    missing = [s for s in sample_cols if s not in all_sample_cols]
    if missing:
        raise RuntimeError(f'  {len(missing)} 主队列样本不在 QC 表内: {missing[:5]}')
    dgla_row = df_all[df_all['Metabolite Name'] == NAME['DGLA']].set_index(
        'Metabolite Name')[all_sample_cols][sample_cols]
    vals = pd.to_numeric(dgla_row.iloc[0], errors='coerce')
    nz = vals[(vals > 0) & vals.notna()]
    print(f'  DGLA 检出 {len(nz)}/{len(sample_cols)} ({len(nz)/len(sample_cols)*100:.1f}%) '
          f'中位 {nz.median():.3g} ng/g, CV={vals.std()/vals.mean()*100:.1f}%')
    mat = pd.concat([mat80, dgla_row])

    cov_idx = build_covariates(sample_cols)
    groups = cov_idx['BMI_group'].values
    X = sm.add_constant(cov_idx[['BMI_overweight', 'age', 'GA_decimal',
                                  'year_A19', 'year_C21', 'batch_2']].astype(float))

    rows = []
    print(f'  {"通路":<18} {"log2FC":>7} {"FC":>5}  {"95% CI":>22}  {"p_OLS":>8}  {"p_Wilc":>8}')
    print('  ' + '-' * 82)
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
        print(f'  {dispname:<18} {beta:+.3f}  {2**beta:.2f}× {ci_str:>22}  '
              f'{p:>8.4f}  {float(p_wil):>8.4f}')

    res = pd.DataFrame(rows)
    out_csv = ROOT / 'results/tables/dgla_pathway_ratios_forest.csv'
    res.to_csv(out_csv, index=False, encoding='utf-8-sig')
    print(f'\n  [OK] {out_csv.relative_to(ROOT)}')

    # ===== Forest plot (与 aa_pathway_forest 同风格) =====
    fig, ax = plt.subplots(figsize=(13, 5.2))
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


    plt.tight_layout()
    out_png = ROOT / 'results/figures/04_Metabolite_boxplots/dgla_pathway_forest.png'
    out_svg = ROOT / 'results/figures/04_Metabolite_boxplots/dgla_pathway_forest.svg'
    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(out_svg, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  [OK] {out_png.relative_to(ROOT)}')
    print(f'  [OK] {out_svg.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
