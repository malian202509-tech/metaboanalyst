"""15c. 16-HDoHE 归属实证检验 — CYP-ω 还是非酶自氧化?

文献争议:
  - panel/LIPID MAPS 注释偏向 CYP-ω
  - 但 VanRollins 2008 (CYP-DHA 主产物 19/20/22-HDHA) + Yin & Porter 2005
    (DHA 自由基氧化 C15 H 脱氢 → 16-OH + 17E 双键重排) 偏向非酶

实证思路:
  跑 8 个 HDoHE 单分子 ANCOVA, 用 effect size + p 值的相对位置看
  16-HDoHE 与 20-HDoHE (CYP) 还是与 7/10/11/13-HDoHE (非酶) 更接近.

输出:
  results/tables/hdohe_singlemolecule_ancova.csv
  results/figures/04_Metabolite_boxplots/hdohe_effect_size_landscape.png
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

# 8 HDoHE in 80 main + 17-HDoHE/4-HDoHE from QC全表
HDOHE_PANEL = {
    # name, full panel string, in_80_main, hypothesized origin
    '14-HDoHE':  ('14-Hydroxy-4Z,7Z,10Z,12E,16Z,19Z-docosahexaenoic acid', True,  'LOX (12-LOX, maresin前体)'),
    '17-HDoHE':  ('17-Hydroxy-4Z,7Z,10Z,13Z,15E,19Z-docosahexaenoic acid', False, 'LOX (15-LOX-1, D-resolvin前体)'),
    '20-HDoHE':  ('20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid', True,  'CYP4 (ω-3, pure)'),
    '16-HDoHE':  ('16-Hydroxy-4Z,7Z,10Z,13Z,17E,19Z-docosahexaenoic acid', True,  'CYP4? OR Non-enzymatic? (争议)'),
    '4-HDoHE':   ('4-Hydroxy-5E,7Z,10Z,13Z,16Z,19Z-docosahexaenoic acid',  False, 'Mixed (5-LOX 副产物 + 非酶)'),
    '7-HDoHE':   ('7-Hydroxy-4Z,8E,10Z,13Z,16Z,19Z-docosahexaenoic acid',  True,  'Non-enzymatic'),
    '10-HDoHE':  ('10-Hydroxy-4Z,7Z,11E,13Z,16Z,19Z-docosahexaenoic acid', True,  'Non-enzymatic'),
    '11-HDoHE':  ('11-Hydroxy-4Z,7Z,9E,13Z,16Z,19Z-docosahexaenoic acid',  True,  'Non-enzymatic'),
    '13-HDoHE':  ('13-Hydroxy-4Z,7Z,10Z,14E,16Z,19Z-docosahexaenoic acid', True,  'Non-enzymatic'),
}


def main():
    print('=== 15c. 16-HDoHE 归属实证检验 ===\n')

    df80 = pd.read_excel(ROOT / 'data/02_preprocessed/ori_n165_filtered80_imputed.xlsx')
    sample_cols = [c for c in df80.columns if c not in META_COLS and c != UNIT_COL]
    mat80 = df80.set_index('Metabolite Name')[sample_cols]

    df_all = pd.read_excel(ROOT / 'data/01_raw/ori_allmet_QC.xlsx')
    all_cols = [c for c in df_all.columns if c not in META_COLS and c != UNIT_COL]

    # Collect all 9 HDoHE rows (raw scale)
    raw_rows = {}
    for short, (full_name, in_80, origin) in HDOHE_PANEL.items():
        if in_80:
            raw_rows[short] = mat80.loc[full_name].astype(float).values
        else:
            row = (df_all[df_all['Metabolite Name'] == full_name]
                   .set_index('Metabolite Name')[all_cols][sample_cols].iloc[0])
            row = pd.to_numeric(row, errors='coerce').astype(float)
            nz = row[row > 0]
            if (~(row > 0)).sum() > 0:
                half_min = float(nz.min()) / 2.0
                row = row.where(row > 0, half_min)
            raw_rows[short] = row.values

    # Covariates (same model: BMI + age + GA + year + batch)
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
    cov_idx = cov.set_index('omx_id').loc[sample_cols]
    groups = cov_idx['BMI_group'].values
    X = sm.add_constant(cov_idx[['BMI_overweight', 'age', 'GA_decimal',
                                  'year_A19', 'year_C21', 'batch_2']].astype(float))

    print(f'{"分子":<12} {"假设来源":<35} {"log2FC":>8} {"FC":>5} '
          f'{"95% CI":>22} {"p_OLS":>8} {"p_Wilc":>8}')
    print('-' * 105)
    rows = []
    for short, (_, _, origin) in HDOHE_PANEL.items():
        y = np.log2(raw_rows[short])
        m = sm.OLS(y, X).fit(cov_type='HC3')
        beta = float(m.params['BMI_overweight'])
        p = float(m.pvalues['BMI_overweight'])
        ci_lo, ci_hi = [float(x) for x in m.conf_int().loc['BMI_overweight'].tolist()]
        _, p_wil = stats.mannwhitneyu(y[groups == '超重肥胖'], y[groups == '正常'],
                                       alternative='two-sided')
        rows.append({
            'short': short, 'hypothesized_origin': origin,
            'log2FC': round(beta, 4),
            'CI_lo': round(ci_lo, 4), 'CI_hi': round(ci_hi, 4),
            'fold_change_linear': round(2 ** beta, 4),
            'p_ols_hc3': round(p, 6),
            'p_wilcoxon': round(float(p_wil), 6),
        })
        ci_str = f'({ci_lo:+.3f},{ci_hi:+.3f})'
        print(f'{short:<12} {origin[:33]:<35} {beta:+.3f}  '
              f'{2**beta:.2f}× {ci_str:>22} {p:>8.4f} {float(p_wil):>8.4f}')

    res = pd.DataFrame(rows)
    res.to_csv(ROOT / 'results/tables/hdohe_singlemolecule_ancova.csv',
               index=False, encoding='utf-8-sig')
    print('\n  ✓ results/tables/hdohe_singlemolecule_ancova.csv')

    # Forest plot
    fig, ax = plt.subplots(figsize=(11, 6.5))
    cmap = {
        'LOX (12-LOX, maresin前体)':         '#2CA02C',
        'LOX (15-LOX-1, D-resolvin前体)':     '#2CA02C',
        'CYP4 (ω-3, pure)':                   '#1F77B4',
        'CYP4? OR Non-enzymatic? (争议)':     '#FF7F0E',
        'Mixed (5-LOX 副产物 + 非酶)':        '#888888',
        'Non-enzymatic':                      '#D62728',
    }
    sorted_rows = sorted(rows, key=lambda r: r['log2FC'])
    ys = np.arange(len(sorted_rows))
    for i, r in enumerate(sorted_rows):
        color = cmap.get(r['hypothesized_origin'], 'black')
        ax.errorbar(r['log2FC'], i,
                    xerr=[[r['log2FC'] - r['CI_lo']], [r['CI_hi'] - r['log2FC']]],
                    fmt='o', color=color, markersize=10, capsize=4, linewidth=2)
        annot = f"{r['short']}  p={r['p_ols_hc3']:.3f}"
        ax.text(r['CI_hi'] + 0.03, i, annot, va='center', fontsize=9.5, color=color)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_yticks(ys)
    ax.set_yticklabels([r['short'] for r in sorted_rows], fontsize=10)
    ax.set_xlabel('log2 fold-change (超重肥胖 vs 正常)  +95% CI', fontsize=10)
    ax.set_title('9 HDoHE 单分子 ANCOVA 效应量 — 看 16-HDoHE 与谁聚类\n'
                 '(绿=LOX SPM前体  蓝=纯 CYP4-ω  橙=争议 16-HDoHE  红=非酶  灰=混合)',
                 fontsize=11.5, pad=10)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_xlim(-0.25, max(r['CI_hi'] for r in sorted_rows) + 0.55)
    plt.tight_layout()
    out_fig = ROOT / 'results/figures/04_Metabolite_boxplots/hdohe_effect_size_landscape.png'
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {out_fig.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
