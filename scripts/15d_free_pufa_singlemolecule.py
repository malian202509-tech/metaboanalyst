"""15d. 游离 PUFA 单分子 ANCOVA — DHA / EPA / AA / LA / ALA / CLA / tLA.

§7B.6 家族 GSEA 显示 Free PUFA 家族 q=0.004 ↓↓ (最强家族信号),
但 family-level 富集是 5 个 PUFA 一起算 mean t. 单分子层面 DHA, EPA, AA 等
具体方向和显著性需要单独算.

5 个 PUFA 在 80 主轨 (从 ori_n165_filtered80_imputed.xlsx, raw 浓度):
  AA, ALA, CLA, LA, tLA
2 个 PUFA 在 ori_allmet_QC.xlsx 全表 (panel 缺):
  EPA (5Z,8Z,11Z,14Z,17Z-Eicosapentaenoic acid)  100% 检出 CV%=5.5%
  DHA (4Z,7Z,10Z,13Z,16Z,19Z-Docosahexaenoic acid) 100% 检出 CV%=17.4%

模型: log2(浓度) ~ BMI + age + GA_decimal + year + batch

输出:
  results/tables/free_pufa_singlemolecule.csv
  results/figures/04_Metabolite_boxplots/free_pufa_forest.png
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
import matplotlib.patches as mpatches

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

PUFA_IN_PANEL = {
    'AA':  'Arachidonic acid',
    'ALA': 'α-Linolenic acid',
    'LA':  'Linoleic acid',
    'CLA': 'Conjugated linoleic acids',
    'tLA': 'Linoelaidic acid',
}
PUFA_FROM_QC = {
    'EPA': '5Z,8Z,11Z,14Z,17Z-Eicosapentaenoic acid',
    'DHA': '4Z,7Z,10Z,13Z,16Z,19Z-Docosahexaenoic acid',
}
# Display order: ω-6 / ω-3 / 其他
ORDER = ['LA', 'AA', 'ALA', 'EPA', 'DHA', 'CLA', 'tLA']
CLASS = {
    'LA':  'ω-6 前体',
    'AA':  'ω-6 主底物 (COX/LOX 入口)',
    'ALA': 'ω-3 前体',
    'EPA': 'ω-3 中游',
    'DHA': 'ω-3 主底物 (SPM 前体)',
    'CLA': '反式 / 共轭',
    'tLA': '反式异构体',
}


def main():
    print('=== 15d. 游离 PUFA 单分子 ANCOVA ===\n')

    df80 = pd.read_excel(ROOT / 'data/02_preprocessed/ori_n165_filtered80_imputed.xlsx')
    sample_cols = [c for c in df80.columns if c not in META_COLS and c != UNIT_COL]
    mat80 = df80.set_index('Metabolite Name')[sample_cols]

    df_all = pd.read_excel(ROOT / 'data/01_raw/ori_allmet_QC.xlsx')
    all_cols = [c for c in df_all.columns if c not in META_COLS and c != UNIT_COL]
    extra = (df_all[df_all['Metabolite Name'].isin(PUFA_FROM_QC.values())]
             .set_index('Metabolite Name')[all_cols][sample_cols])
    extra = extra.apply(pd.to_numeric, errors='coerce').astype(float)

    # 检出 / 插补
    print('  165 主队列实际检出:')
    for k, v in PUFA_IN_PANEL.items():
        row = pd.to_numeric(mat80.loc[v], errors='coerce')
        n_pos = ((row > 0) & row.notna()).sum()
        nz = row[row > 0]
        print(f'    {k:<5} (80 主轨): {n_pos:>3}/165 ({n_pos/165*100:5.1f}%) 中位 {nz.median():.3f} ng/g')
    for k, v in PUFA_FROM_QC.items():
        row = extra.loc[v]
        n_pos = ((row > 0) & row.notna()).sum()
        nz = row[row > 0]
        print(f'    {k:<5} (QC 全表): {n_pos:>3}/165 ({n_pos/165*100:5.1f}%) 中位 {nz.median():.3f} ng/g')
        if (~(row > 0)).sum() > 0:
            half_min = float(nz.min()) / 2.0
            extra.loc[v] = row.where(row > 0, half_min).astype(float).values

    full = pd.concat([mat80, extra])

    # Covariates
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

    all_pufa = {**PUFA_IN_PANEL, **PUFA_FROM_QC}
    rows = []
    print('\n' + '=' * 105)
    print(f'{"PUFA":<5} {"class":<28} {"log2FC":>7} {"FC":>5} {"95% CI":>22} '
          f'{"p_OLS":>8} {"p_Wilc":>8} dir')
    print('-' * 105)
    for k in ORDER:
        v = all_pufa[k]
        y = np.log2(full.loc[v].astype(float).values)
        m = sm.OLS(y, X).fit(cov_type='HC3')
        beta = float(m.params['BMI_overweight'])
        p = float(m.pvalues['BMI_overweight'])
        ci_lo, ci_hi = [float(x) for x in m.conf_int().loc['BMI_overweight'].tolist()]
        _, p_wil = stats.mannwhitneyu(y[groups == '超重肥胖'], y[groups == '正常'],
                                       alternative='two-sided')
        med_norm = float(np.median(full.loc[v].astype(float).values[groups == '正常']))
        med_over = float(np.median(full.loc[v].astype(float).values[groups == '超重肥胖']))
        direction = '↑' if beta > 0 else '↓'
        rows.append({
            'PUFA': k, 'class': CLASS[k],
            'median_normal_ng_g': round(med_norm, 1),
            'median_overweight_ng_g': round(med_over, 1),
            'log2FC': round(beta, 4),
            'CI_lo': round(ci_lo, 4), 'CI_hi': round(ci_hi, 4),
            'fold_change_linear': round(2.0 ** beta, 4),
            'p_ols_hc3': round(p, 6),
            'p_wilcoxon': round(float(p_wil), 6),
            'direction': direction,
        })
        ci_str = f'({ci_lo:+.3f},{ci_hi:+.3f})'
        sig = ' *' if p < 0.05 else ('   ' if p >= 0.10 else ' .')
        print(f'{k:<5} {CLASS[k]:<28} {beta:+.3f}  {2**beta:.2f}× '
              f'{ci_str:>22} {p:>8.4f}{sig} {float(p_wil):>8.4f}   {direction}')

    res = pd.DataFrame(rows)
    res.to_csv(ROOT / 'results/tables/free_pufa_singlemolecule.csv',
               index=False, encoding='utf-8-sig')
    print('\n  ✓ results/tables/free_pufa_singlemolecule.csv')

    # Forest plot
    fig, ax = plt.subplots(figsize=(11, 6))
    cmap = {
        'ω-6 前体':                    '#D62728',  # 红 (ω-6, 促炎倾向)
        'ω-6 主底物 (COX/LOX 入口)':    '#FF7F0E',  # 橙
        'ω-3 前体':                    '#2CA02C',  # 绿
        'ω-3 中游':                    '#17BECF',
        'ω-3 主底物 (SPM 前体)':        '#1F77B4',  # 蓝
        '反式 / 共轭':                  '#888888',
        '反式异构体':                   '#aaaaaa',
    }
    ys = np.arange(len(res))[::-1]
    for i, (_, r) in enumerate(res.iterrows()):
        y = ys[i]
        color = cmap.get(r['class'], 'black')
        ax.errorbar(r['log2FC'], y,
                    xerr=[[r['log2FC'] - r['CI_lo']], [r['CI_hi'] - r['log2FC']]],
                    fmt='o', color=color, markersize=11, capsize=4, linewidth=2)
        sig_mark = ' *' if r['p_ols_hc3'] < 0.05 else ''
        ax.text(r['CI_hi'] + 0.02, y, f"{r['PUFA']}  p={r['p_ols_hc3']:.3f}{sig_mark}",
                va='center', fontsize=10, color=color)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_yticks(ys)
    ax.set_yticklabels([r['PUFA'] for _, r in res.iterrows()], fontsize=11)
    ax.set_xlabel('log2 fold-change (超重肥胖 vs 正常)  ±95% CI', fontsize=10)
    ax.set_title('7 个游离 PUFA 单分子 ANCOVA — DHA 浓度组间对比\n'
                 '(* = raw p<0.05; n=165; 调整 age/GA/year/batch)',
                 fontsize=11.5, pad=10)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    legend_handles = [mpatches.Patch(color=c, label=l) for l, c in cmap.items()]
    ax.legend(handles=legend_handles, loc='lower left', fontsize=8.5,
              frameon=True, ncol=2)
    plt.tight_layout()
    out_fig = ROOT / 'results/figures/04_Metabolite_boxplots/free_pufa_forest.png'
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {out_fig.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
