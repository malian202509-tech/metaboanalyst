"""15b. DHA 通路最终图景 (2026-05-19).

补回供应商剔除的 17-HDoHE 与 4-HDoHE + 16-HDoHE 重归非酶后, 把 DHA 池的氧化
按 4 条通路 + 1 条对照重新呈现:

  ★ 主比值 (5 个, 都画 boxplot):
    1. CYP4-ω 严格版:    20-HDoHE / DHA              ω-3 (CYP4 偏好位置)
    2. LOX SPM 前体:     (14+17-HDoHE) / DHA         12-LOX maresin + 15-LOX D-resolvin
    3. LOX 15 单酶:      17-HDoHE / DHA              15-LOX-1 (D-resolvin 前体)
    4. 自氧化主指数:     (7+10+11+13+16-HDoHE) / DHA 最干净非酶池 (16-HDoHE 重归后)
    5. 5-LOX 对照:       4-HDoHE 单分子              5-LOX 副产物 + 非酶混合 (unchanged)

  [ref] supplementary 行 (不画, 仅入表):
    A. 自氧化 v1 旧版:   (7+10+11+13-HDoHE) / DHA    无 4 无 16
    B. 自氧化 v2:        (4+7+10+11+13-HDoHE) / DHA  含 4 无 16
    C. 自氧化 v3:        (4+7+10+11+13+16) / DHA     全收
    D. 旧 CYP 混合:      (16+20-HDoHE) / DHA         16 归 CYP 旧版

数据来源 + 方法学与 15_enzyme_ratios.py 完全一致:
  - 67 个 80 主轨 metabolite: ori_n165_filtered80_imputed.xlsx (QRILC 插补 raw)
  - 17/4-HDoHE + DHA: ori_allmet_QC.xlsx 抽取, 17-HDoHE 半最小值补 38 个 0
  - raw 浓度求和 → log2 → ANCOVA `~ BMI + age + GA + year + batch`
  - OLS+HC3 + Wilcoxon 互校

输出:
  results/tables/lox_dha_and_autoox_corrected.csv   (9 行: 5 主 + 4 supp)
  results/figures/04_Metabolite_boxplots/lox_dha_and_autoox_v2.png  (1×5 网格)
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
BMI_COLORS = {'正常': '#67A9CF', '超重肥胖': '#EF8A62'}
# 通路色与 15_enzyme_ratios.py 一致
PATH_COLORS = {
    'CYP': '#1F77B4',         # 蓝
    'LOX-SPM': '#2CA02C',     # 绿
    'LOX-15': '#2CA02C',
    'Non-enz': '#D62728',     # 红
    '5-LOX-ref': '#888888',   # 灰
}

HDoHE_80 = {
    '14-HDoHE': '14-Hydroxy-4Z,7Z,10Z,12E,16Z,19Z-docosahexaenoic acid',
    '7-HDoHE':  '7-Hydroxy-4Z,8E,10Z,13Z,16Z,19Z-docosahexaenoic acid',
    '10-HDoHE': '10-Hydroxy-4Z,7Z,11E,13Z,16Z,19Z-docosahexaenoic acid',
    '11-HDoHE': '11-Hydroxy-4Z,7Z,9E,13Z,16Z,19Z-docosahexaenoic acid',
    '13-HDoHE': '13-Hydroxy-4Z,7Z,10Z,14E,16Z,19Z-docosahexaenoic acid',
    '16-HDoHE': '16-Hydroxy-4Z,7Z,10Z,13Z,17E,19Z-docosahexaenoic acid',
    '20-HDoHE': '20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid',
}
EXTRA = {
    '17-HDoHE': '17-Hydroxy-4Z,7Z,10Z,13Z,15E,19Z-docosahexaenoic acid',
    '4-HDoHE':  '4-Hydroxy-5E,7Z,10Z,13Z,16Z,19Z-docosahexaenoic acid',
    'DHA':      '4Z,7Z,10Z,13Z,16Z,19Z-Docosahexaenoic acid',
}


def main():
    print('=== 15b. DHA 通路最终图景 (16-HDoHE 重归非酶后) ===\n')

    df80 = pd.read_excel(ROOT / 'data/02_preprocessed/ori_n165_filtered80_imputed.xlsx')
    sample_cols = [c for c in df80.columns if c not in META_COLS and c != UNIT_COL]
    mat80 = df80.set_index('Metabolite Name')[sample_cols]
    for v in HDoHE_80.values():
        assert v in mat80.index, f'MISSING in 80 main: {v}'

    df_all = pd.read_excel(ROOT / 'data/01_raw/ori_allmet_QC.xlsx')
    all_cols = [c for c in df_all.columns if c not in META_COLS and c != UNIT_COL]
    extra = (df_all[df_all['Metabolite Name'].isin(EXTRA.values())]
             .set_index('Metabolite Name')[all_cols][sample_cols])
    extra = extra.apply(pd.to_numeric, errors='coerce').astype(float)

    print('  165 主队列实际检出:')
    for k, v in EXTRA.items():
        row = extra.loc[v]
        n_pos = ((row > 0) & row.notna()).sum()
        nz = row[row > 0]
        print(f'    {k:<10}: {n_pos:>3}/165 ({n_pos/165*100:5.1f}%)  '
              f'中位 {nz.median():.3f} ng/g')

    for k, v in EXTRA.items():
        row = extra.loc[v]
        nz = row[row > 0]
        n_imp = (~(row > 0)).sum()
        if n_imp > 0:
            half_min = float(nz.min()) / 2.0
            extra.loc[v] = row.where(row > 0, half_min).astype(float).values
            print(f'    半最小值插补 {k}: {n_imp} 个 0 → {half_min:.4f}')

    full = pd.concat([mat80, extra])

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

    # 5 主比值 + 4 supplementary
    ratios = [
        # 主比值
        ('main', 'CYP', 'CYP4-ω 严格: 20-HDoHE / DHA',
         'CYP4 ω-3 羟化通量 (DHA 系, 严格版; 16-HDoHE 已剔除)',
         [HDoHE_80['20-HDoHE']], [EXTRA['DHA']]),
        ('main', 'LOX-SPM', 'LOX SPM 前体: (14+17-HDoHE) / DHA',
         '12-LOX (maresin) + 15-LOX (D-resolvin) 前体生成通量',
         [HDoHE_80['14-HDoHE'], EXTRA['17-HDoHE']], [EXTRA['DHA']]),
        ('main', 'LOX-15', 'LOX 15 单酶: 17-HDoHE / DHA',
         '15-LOX-1 单独通量 (D-resolvin 前体)',
         [EXTRA['17-HDoHE']], [EXTRA['DHA']]),
        ('main', 'Non-enz', '自氧化主指数: (7+10+11+13+16-HDoHE) / DHA',
         'DHA 池非酶自氧化通量 (最干净版; 16-HDoHE 现归非酶)',
         [HDoHE_80['7-HDoHE'], HDoHE_80['10-HDoHE'], HDoHE_80['11-HDoHE'],
          HDoHE_80['13-HDoHE'], HDoHE_80['16-HDoHE']], [EXTRA['DHA']]),
        ('main', '5-LOX-ref', '5-LOX 通路对照: 4-HDoHE / DHA',
         '5-LOX 副产物 (混合非酶); 对照 5-LOX 通路是否激活',
         [EXTRA['4-HDoHE']], [EXTRA['DHA']]),
        # Supplementary (入表不画)
        ('supp', '-', '[supp] Auto-ox v1 旧版: (7+10+11+13-HDoHE) / DHA',
         '旧版自氧化指数 (无 4 无 16)',
         [HDoHE_80['7-HDoHE'], HDoHE_80['10-HDoHE'],
          HDoHE_80['11-HDoHE'], HDoHE_80['13-HDoHE']], [EXTRA['DHA']]),
        ('supp', '-', '[supp] Auto-ox v2: (4+7+10+11+13-HDoHE) / DHA',
         '含 4 无 16',
         [EXTRA['4-HDoHE'], HDoHE_80['7-HDoHE'], HDoHE_80['10-HDoHE'],
          HDoHE_80['11-HDoHE'], HDoHE_80['13-HDoHE']], [EXTRA['DHA']]),
        ('supp', '-', '[supp] Auto-ox v3: (4+7+10+11+13+16-HDoHE) / DHA',
         '全收',
         [EXTRA['4-HDoHE'], HDoHE_80['7-HDoHE'], HDoHE_80['10-HDoHE'],
          HDoHE_80['11-HDoHE'], HDoHE_80['13-HDoHE'], HDoHE_80['16-HDoHE']],
         [EXTRA['DHA']]),
        ('supp', '-', '[supp] 旧 CYP 混合: (16+20-HDoHE) / DHA',
         '修订前混合版 (16-HDoHE 旧归 CYP)',
         [HDoHE_80['16-HDoHE'], HDoHE_80['20-HDoHE']], [EXTRA['DHA']]),
    ]

    print('\n' + '=' * 105)
    print(f'{"tier":<5} {"path":<10} {"ratio":<50} {"log2FC":>7} {"FC":>5} '
          f'{"p_OLS":>8} {"p_Wilc":>8}  dir')
    print('-' * 105)
    ratio_matrix = pd.DataFrame(index=sample_cols)
    rows = []
    for tier, path, name, rationale, num_list, den_list in ratios:
        num_sum = full.loc[num_list].sum(axis=0).astype(float).values
        den_sum = full.loc[den_list].sum(axis=0).astype(float).values
        log2_r = np.log2(num_sum) - np.log2(den_sum)
        ratio_matrix[name] = log2_r
        m = sm.OLS(log2_r, X).fit(cov_type='HC3')
        beta = float(m.params['BMI_overweight'])
        p = float(m.pvalues['BMI_overweight'])
        ci_lo, ci_hi = [float(x) for x in m.conf_int().loc['BMI_overweight'].tolist()]
        _, p_wil = stats.mannwhitneyu(log2_r[groups == '超重肥胖'],
                                       log2_r[groups == '正常'],
                                       alternative='two-sided')
        p_wil = float(p_wil)
        direction = '↑' if beta > 0 else '↓'
        rows.append({
            'tier': tier, 'pathway': path, 'ratio': name, 'rationale': rationale,
            'log2FC': round(beta, 4), 'CI_lo': round(ci_lo, 4),
            'CI_hi': round(ci_hi, 4),
            'fold_change_linear': round(2.0 ** beta, 4),
            'p_ols_hc3': round(p, 6), 'p_wilcoxon': round(p_wil, 6),
            'direction': direction,
        })
        print(f'{tier:<5} {path:<10} {name[:48]:<50} {beta:+.3f}  {2**beta:.2f}× '
              f'{p:>8.4f} {p_wil:>8.4f}   {direction}')

    res = pd.DataFrame(rows)
    res.to_csv(ROOT / 'results/tables/lox_dha_and_autoox_corrected.csv',
               index=False, encoding='utf-8-sig')
    print('\n  ✓ results/tables/lox_dha_and_autoox_corrected.csv')

    # 1×5 boxplot grid (主比值)
    main_rows = res[res['tier'] == 'main'].reset_index(drop=True)
    fig, axes = plt.subplots(1, 5, figsize=(22, 6))
    grp_order = ['正常', '超重肥胖']
    for i, r in main_rows.iterrows():
        ax = axes[i]
        log2_r = ratio_matrix[r['ratio']].values
        data = [log2_r[groups == g] for g in grp_order]
        bp = ax.boxplot(data, positions=[0, 1], widths=0.55, patch_artist=True,
                        boxprops=dict(linewidth=1.0),
                        medianprops=dict(color='black', linewidth=1.5),
                        whiskerprops=dict(linewidth=1.0),
                        flierprops=dict(marker='', markersize=0))
        for b, g in zip(bp['boxes'], grp_order):
            b.set_facecolor(BMI_COLORS[g])
            b.set_alpha(0.5)
            b.set_edgecolor(BMI_COLORS[g])
        rng = np.random.default_rng(42)
        for xp, g in enumerate(grp_order):
            yv = log2_r[groups == g]
            xj = rng.normal(xp, 0.07, size=len(yv))
            ax.scatter(xj, yv, s=18, c=BMI_COLORS[g], edgecolors='black',
                       linewidths=0.4, alpha=0.6, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'{g}\n(n={int((groups==g).sum())})'
                            for g in grp_order], fontsize=10)
        ax.set_ylabel('log2 (ratio)', fontsize=10)
        sig_mark = ' [sig]' if r['p_ols_hc3'] < 0.05 else ''
        title = (f"[{r['pathway']}] {r['ratio']}\n"
                 f"log2FC={r['log2FC']:+.3f}  FC={r['fold_change_linear']:.2f}×  "
                 f"p_OLS={r['p_ols_hc3']:.4f}  p_Wilc={r['p_wilcoxon']:.4f} {sig_mark}\n"
                 f"{r['rationale']}")
        color = PATH_COLORS.get(r['pathway'], 'black')
        ax.set_title(title, fontsize=9.5, pad=8, color=color)
        ax.spines['top'].set_color(color)
        ax.spines['top'].set_linewidth(2.5)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
    fig.suptitle('修订后 DHA 通路最终图景 — 4 条独立氧化通路 + 5-LOX 对照 '
                 '(n=165, ANCOVA-adjusted for age/GA/year/batch)',
                 fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    fig_path = ROOT / 'results/figures/04_Metabolite_boxplots/lox_dha_and_autoox_v2.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {fig_path.relative_to(ROOT)}')

    # Single-metabolite ANCOVA for 17/4 HDoHE (重复 15c 一行式输出便于对照)
    print('\n=== 单分子 ANCOVA: 17/4-HDoHE (log2 raw) — 检查方向 ===')
    for k, v in [('17-HDoHE', EXTRA['17-HDoHE']), ('4-HDoHE', EXTRA['4-HDoHE'])]:
        y = np.log2(full.loc[v].astype(float).values)
        m = sm.OLS(y, X).fit(cov_type='HC3')
        beta = float(m.params['BMI_overweight'])
        p = float(m.pvalues['BMI_overweight'])
        ci_lo, ci_hi = [float(x) for x in m.conf_int().loc['BMI_overweight'].tolist()]
        _, p_wil = stats.mannwhitneyu(y[groups == '超重肥胖'], y[groups == '正常'],
                                       alternative='two-sided')
        print(f'  {k}: log2FC={beta:+.3f}  FC={2**beta:.2f}×  '
              f'CI=({ci_lo:+.3f},{ci_hi:+.3f})  '
              f'p_OLS={p:.4f}  p_Wilc={float(p_wil):.4f}')


if __name__ == '__main__':
    main()
