"""15e. AA 池通路全景 — panel 内所有 AA 衍生物 ANCOVA 结果汇总.

按 7 条 AA 通路组织 panel 内每个分子的方向 + p + log2FC:
  1. AA 母体 (Free PUFA)
  2. AA → COX 通路 (PG / TX / HHT)
  3. AA → LOX 通路 (HETE / HETrE 经 ETrE 中间体)
  4. AA → CYP-Epo (EpETrE, 环氧化)
  5. AA → CYP-Epo → sEH (DiHETrE, 二醇)
  6. AA → CYP4 ω-羟化 (16/18-HETE)
  7. AA → 非酶氧化 (11/8-HETE, racemic)
  8. AA → 内源大麻素 (AEA, NAPE-PLD)

数据源:
  ancova_main_80.csv (67 特征)
  ancova_main_50.csv (73 特征, 含 LTB4/LTE4 等 80 主轨没收的)
  metabolite_family_map.csv (substrate 与 enzyme 注释)

输出:
  results/tables/aa_pool_pathway_landscape.csv
  results/figures/04_Metabolite_boxplots/aa_pool_forest.png
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
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

# AA 池通路分组 (顺序固定: COX → LOX → CYP-Epo → sEH → CYP4 → 非酶 → 内源大麻素 → 母体)
PATHWAY_ORDER = [
    'AA 母体', 'COX (PG/TX)', 'LOX (HETE/HETrE/LT)',
    'CYP-Epo (EpETrE)', 'sEH (DiHETrE)',
    'CYP4 ω-羟化 (HETE-ω)', '非酶氧化 (auto-ox HETE)',
    '内源大麻素 (AEA)',
]
PATH_COLORS = {
    'AA 母体':             '#666666',
    'COX (PG/TX)':         '#D62728',
    'LOX (HETE/HETrE/LT)': '#2CA02C',
    'CYP-Epo (EpETrE)':    '#1F77B4',
    'sEH (DiHETrE)':       '#9467BD',
    'CYP4 ω-羟化 (HETE-ω)': '#17BECF',
    '非酶氧化 (auto-ox HETE)': '#FF7F0E',
    '内源大麻素 (AEA)':      '#bcbd22',
}


def classify_pathway(row):
    """Classify each AA-pool metabolite into one of 8 pathways."""
    fmain = row['family_main']
    fsub = row['family_sub']
    enz = row['enzyme']
    sub = row['substrate']
    sl = row['short_label']
    # 1. Free PUFA = AA 母体
    if fmain == 'Free PUFA' and sl == 'AA':
        return 'AA 母体'
    # 2. COX 通路
    if fmain == 'AA-COX':
        return 'COX (PG/TX)'
    # 3. LOX (HETE/HETrE + LT)
    if fmain == 'AA-LOX':
        if sl in ['11-HETE', '8-HETE']:
            return '非酶氧化 (auto-ox HETE)'
        return 'LOX (HETE/HETrE/LT)'
    # 4-6. CYP/sEH 家族细分
    if fmain == 'AA-CYP/sEH':
        if 'EpETrE' in fsub:
            return 'CYP-Epo (EpETrE)'
        if 'DiHETrE' in fsub:
            return 'sEH (DiHETrE)'
        if 'HETE-ω' in fsub:
            return 'CYP4 ω-羟化 (HETE-ω)'
    # 8. AEA
    if fmain == 'Endocannabinoid' and sl == 'AEA':
        return '内源大麻素 (AEA)'
    return None


def main():
    print('=== 15e. AA 池通路全景 ===\n')

    fmap = pd.read_csv(ROOT / 'data/02_preprocessed/metabolite_family_map.csv')
    anc80 = pd.read_csv(ROOT / 'results/tables/ancova_main_80.csv')
    anc50 = pd.read_csv(ROOT / 'results/tables/ancova_main_50.csv')

    # 标记 AA 池成员 (substrate = AA 或 AA→ETrE; 加上 AEA endo; 加上母体 AA)
    aa_pool = fmap[
        (fmap['substrate'].isin(['AA', 'AA→ETrE'])) |
        (fmap['short_label'] == 'AEA')
    ].copy()
    aa_pool['pathway'] = aa_pool.apply(classify_pathway, axis=1)
    unclass = aa_pool[aa_pool['pathway'].isna()]
    if len(unclass) > 0:
        print(f'  ⚠ {len(unclass)} 个 AA 池分子未分类:')
        for _, r in unclass.iterrows():
            sl = r['short_label']
            fm = r['family_main']
            fs = r['family_sub']
            print(f'    {sl:<20} fmain={fm} fsub={fs}')

    # 合并 ANCOVA 结果 (优先 80, 缺失补 50)
    anc80_sub = anc80[['Metabolite Name', 'log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi',
                       'p_ols_hc3', 'p_limma', 'p_wilcoxon', 'is_robust_to_outlier']].copy()
    anc80_sub['source_track'] = '80_main'
    anc50_sub = anc50[['Metabolite Name', 'log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi',
                       'p_ols_hc3', 'p_limma', 'p_wilcoxon', 'is_robust_to_outlier']].copy()
    anc50_sub['source_track'] = '50_explor'
    in80 = set(anc80['Metabolite Name'])
    extra50 = anc50_sub[~anc50_sub['Metabolite Name'].isin(in80)]
    anc_all = pd.concat([anc80_sub, extra50], ignore_index=True)

    merged = aa_pool.merge(anc_all, on='Metabolite Name', how='left')
    # AA 母体 (Free PUFA) 不在 ancova_main 表 (Free PUFA 家族不参与 limma 主分析)
    # 从 free_pufa_singlemolecule.csv 拿
    pufa = pd.read_csv(ROOT / 'results/tables/free_pufa_singlemolecule.csv')
    aa_row = pufa[pufa['PUFA'] == 'AA'].iloc[0]
    aa_idx = merged[merged['short_label'] == 'AA'].index
    if len(aa_idx) > 0:
        i = aa_idx[0]
        merged.loc[i, 'log2FC'] = aa_row['log2FC']
        merged.loc[i, 'log2FC_CI_lo'] = aa_row['CI_lo']
        merged.loc[i, 'log2FC_CI_hi'] = aa_row['CI_hi']
        merged.loc[i, 'p_ols_hc3'] = aa_row['p_ols_hc3']
        merged.loc[i, 'p_wilcoxon'] = aa_row['p_wilcoxon']
        merged.loc[i, 'source_track'] = 'free_pufa'

    miss = merged['log2FC'].isna().sum()
    if miss > 0:
        print(f'  ⚠ {miss} 个分子缺 ANCOVA 数据:')
        print(merged[merged['log2FC'].isna()][['short_label', 'pathway']])

    # Sort by pathway order, then by log2FC
    merged['_path_ord'] = merged['pathway'].map(
        {p: i for i, p in enumerate(PATHWAY_ORDER)})
    merged = merged.sort_values(['_path_ord', 'log2FC'],
                                 ascending=[True, False]).reset_index(drop=True)

    out_cols = ['short_label', 'pathway', 'family_main', 'family_sub', 'enzyme',
                'enzyme_evidence_level', 'non_enzymatic_contribution',
                'in_main_80', 'log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi',
                'p_ols_hc3', 'p_limma', 'p_wilcoxon', 'is_robust_to_outlier',
                'source_track']
    out = merged[out_cols].copy()
    out.to_csv(ROOT / 'results/tables/aa_pool_pathway_landscape.csv',
               index=False, encoding='utf-8-sig')
    print(f'  ✓ results/tables/aa_pool_pathway_landscape.csv  (n={len(out)} 分子)')

    # 按通路汇总
    print('\n  === 按通路汇总 ===')
    print(f'{"通路":<28} {"n":>3}  {"↑ raw p<0.05":>14} {"↓ raw p<0.05":>14} '
          f'{"↑ 趋势(全体)":>14}')
    print('-' * 90)
    for path in PATHWAY_ORDER:
        sub = merged[merged['pathway'] == path]
        if len(sub) == 0:
            continue
        n = len(sub)
        sig_up = ((sub['p_ols_hc3'] < 0.05) & (sub['log2FC'] > 0)).sum()
        sig_dn = ((sub['p_ols_hc3'] < 0.05) & (sub['log2FC'] < 0)).sum()
        all_up = (sub['log2FC'] > 0).sum()
        print(f'{path:<28} {n:>3}  {sig_up:>14} {sig_dn:>14} '
              f'{f"{all_up}/{n}":>14}')

    # ===== Forest plot =====
    fig, ax = plt.subplots(figsize=(12, max(8, 0.32 * len(merged))))
    ys = np.arange(len(merged))[::-1]
    for i, r in merged.iterrows():
        if pd.isna(r['log2FC']):
            continue
        color = PATH_COLORS.get(r['pathway'], 'black')
        y = ys[i]
        ax.errorbar(r['log2FC'], y,
                    xerr=[[r['log2FC'] - r['log2FC_CI_lo']],
                          [r['log2FC_CI_hi'] - r['log2FC']]],
                    fmt='o', color=color, markersize=7, capsize=3, linewidth=1.5)
        sig = ' *' if r['p_ols_hc3'] < 0.05 else ''
        trk = ' [50]' if r['source_track'] == '50_explor' else ''
        ax.text(r['log2FC_CI_hi'] + 0.04, y,
                f"{r['short_label']}  p={r['p_ols_hc3']:.3f}{sig}{trk}",
                va='center', fontsize=8, color=color)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_yticks(ys)
    ax.set_yticklabels([r['short_label'] for _, r in merged.iterrows()], fontsize=8)
    ax.set_xlabel('log2 fold-change (超重肥胖 vs 正常)  ±95% CI', fontsize=10)
    ax.set_title('AA 池通路全景 — panel 内所有 AA 衍生物 ANCOVA 结果\n'
                 '(* = raw p<0.05; [50] = 50 探索轨; n=165, 调整 age/GA/year/batch)',
                 fontsize=11.5, pad=10)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    handles = [mpatches.Patch(color=PATH_COLORS[p], label=p) for p in PATHWAY_ORDER]
    ax.legend(handles=handles, loc='lower right', fontsize=8.5, frameon=True, ncol=1)
    ax.set_xlim(-1.8, 2.0)
    plt.tight_layout()
    out_fig = ROOT / 'results/figures/04_Metabolite_boxplots/aa_pool_forest.png'
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'\n  ✓ {out_fig.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
