"""17. ORA 候选富集分析 (Over-Representation Analysis, 超几何检验).

设计:
  - 用候选 (差异代谢物) vs 全集 vs 各 gene set, 做超几何检验
  - 三个 gene set 维度互补:
      (a) family_main (7 大家族, 与 scripts/16 同源)
      (b) substrate (上游底物: AA / LA / DHA / EPA / DGLA / DPA / ALA / 内源大麻素前体)
      (c) enzyme_system (下游酶: COX / LOX / CYP / sEH / NAPE-PLD / Non-enzymatic)
  - 双轨各跑一次: 80 主轨 (候选 2) + 50 探索轨 (候选 4)

数学事实 (函数说明里也声明):
  - 80 轨候选 n=2: 命中 1 家族的 hypergeom p ∈ [0.10, 0.30]
  - 50 轨候选 n=4: 命中 1 家族的 hypergeom p ∈ [0.20, 0.50]
  - BH-q<0.05 在该候选数下几乎不可能, 但 ORA 能告诉:
      "哪些家族被命中" + "命中是否高于随机期望 (enrichment ratio)"
  - 与 GSEA-like (scripts/16) 互补: GSEA 看家族整体方向, ORA 看候选集中度

输入:
  results/tables/ancova_main_{80,50}.csv (全集)
  results/tables/diff_candidates_{80,50}.csv (候选)
  data/02_preprocessed/metabolite_family_map.csv (family_main / substrate / enzyme)

输出:
  results/tables/ora_candidate_enrichment.csv (双轨 × 三维度合一表)
  results/figures/05_Enrichment/ora_candidate_bubble.png (按 enrichment_ratio + p 排布)
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / 'results' / 'tables'
PREP_DIR = ROOT / 'data' / '02_preprocessed'
FIG_DIR = ROOT / 'results' / 'figures' / '05_Enrichment'
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 三个 gene set 维度
DIMENSIONS = ['family_main', 'substrate', 'enzyme']


def hypergeom_ora(k, K, M, N):
    """超几何 over-representation:
      k = 该 set 中候选数
      K = 该 set 大小 (全集中)
      M = 候选总数
      N = 全集大小
    H0: 候选与 set 独立 → P(X >= k)
    """
    if k == 0:
        return 1.0
    return float(stats.hypergeom.sf(k - 1, N, K, M))


def run_ora_for_track(tag, fmap):
    """对一个轨道跑三个维度的 ORA, 返回长表 (每行 = 一个 set × 一个维度)."""
    anc = pd.read_csv(TABLES / f'ancova_main_{tag}.csv', encoding='utf-8-sig')
    cand = pd.read_csv(TABLES / f'diff_candidates_{tag}.csv', encoding='utf-8-sig')

    bg_set = set(anc['Metabolite Name'])
    cand_set = set(cand['Metabolite Name'])
    N = len(bg_set)
    M = len(cand_set)

    fmap_in_bg = fmap[fmap['Metabolite Name'].isin(bg_set)].copy()

    print(f'\n=== 轨道 {tag}%: N={N}, M={M} 候选 ===')

    rows = []
    for dim in DIMENSIONS:
        gsets = sorted([v for v in fmap_in_bg[dim].dropna().unique() if v != '-'])
        for gs in gsets:
            members = set(fmap_in_bg.loc[fmap_in_bg[dim] == gs, 'Metabolite Name'])
            K = len(members)
            hits = members & cand_set
            k = len(hits)
            p = hypergeom_ora(k, K, M, N)
            expect = M * K / N if N > 0 else 0
            er = (k / expect) if expect > 0 else np.nan
            rows.append({
                'track': tag,
                'dimension': dim,
                'gene_set': gs,
                'set_size_K': K,
                'candidates_in_set_k': k,
                'total_candidates_M': M,
                'background_N': N,
                'expected_k': round(expect, 3),
                'enrichment_ratio': round(er, 3) if pd.notna(er) else np.nan,
                'p_hypergeom': round(p, 6),
                'candidates_hit': '; '.join(sorted(hits)),
            })
    return pd.DataFrame(rows)


def main():
    print('=== 17. ORA 候选富集 (超几何) ===')

    fmap = pd.read_csv(PREP_DIR / 'metabolite_family_map.csv', encoding='utf-8-sig')
    print(f'  家族映射: {len(fmap)} 个特征; 3 维度 (family_main / substrate / enzyme)')

    all_rows = []
    for tag in ['80', '50']:
        sub = run_ora_for_track(tag, fmap)
        # BH 校正 (per track + per dimension, 这样 family_main 7 个 q 独立校正)
        for dim in DIMENSIONS:
            mask = sub['dimension'] == dim
            pvals = sub.loc[mask, 'p_hypergeom'].values
            sub.loc[mask, 'q_bh'] = np.round(multipletests(pvals, method='fdr_bh')[1], 6)
        sub = sub.sort_values(['dimension', 'p_hypergeom']).reset_index(drop=True)
        all_rows.append(sub)

        # 打印 hit 的部分
        hit_only = sub[sub['candidates_in_set_k'] > 0]
        if len(hit_only) > 0:
            print(f'  --- 命中的 set ({tag}%) ---')
            print(hit_only[['dimension', 'gene_set', 'set_size_K',
                            'candidates_in_set_k', 'expected_k',
                            'enrichment_ratio', 'p_hypergeom', 'q_bh',
                            'candidates_hit']].to_string(index=False))

    res = pd.concat(all_rows, ignore_index=True)
    out_csv = TABLES / 'ora_candidate_enrichment.csv'
    res.to_csv(out_csv, index=False, encoding='utf-8-sig')
    print(f'\n  ✓ {out_csv.relative_to(ROOT)}')

    # === Bubble plot: 只画命中的 set, 双轨并列 ===
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax_idx, tag in enumerate(['80', '50']):
        ax = axes[ax_idx]
        sub = res[(res['track'] == tag) & (res['candidates_in_set_k'] > 0)].copy()
        if len(sub) == 0:
            ax.text(0.5, 0.5, f'No hits in {tag}% track', ha='center', va='center',
                    transform=ax.transAxes)
            continue
        sub = sub.sort_values(['dimension', 'enrichment_ratio'], ascending=[True, True]).reset_index(drop=True)

        # 显示标签: dim + gene_set
        labels = [f'[{r["dimension"][:3]}] {r["gene_set"]}' for _, r in sub.iterrows()]
        y_pos = np.arange(len(sub))

        sizes = (sub['set_size_K'] ** 1.0) * 30
        neg_log_p = -np.log10(sub['p_hypergeom'].clip(1e-4))
        sc = ax.scatter(sub['enrichment_ratio'], y_pos,
                        s=sizes, c=neg_log_p, cmap='YlOrRd',
                        edgecolors='black', linewidths=0.6, alpha=0.85,
                        vmin=0, vmax=max(1.5, neg_log_p.max()))
        for i, r in sub.iterrows():
            ax.text(r['enrichment_ratio'] + 0.1, i,
                    f' k={r["candidates_in_set_k"]}/K={r["set_size_K"]}, '
                    f'p={r["p_hypergeom"]:.3f}',
                    va='center', fontsize=8.5)
        ax.axvline(1.0, color='gray', lw=1.0, linestyle='--')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8.5)
        ax.set_xlabel('Enrichment ratio (observed / expected)', fontsize=10)
        ax.set_title(f'{tag}% track  (M={int(sub["total_candidates_M"].iloc[0])}, '
                     f'N={int(sub["background_N"].iloc[0])})',
                     fontsize=11, fontweight='bold')
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        cb = plt.colorbar(sc, ax=ax, pad=0.02, shrink=0.7)
        cb.set_label('-log10 p', fontsize=9)

    fig.suptitle('ORA candidate enrichment (hypergeometric) — hits only',
                 fontsize=13, fontweight='bold', y=0.99)
    plt.subplots_adjust(left=0.18, right=0.96, top=0.92, bottom=0.10, wspace=0.55)
    fig_path = FIG_DIR / 'ora_candidate_bubble.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {fig_path.relative_to(ROOT)}')

    # 简要解读
    print('\n=== 解读提示 ===')
    print('  - enrichment_ratio > 1 = 候选在该 set 中过表达 (即使 p 不显著也有方向意义)')
    print('  - p_hypergeom 在候选 n=2-4 下数学上很难 <0.05 (统计功效问题), 不是数据问题')
    print('  - 与 scripts/16 GSEA-like 互校: GSEA 看家族整体 t 方向, ORA 看候选命中集中度')


if __name__ == '__main__':
    main()
