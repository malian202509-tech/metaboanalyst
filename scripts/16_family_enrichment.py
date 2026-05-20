"""16. GSEA-like 家族富集 (SOP §5.3 改造版).

候选数太少 (80轨 2, 50轨 4), 传统 ORA 不可行 → 用 ranked-list 富集:
  - 把 67 个特征按 t_limma 排序 (t 大 = 超重肥胖↑, t 小 = ↓)
  - 对每个家族 (7 类) 做两种互补检验:
      (a) Mann-Whitney U: 家族内 t vs 家族外 t (one-sided ↑ 与 ↓ 各 one-sided)
      (b) Permutation: 随机抽相同大小的特征集 10000 次, 比较 mean_rank 偏离度
  - leading-edge: 家族内 top-3 t 最大或最小的成员

输出:
  results/tables/family_enrichment_ranked.csv
  results/figures/05_Enrichment/family_enrichment_bubble.png
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

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

N_PERMUTATIONS = 10000
RNG = np.random.default_rng(42)

FAMILY_ORDER = [
    'Free PUFA', 'Endocannabinoid', 'LA-oxylipin',
    'AA-COX', 'AA-LOX', 'AA-CYP/sEH', 'ω-3 PUFA oxylipins',
]
FAMILY_COLORS = {
    'Free PUFA':            '#A0A0A0',
    'Endocannabinoid':      '#9E9AC8',
    'LA-oxylipin':          '#74C476',
    'AA-COX':               '#FB6A4A',
    'AA-LOX':               '#FD8D3C',
    'AA-CYP/sEH':           '#6BAED6',
    'ω-3 PUFA oxylipins': '#41AB5D',
}


def permutation_p(observed_mean_rank, family_size, total_n, n_perm=N_PERMUTATIONS):
    """随机抽 family_size 个特征的 mean rank, 比较 |observed - 中位 rank| 是否极端 (双侧).

    Returns: (p_two_sided, p_up_directional, p_down_directional)
    """
    median_rank = (total_n + 1) / 2.0
    obs_dev = observed_mean_rank - median_rank   # 正 = 偏向高 rank (高 t)

    perm_means = np.zeros(n_perm)
    for i in range(n_perm):
        idx = RNG.choice(total_n, size=family_size, replace=False)
        perm_means[i] = np.mean(idx + 1)  # rank 从 1 起算

    perm_devs = perm_means - median_rank
    p_two = (np.sum(np.abs(perm_devs) >= abs(obs_dev) - 1e-12) + 1) / (n_perm + 1)
    p_up = (np.sum(perm_devs >= obs_dev - 1e-12) + 1) / (n_perm + 1)
    p_down = (np.sum(perm_devs <= obs_dev + 1e-12) + 1) / (n_perm + 1)
    return p_two, p_up, p_down


def main():
    print('=== 16. GSEA-like 家族富集 ===\n')

    # 加载 80 主轨 ANCOVA + 家族映射
    anc = pd.read_csv(TABLES / 'ancova_main_80.csv', encoding='utf-8-sig')
    fmap = pd.read_csv(PREP_DIR / 'metabolite_family_map.csv', encoding='utf-8-sig')
    df = anc.merge(fmap[['Metabolite Name', 'short_label', 'family_main']],
                    on='Metabolite Name', how='inner')
    df = df.sort_values('t_limma', ascending=False).reset_index(drop=True)
    df['rank'] = df.index + 1  # 1 = 最高 t (最 ↑)
    n_total = len(df)
    print(f'  排序 {n_total} 个特征 (按 t_limma 降序; rank=1 → 最 ↑)')

    rows = []
    for fam in FAMILY_ORDER:
        sub = df[df['family_main'] == fam]
        n = len(sub)
        if n == 0:
            continue

        ranks_in = sub['rank'].values
        t_in = sub['t_limma'].values
        t_out = df.loc[~df['family_main'].eq(fam), 't_limma'].values

        mean_rank = float(np.mean(ranks_in))
        mean_t = float(np.mean(t_in))
        median_t = float(np.median(t_in))
        n_up = int((t_in > 0).sum())
        n_down = int((t_in < 0).sum())

        # Mann-Whitney (双侧)
        try:
            _, p_mwu = stats.mannwhitneyu(t_in, t_out, alternative='two-sided')
            p_mwu = float(p_mwu)
        except ValueError:
            p_mwu = np.nan

        # Mann-Whitney 单侧 (上侧)
        try:
            _, p_mwu_up = stats.mannwhitneyu(t_in, t_out, alternative='greater')
            p_mwu_up = float(p_mwu_up)
            _, p_mwu_down = stats.mannwhitneyu(t_in, t_out, alternative='less')
            p_mwu_down = float(p_mwu_down)
        except ValueError:
            p_mwu_up = p_mwu_down = np.nan

        # Permutation
        p_perm_two, p_perm_up, p_perm_down = permutation_p(mean_rank, n, n_total)

        # Leading-edge: top-3 + bottom-3 by t (内部)
        sub_sorted = sub.sort_values('t_limma', ascending=False)
        top3 = sub_sorted.head(3)
        leading_up = '; '.join(
            f'{r["short_label"]}(t={r["t_limma"]:+.2f})' for _, r in top3.iterrows()
            if r['t_limma'] > 0
        )
        bot3 = sub_sorted.tail(3).iloc[::-1]
        leading_down = '; '.join(
            f'{r["short_label"]}(t={r["t_limma"]:+.2f})' for _, r in bot3.iterrows()
            if r['t_limma'] < 0
        )

        rows.append({
            'family':              fam,
            'n_members':           n,
            'mean_rank':           round(mean_rank, 2),
            'median_rank':         float(np.median(ranks_in)),
            'mean_t_limma':        round(mean_t, 4),
            'median_t_limma':      round(median_t, 4),
            'n_up_t>0':            n_up,
            'n_down_t<0':          n_down,
            'frac_up':             round(n_up / n, 3),
            'p_mwu_two_sided':     round(p_mwu, 6) if pd.notna(p_mwu) else np.nan,
            'p_mwu_up':            round(p_mwu_up, 6) if pd.notna(p_mwu_up) else np.nan,
            'p_mwu_down':          round(p_mwu_down, 6) if pd.notna(p_mwu_down) else np.nan,
            'p_perm_two_sided':    round(p_perm_two, 6),
            'p_perm_up':           round(p_perm_up, 6),
            'p_perm_down':         round(p_perm_down, 6),
            'leading_edge_up':     leading_up,
            'leading_edge_down':   leading_down,
        })

    res = pd.DataFrame(rows)

    # BH 校正 (两个主要 p: MWU 双侧 + Perm 双侧)
    from statsmodels.stats.multitest import multipletests
    for col_p, col_q in [('p_mwu_two_sided', 'q_mwu_bh'),
                          ('p_perm_two_sided', 'q_perm_bh')]:
        res[col_q] = multipletests(res[col_p].fillna(1.0), method='fdr_bh')[1].round(6)

    # 打印
    print('\n  === 家族富集结果 (按 mean_t_limma 降序) ===')
    print(res.sort_values('mean_t_limma', ascending=False)[
        ['family', 'n_members', 'frac_up', 'mean_t_limma',
         'p_mwu_two_sided', 'p_perm_two_sided', 'q_mwu_bh', 'q_perm_bh']
    ].to_string(index=False))

    out_csv = TABLES / 'family_enrichment_ranked.csv'
    res.to_csv(out_csv, index=False, encoding='utf-8-sig')
    print(f'\n  ✓ {out_csv.relative_to(ROOT)}')

    # === Bubble plot ===
    fig, ax = plt.subplots(figsize=(11, 6.5))

    # 按 mean_t 排序 (上至下 = ↑ 至 ↓)
    plot_df = res.sort_values('mean_t_limma').reset_index(drop=True)
    y_pos = np.arange(len(plot_df))

    # bubble size ∝ n_members; color = -log10 p_perm_two
    sizes = (plot_df['n_members'] ** 1.5) * 12
    neg_log_p = -np.log10(plot_df['p_perm_two_sided'].clip(1e-4))
    sc = ax.scatter(plot_df['mean_t_limma'], y_pos,
                    s=sizes, c=neg_log_p, cmap='YlOrRd',
                    edgecolors='black', linewidths=0.6, alpha=0.85,
                    vmin=0, vmax=max(2.5, neg_log_p.max()))

    # 标记每个家族颜色块
    for i, r in plot_df.iterrows():
        ax.text(-2.0, i, r['family'], ha='right', va='center',
                fontsize=10.5, fontweight='bold',
                color=FAMILY_COLORS.get(r['family'], 'black'))
        # 标注 n_members + frac_up
        ax.text(r['mean_t_limma'] + 0.05, i, f' n={r["n_members"]}, ↑{r["frac_up"]:.0%}',
                ha='left', va='center', fontsize=8.5)

    ax.axvline(0, color='gray', lw=1.0, linestyle='--')
    ax.set_yticks([])
    ax.set_xlabel('Mean t_limma (per family)  →  positive = ↑ in overweight/obese', fontsize=11)
    ax.set_title('Family-level enrichment: ranked-list of 67 oxylipins by t_limma\n'
                 '(bubble size ∝ n_members; color = -log10 permutation p, 10k perms)',
                 fontsize=12, fontweight='bold', pad=10)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_xlim(-2.8, 1.5)
    ax.set_ylim(-0.8, len(plot_df) - 0.2)

    cb = plt.colorbar(sc, ax=ax, pad=0.02, shrink=0.85)
    cb.set_label('-log10 (permutation p, two-sided)', fontsize=10)

    plt.subplots_adjust(left=0.02, right=0.96, top=0.90, bottom=0.10)
    fig_path = FIG_DIR / 'family_enrichment_bubble.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {fig_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
