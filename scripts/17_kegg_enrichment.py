"""17. KEGG 通路富集分析 (ORA + GSEA-like).

数据现状 (实测):
  - 80 主轨 67 特征中 40 有 KEGG ID (59.7%)
  - 80 主轨 2 个候选 (12-HETrE, 20-HDoHE) 都无 KEGG ID
  - 50 探索轨 4 个候选只有 1 个 (5,6-DiHETrE C14772) 有 KEGG ID
  → 直接用 supplier 的 KEGG ID 做 ORA 几乎无功效

补救策略:
  1. 基于化学结构 + 已知 KEGG 通路定义, 把全部 73 特征手工映射到所属 KEGG 通路
     - hsa00590: Arachidonic acid metabolism (AA 派生物主战场)
     - hsa00591: Linoleic acid metabolism (LA 派生物)
     - hsa00592: alpha-Linolenic acid metabolism (ALA/DHA/EPA/DPA 派生物)
     - hsa04923: Regulation of lipolysis in adipocytes (内源大麻素 + 部分 PG)
     - hsa01040: Biosynthesis of unsaturated fatty acids (PUFA 前体)
  2. ORA: 超几何检验 (候选 vs 全集 vs 通路)
  3. GSEA-like: 全集按 t_limma 排序, 每通路做 permutation test (10000 次)
  4. Bubble plot

局限性 (论文 Methods 必须明示):
  - 氧化脂质 panel 高度集中于 hsa00590 (~85% 特征), KEGG 区分度天然较差
  - 推荐互校用 SOP §5.1 自定义 7 家族富集 (scripts/16) 提供机制颗粒度
  - 候选数 (n=2-4) 在 ORA 框架下几乎无统计功效, GSEA 是更适合的方法

输入:
  results/tables/ancova_main_80.csv
  results/tables/diff_candidates_80.csv
  results/tables/diff_candidates_50.csv
  data/02_preprocessed/metabolite_family_map.csv (用 short_label / family_main 辅助标注)

输出:
  results/tables/kegg_enrichment_ora.csv
  results/tables/kegg_enrichment_gsea.csv
  results/tables/kegg_pathway_membership.csv (每代谢物 → 所属通路的审计表)
  results/figures/05_Enrichment/kegg_enrichment_bubble.png
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

N_PERM = 10000
RNG = np.random.default_rng(42)

# === KEGG 通路定义 ===
# 基于 KEGG REST API + 化学结构归属手工编码 (Homo sapiens 通路)
# 每个 metabolite 可属多通路 (cross-membership allowed)
KEGG_PATHWAYS = {
    'hsa00590': {
        'name': 'Arachidonic acid metabolism',
        'short': 'AA metabolism',
        'members': [
            # Free AA
            'Arachidonic acid',
            # COX 通路: PG / TX / PGI 衍生物
            'Thromboxane B2',
            'Prostaglandin A2', 'Prostaglandin D2', 'Prostaglandin E2', 'Prostaglandin F2α',
            'Prostaglandin J2', 'δ-12-Prostaglandin J2',
            '11β-Prostaglandin E2', '15-Deoxy-δ-12,14-prostaglandin D2',
            '11β-13,14-Dihydro-15-keto prostaglandinF2α',
            '15-Keto prostaglandin F2α', '6,15-Diketo-13,14-dihydro-prostaglandin F1α',
            '13,14-Dihydro-15-keto prostaglandin E2', '13,14-Dihydro-15-keto-pgf2α',
            '15-Keto prostaglandin E2',
            '12S-Hydroxy-5Z,8E,10E-heptadecatrienoic acid',
            '1α,1β-Dihomo prostaglandin E2',
            # LOX 通路: HETE / LT
            '11-Hydroxy-5Z,8Z,11E,14Z-eicosatetraenoic acid',
            '15-Hydroxy-5Z,8Z,11Z,13E-eicosatetraenoic acid',
            '8-Hydroxy-5Z,9E,11Z,14Z-eicosatetraenoic acid',
            '12-Hydroxy-5Z,8Z,10E,14Z-eicosatetraenoic acid',
            'Leukotriene E4', 'Leukotriene B4',
            # CYP/sEH 通路: EpETrE / DiHETrE / 16-18-HETE-ω
            '16-Hydroxy-5Z,8Z,11Z,14Z-eicosatetraenoic acid',
            '18-Hydroxy-5Z,8Z,11Z,14Z-eicosatetraenoic acid',
            '8,9-Dpoxy-5Z,11Z,14Z-eicosatrienoic acid',
            '11,12-Epoxy-5Z,8Z,14Z-eicosatrienoic acid',
            '5,6-DiHydroxy-8Z,11Z,14Z-eicosatrienoic acid',
            '8,9-DiHydroxy-5Z,11Z,14Z-eicosatrienoic acid',
            '11,12-DiHydroxy-5Z,8Z,14Z-eicosatrienoic acid',
            '14,15-DiHydroxy-5Z,8Z,11Z-eicosatrienoic acid',
            # HETrE (AA-derived monohydroxy)
            '8-Hydroxy-9E,11Z,14Z-eicosatrienoic acid',
            '12-Hydroxy-8Z,10E,14Z-eicosatrienoic acid',
            '15-Hydroxy-8Z,11Z,13E-eicosatrienoic acid',
            '15-Hydroxy-11Z,13E-eicosadienoic acid',
        ],
    },
    'hsa00591': {
        'name': 'Linoleic acid metabolism',
        'short': 'LA metabolism',
        'members': [
            'Linoleic acid', 'Linoelaidic acid', 'Conjugated linoleic acids',
            'Linoleoyl ethanolamide',  # LEA (LA → AEA-类 / NAPE-PLD)
            # HODE / oxoODE / HpODE
            '9-Hydroxy-10E,12Z-octadecadienoic acid',
            '13-Hydroxy-9Z,11E-octadecadienoic acid',
            '9-Oxo-10E,12Z-octadecadienoic acid',
            '13-Oxo-9Z,11E-octadecadienoicacid',
            '9-Hydroperoxy-10E,12E-octadecadienoic acid',
            '13-Hydroperoxy-9Z,11E-octadecadienoic acid',
            # EpOME / DiHOME (LA-CYP-sEH)
            '9,10-Epoxy-12Z-octadecenoic acid',
            '12,13-Epoxy-9Z-octadecenoic acid',
            '9,10-DiHydroxy-12Z-octadecenoic acid',
            '12,13 -DiHydroxy-9Z-octadecenoic acid',
        ],
    },
    'hsa00592': {
        'name': 'alpha-Linolenic acid metabolism',
        'short': 'α-Linolenic / ω3 metabolism',
        'members': [
            'α-Linolenic acid',
            # EPA 及其衍生物
            '12-Hydroxy-5,8,10,14,17-eicosapentaenoic acid',
            '15-Hydroperoxy-5,8,11,14,17-eicosapentaenoic acid',
            '11-Hydroxy- 5Z,8Z,12E,14Z,17Z-eicosapentaenoic acid',
            'Prostaglandin F3α', 'Thromboxane B3',
            # DHA 衍生物 (HDoHE)
            '14-Hydroxy-4Z,7Z,10Z,12E,16Z,19Z-docosahexaenoic acid',
            '8-Hydroxy-4Z,6E,10Z,13Z,16Z,19Z-docosahexaenoic acid',
            '10-Hydroxy-4Z,7Z,11E,13Z,16Z,19Z-docosahexaenoic acid',
            '11-Hydroxy-4Z,7Z,9E,13Z,16Z,19Z-docosahexaenoic acid',
            '13-Hydroxy-4Z,7Z,10Z,14E,16Z,19Z-docosahexaenoic acid',
            '7-Hydroxy-4Z,8E,10Z,13Z,16Z,19Z-docosahexaenoic acid',
            '20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid',
            '16-Hydroxy-4Z,7Z,10Z,13Z,17E,19Z-docosahexaenoic acid',
            # DPA 衍生物
            '19(20)-Epoxy-4Z,7Z,10Z,13Z,16Z-docosapentaenoic acid',
            '19,20-DiHydroxy-4Z,7Z,10Z,13Z,16Z-docosapentaenoic acid',
            # PGE1 系列 (DGLA, KEGG 部分归 hsa00592, 部分归 hsa00590)
            'Prostaglandin E1', 'Prostaglandin F1α', 'Prostaglandin D1',
            '15-Keto prostaglandin E1', '13,14-Dihydro-15-keto Prostaglandin E1',
            '13,14-Dihydro-15-keto Prostaglandin F1α',
        ],
    },
    'hsa04923': {
        'name': 'Regulation of lipolysis in adipocytes',
        'short': 'Adipocyte lipolysis (signaling)',
        'members': [
            'Anandamide',  # AEA, 内源大麻素 → CB1R → 抑制脂解
            'Arachidonic acid',
            'Prostaglandin E2',
        ],
    },
    'hsa01040': {
        'name': 'Biosynthesis of unsaturated fatty acids',
        'short': 'PUFA biosynthesis',
        'members': [
            'Arachidonic acid', 'Linoleic acid', 'α-Linolenic acid',
            'Linoelaidic acid', 'Conjugated linoleic acids',
        ],
    },
}


def hypergeometric_test(k, M, K, N):
    """超几何检验 (单侧 over-representation):
    k = 通路中的候选数
    M = 候选总数 (drawn)
    K = 通路成员数 (success in population)
    N = 全集大小
    H0: 候选与通路独立 → P(X >= k)
    """
    if k == 0:
        return 1.0
    p = stats.hypergeom.sf(k - 1, N, K, M)
    return float(p)


def permutation_gsea(observed_mean_rank, family_size, total_n, n_perm=N_PERM):
    """与 16 脚本相同方法: 随机抽 family_size 个特征, 比较 |mean_rank - 中位 rank| 偏离度."""
    median_rank = (total_n + 1) / 2.0
    obs_dev = observed_mean_rank - median_rank
    perm_means = np.array([
        np.mean(RNG.choice(total_n, size=family_size, replace=False) + 1)
        for _ in range(n_perm)
    ])
    perm_devs = perm_means - median_rank
    p_two = (np.sum(np.abs(perm_devs) >= abs(obs_dev) - 1e-12) + 1) / (n_perm + 1)
    p_up = (np.sum(perm_devs >= obs_dev - 1e-12) + 1) / (n_perm + 1)
    p_down = (np.sum(perm_devs <= obs_dev + 1e-12) + 1) / (n_perm + 1)
    return p_two, p_up, p_down


def main():
    print('=== 17. KEGG 通路富集 (ORA + GSEA) ===\n')

    # 加载数据
    anc = pd.read_csv(TABLES / 'ancova_main_80.csv', encoding='utf-8-sig')
    cand80 = pd.read_csv(TABLES / 'diff_candidates_80.csv', encoding='utf-8-sig')
    fmap = pd.read_csv(PREP_DIR / 'metabolite_family_map.csv', encoding='utf-8-sig')

    bg_set = set(anc['Metabolite Name'])      # 67 个全集
    cand_set = set(cand80['Metabolite Name'])  # 2 个候选

    # === 通路成员归属审计 ===
    membership = []
    for m in anc['Metabolite Name']:
        in_pw = []
        for pid, pdef in KEGG_PATHWAYS.items():
            if m in pdef['members']:
                in_pw.append(pid)
        membership.append({
            'metabolite': m,
            'short_label': fmap.set_index('Metabolite Name').loc[m, 'short_label'] if m in fmap['Metabolite Name'].values else m,
            'family_main': fmap.set_index('Metabolite Name').loc[m, 'family_main'] if m in fmap['Metabolite Name'].values else '-',
            'in_kegg_pathways': '; '.join(in_pw) if in_pw else '(none)',
            'is_candidate_80': m in cand_set,
        })
    mb_df = pd.DataFrame(membership)
    mb_df.to_csv(TABLES / 'kegg_pathway_membership.csv', index=False, encoding='utf-8-sig')
    print(f'  ✓ {(TABLES / "kegg_pathway_membership.csv").relative_to(ROOT)}')

    n_unmapped = int((mb_df['in_kegg_pathways'] == '(none)').sum())
    print(f'  全 67 特征中: 已映射 {len(mb_df) - n_unmapped}, 未映射 {n_unmapped}')

    # === ORA (超几何检验) ===
    print(f'\n--- ORA: 候选 (n={len(cand_set)}) vs 全集 (n={len(bg_set)}) ---')
    bg_in_pw = {pid: [m for m in pdef['members'] if m in bg_set]
                for pid, pdef in KEGG_PATHWAYS.items()}
    cand_in_pw = {pid: [m for m in pdef['members'] if m in cand_set]
                  for pid, pdef in KEGG_PATHWAYS.items()}

    ora_rows = []
    for pid, pdef in KEGG_PATHWAYS.items():
        K = len(bg_in_pw[pid])   # 通路成员数 (在全集中)
        k = len(cand_in_pw[pid])  # 通路中候选数
        M = len(cand_set)        # 候选总数
        N = len(bg_set)          # 全集
        if K == 0:
            continue
        p_ora = hypergeometric_test(k, M, K, N)
        # Enrichment ratio
        expect = M * K / N
        er = k / expect if expect > 0 else np.nan
        ora_rows.append({
            'pathway_id': pid,
            'pathway_name': pdef['name'],
            'pathway_size_K': K,
            'candidates_in_pathway_k': k,
            'all_candidates_M': M,
            'background_N': N,
            'expected_k': round(expect, 3),
            'enrichment_ratio': round(er, 3),
            'p_hypergeometric': round(p_ora, 6),
            'candidates_hit': '; '.join(cand_in_pw[pid]) if k > 0 else '',
        })

    ora_df = pd.DataFrame(ora_rows)
    ora_df['q_bh'] = multipletests(ora_df['p_hypergeometric'], method='fdr_bh')[1].round(6)
    ora_df = ora_df.sort_values('p_hypergeometric').reset_index(drop=True)
    ora_df.to_csv(TABLES / 'kegg_enrichment_ora.csv', index=False, encoding='utf-8-sig')
    print(ora_df[['pathway_id', 'pathway_name', 'pathway_size_K',
                  'candidates_in_pathway_k', 'enrichment_ratio',
                  'p_hypergeometric', 'q_bh']].to_string(index=False))
    print(f'  ✓ {(TABLES / "kegg_enrichment_ora.csv").relative_to(ROOT)}')

    # === GSEA-like ===
    print(f'\n--- GSEA-like: 全 67 按 t_limma 排序 ---')
    rank_df = anc.sort_values('t_limma', ascending=False).reset_index(drop=True)
    rank_df['rank'] = rank_df.index + 1
    n_total = len(rank_df)
    rank_lookup = dict(zip(rank_df['Metabolite Name'], rank_df['rank']))
    t_lookup = dict(zip(rank_df['Metabolite Name'], rank_df['t_limma']))

    gsea_rows = []
    for pid, pdef in KEGG_PATHWAYS.items():
        members_in_bg = bg_in_pw[pid]
        K = len(members_in_bg)
        if K == 0:
            continue
        ranks = np.array([rank_lookup[m] for m in members_in_bg])
        ts = np.array([t_lookup[m] for m in members_in_bg])
        mean_rank = float(np.mean(ranks))
        mean_t = float(np.mean(ts))
        n_up = int((ts > 0).sum())

        # MWU vs 通路外特征
        out_ts = np.array([t_lookup[m] for m in bg_set if m not in set(members_in_bg)])
        try:
            _, p_mwu = stats.mannwhitneyu(ts, out_ts, alternative='two-sided')
            p_mwu = float(p_mwu)
        except ValueError:
            p_mwu = np.nan

        # Permutation
        p_perm_two, p_perm_up, p_perm_down = permutation_gsea(mean_rank, K, n_total)

        # Leading-edge: top-3 by |t|
        sub = pd.DataFrame({'m': members_in_bg, 't': ts}).sort_values('t', ascending=False)
        leading_up = '; '.join(f'{r["m"][:35]}(t={r["t"]:+.2f})' for _, r in sub.head(3).iterrows() if r['t'] > 0)
        leading_down = '; '.join(f'{r["m"][:35]}(t={r["t"]:+.2f})' for _, r in sub.tail(3).iloc[::-1].iterrows() if r['t'] < 0)

        gsea_rows.append({
            'pathway_id': pid,
            'pathway_name': pdef['name'],
            'pathway_size_K': K,
            'mean_rank': round(mean_rank, 2),
            'mean_t_limma': round(mean_t, 4),
            'n_up_t>0': n_up,
            'frac_up': round(n_up / K, 3),
            'p_mwu_two_sided': round(p_mwu, 6) if pd.notna(p_mwu) else np.nan,
            'p_perm_two_sided': round(p_perm_two, 6),
            'p_perm_up': round(p_perm_up, 6),
            'p_perm_down': round(p_perm_down, 6),
            'leading_edge_up': leading_up,
            'leading_edge_down': leading_down,
        })

    gsea_df = pd.DataFrame(gsea_rows)
    gsea_df['q_mwu_bh'] = multipletests(gsea_df['p_mwu_two_sided'].fillna(1.0),
                                         method='fdr_bh')[1].round(6)
    gsea_df['q_perm_bh'] = multipletests(gsea_df['p_perm_two_sided'], method='fdr_bh')[1].round(6)
    gsea_df = gsea_df.sort_values('mean_t_limma', ascending=False).reset_index(drop=True)
    gsea_df.to_csv(TABLES / 'kegg_enrichment_gsea.csv', index=False, encoding='utf-8-sig')
    print(gsea_df[['pathway_id', 'pathway_name', 'pathway_size_K',
                   'mean_t_limma', 'frac_up',
                   'p_mwu_two_sided', 'p_perm_two_sided',
                   'q_mwu_bh', 'q_perm_bh']].to_string(index=False))
    print(f'  ✓ {(TABLES / "kegg_enrichment_gsea.csv").relative_to(ROOT)}')

    # === Bubble plot (基于 GSEA, 因 ORA 候选太少几乎无信号) ===
    fig, ax = plt.subplots(figsize=(11, 4.5))
    plot_df = gsea_df.sort_values('mean_t_limma').reset_index(drop=True)
    y_pos = np.arange(len(plot_df))
    sizes = (plot_df['pathway_size_K'] ** 1.3) * 10
    neg_log_p = -np.log10(plot_df['p_perm_two_sided'].clip(1e-4))
    sc = ax.scatter(plot_df['mean_t_limma'], y_pos,
                    s=sizes, c=neg_log_p, cmap='YlOrRd',
                    edgecolors='black', linewidths=0.6, alpha=0.85,
                    vmin=0, vmax=max(2.5, neg_log_p.max()))
    for i, r in plot_df.iterrows():
        ax.text(-1.8, i, f'{r["pathway_id"]}\n{r["pathway_name"][:30]}',
                ha='right', va='center', fontsize=9.5, fontweight='bold')
        ax.text(r['mean_t_limma'] + 0.04, i,
                f' K={r["pathway_size_K"]}, ↑{r["frac_up"]:.0%}',
                ha='left', va='center', fontsize=8.5)

    ax.axvline(0, color='gray', lw=1.0, linestyle='--')
    ax.set_yticks([])
    ax.set_xlabel('Mean t_limma (per KEGG pathway)  →  positive = ↑ in overweight/obese', fontsize=10.5)
    ax.set_title('KEGG pathway enrichment (GSEA-like, 80% main track, 10k perms)\n'
                 'bubble size ∝ pathway size; color = -log10 permutation p',
                 fontsize=11.5, fontweight='bold', pad=10)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_xlim(-2.5, 1.5)
    ax.set_ylim(-0.7, len(plot_df) - 0.3)

    cb = plt.colorbar(sc, ax=ax, pad=0.02, shrink=0.85)
    cb.set_label('-log10 (permutation p, two-sided)', fontsize=9.5)

    plt.subplots_adjust(left=0.02, right=0.96, top=0.88, bottom=0.13)
    fig_path = FIG_DIR / 'kegg_enrichment_bubble.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'\n  ✓ {fig_path.relative_to(ROOT)}')

    print('\n=== 局限性 (论文 Methods 必须写) ===')
    print('  - 氧化脂质 panel 高度集中于 hsa00590 (~85% 特征), KEGG 区分度天然差')
    print('  - ORA 在候选 n=2 下几乎无统计功效 → 主要参考 GSEA-like 结果')
    print('  - 推荐配合 scripts/16 自定义 7 家族富集 (基于酶系归类) 互校')


if __name__ == '__main__':
    main()
