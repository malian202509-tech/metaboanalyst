"""20. 本地 MSEA QEA (GlobalTest) 通路富集分析.

Python 自实现 MetaboAnalyst MSEA QEA 的底层算法 (Goeman et al. 2004,
Bioinformatics), 不依赖 R/MetaboAnalystR/MetaboAnalyst web.

  输入: results/metaboanalyst/concentration_residualized_80.csv
        (165 样本 × 67 代谢物, 已 ComBat + 协变量残差化, 标签 Normal/Overweight_Obese)
  通路库: 脚本 17 KEGG_PATHWAYS (5 条手工策划的 KEGG 通路, hsa00590/591/592/
        04923/01040)
  算法: GlobalTest Q 统计量 + 10000 次 permutation; BH-FDR 校正

GlobalTest QEA 数学:
  对每条通路 P (含 K 个代谢物):
    X (n×K) = 该通路成员的浓度子矩阵
    Y (n,)  = BMI 二分类 (Overweight_Obese=1, Normal=0), 均值 μ̂
    T (K,)  = X' (Y - μ̂) / sqrt(μ̂(1-μ̂))  — 标准化后的 X-Y 协方差向量
    Q_obs   = ||T||² / K                  — 平均每代谢物的 score² (二次型)
    permutation 打乱 Y 10000 次, 重算 Q → null distribution
    p_perm  = (#{Q_perm ≥ Q_obs} + 1) / (n_perm + 1)

  与 Subramanian GSEA 的关键区别:
    - GSEA: ranked list + ES + KS-like 检验, 假设代谢物独立
    - QEA:  浓度矩阵 + 二次型统计, 自动吸收代谢物间的协方差结构

输出 (results/tables/):
  msea_qea_local_80.csv           5 通路 × QEA 结果 (residualized 输入, 主结果)
  msea_qea_local_80_raw.csv       5 通路 × QEA 结果 (未残差化 ComBat 输入, 对照)
  msea_qea_compare_residvsraw.md  残差化 vs 未残差化对比小结
输出 (results/figures/05_Enrichment/):
  msea_qea_bubble_local.png       bubble plot (仿 MetaboAnalyst 风格)
"""
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))

# 直接 import 脚本 17 的 KEGG_PATHWAYS + 脚本 19 的 asciify_name
# (避免重复定义; importlib 不会执行 __name__ == '__main__' 块)
import importlib.util


def _load_module(modname, path):
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s17 = _load_module('s17', ROOT / 'scripts' / '17_kegg_enrichment.py')
s19 = _load_module('s19', ROOT / 'scripts' / '19_metaboanalyst_export.py')

# 浓度表列名已 ASCII 化 (脚本 19 的修复), 因此通路成员匹配也需 ASCII 化
asciify_name = s19.asciify_name
KEGG_PATHWAYS = {
    pid: {**pdef, 'members': [asciify_name(m) for m in pdef['members']]}
    for pid, pdef in s17.KEGG_PATHWAYS.items()
}

MA_DIR = ROOT / 'results' / 'metaboanalyst'
TABLES = ROOT / 'results' / 'tables'
FIG_DIR = ROOT / 'results' / 'figures' / '05_Enrichment'
FIG_DIR.mkdir(parents=True, exist_ok=True)

N_PERM = 10000
RNG = np.random.default_rng(42)


def globaltest_qea(X, Y, n_perm=N_PERM, rng=None):
    """GlobalTest QEA 统计量 + permutation p.

    X: (n, K) 浓度子矩阵, 一行一样本, 一列一代谢物
    Y: (n,)   二分类标签 (0/1)
    return: dict 含 Q_obs, p_perm, T_vec (per-metabolite contribution)
    """
    if rng is None:
        rng = np.random.default_rng(42)
    n, K = X.shape
    Y = np.asarray(Y).astype(float)
    mu = Y.mean()
    sd = np.sqrt(mu * (1 - mu))
    if sd < 1e-12:
        return {'Q_obs': np.nan, 'p_perm': np.nan, 'T_vec': np.full(K, np.nan)}

    # 中心化 Y, 然后 T_k = sum_i X_ik (Y_i - μ̂) / sd
    Y_c = (Y - mu) / sd
    # 也对 X 中心化 (移除每个代谢物的均值, 增加数值稳定; 不影响 Q 大小排序)
    Xc = X - X.mean(axis=0, keepdims=True)
    T_obs = Xc.T @ Y_c          # (K,)
    Q_obs = float(np.dot(T_obs, T_obs) / K)

    # Permutation: 打乱 Y_c (在 H0 下与 X 独立)
    Q_perm = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        Y_perm = rng.permutation(Y_c)
        T_perm = Xc.T @ Y_perm
        Q_perm[i] = np.dot(T_perm, T_perm) / K

    p_perm = (np.sum(Q_perm >= Q_obs - 1e-12) + 1) / (n_perm + 1)
    return {
        'Q_obs': Q_obs,
        'p_perm': float(p_perm),
        'T_vec': T_obs,
        'Q_perm_median': float(np.median(Q_perm)),
        'Q_perm_99pct': float(np.quantile(Q_perm, 0.99)),
    }


def per_metabolite_direction(X, Y, mu):
    """每个代谢物在 BMI 方向上的 effect (用于 leading-edge 标注).

    简化版: pearson r (X 列 vs Y) — 因为已残差化, 等价于 BMI 净效应方向.
    return: (K,) 的 r 数组
    """
    Y_c = Y - Y.mean()
    Xc = X - X.mean(axis=0)
    num = Xc.T @ Y_c
    denom = np.sqrt((Xc ** 2).sum(axis=0) * (Y_c ** 2).sum())
    return num / np.where(denom > 0, denom, 1.0)


def run_qea(conc_path, label_name):
    """对一个浓度表跑 5 通路 QEA, 返回 DataFrame."""
    df = pd.read_csv(conc_path, encoding='utf-8-sig')
    print(f'\n--- {label_name} ---')
    print(f'  输入: {conc_path.relative_to(ROOT)}  shape={df.shape}')

    # Sample / Label / 67 metabolites
    Y = (df['Label'] == 'Overweight_Obese').astype(int).values
    bg_metabolites = [c for c in df.columns if c not in ('Sample', 'Label')]
    X_full = df[bg_metabolites].values
    print(f'  样本: {len(Y)} (Overweight_Obese={Y.sum()}, Normal={len(Y)-Y.sum()}); '
          f'代谢物全集: {len(bg_metabolites)}')

    rows = []
    for pid, pdef in KEGG_PATHWAYS.items():
        # 通路成员 ∩ 全集
        in_bg = [m for m in pdef['members'] if m in bg_metabolites]
        K = len(in_bg)
        if K < 2:
            print(f'  跳过 {pid} ({pdef["short"]}): 只有 {K} 个成员在全集中')
            continue
        idx = [bg_metabolites.index(m) for m in in_bg]
        X_path = X_full[:, idx]

        rng = np.random.default_rng(42)
        res = globaltest_qea(X_path, Y, n_perm=N_PERM, rng=rng)

        # Leading edge: top-3 by |r|
        r_per_met = per_metabolite_direction(X_path, Y, Y.mean())
        order = np.argsort(-np.abs(r_per_met))[:3]
        leading = '; '.join(
            f'{in_bg[i][:36]}(r={r_per_met[i]:+.3f})' for i in order
        )

        rows.append({
            'pathway_id': pid,
            'pathway_name': pdef['name'],
            'pathway_size_K': K,
            'Q_obs': round(res['Q_obs'], 4),
            'Q_perm_median': round(res['Q_perm_median'], 4),
            'Q_perm_99pct': round(res['Q_perm_99pct'], 4),
            'p_perm': round(res['p_perm'], 6),
            'n_up (r>0)': int((r_per_met > 0).sum()),
            'frac_up': round(float((r_per_met > 0).mean()), 3),
            'leading_top3_by_|r|': leading,
        })

    out = pd.DataFrame(rows)
    out['q_bh'] = multipletests(out['p_perm'], method='fdr_bh')[1].round(6)
    out = out.sort_values('p_perm').reset_index(drop=True)
    return out


def bubble_plot(df, fig_path, title_suffix=''):
    """仿 MetaboAnalyst QEA bubble: x=enrichment ratio (Q_obs/Q_perm_median), y=pathway, size=K, color=-log10 p."""
    df = df.copy()
    df['enrichment_ratio'] = df['Q_obs'] / df['Q_perm_median'].replace(0, np.nan)
    df = df.sort_values('enrichment_ratio', ascending=True).reset_index(drop=True)
    n = len(df)

    fig, ax = plt.subplots(figsize=(10, max(3.5, n * 0.9)))
    sizes = (df['pathway_size_K'] ** 1.3) * 18
    neg_log_p = -np.log10(df['p_perm'].clip(1 / (N_PERM + 1)))
    sc = ax.scatter(df['enrichment_ratio'], np.arange(n),
                    s=sizes, c=neg_log_p, cmap='YlOrRd',
                    edgecolors='black', linewidths=0.7, alpha=0.9,
                    vmin=0, vmax=max(2.0, neg_log_p.max()))

    for i, r in df.iterrows():
        ax.text(0.02, i, f'{r["pathway_id"]}\n{r["pathway_name"][:35]}',
                ha='right', va='center', fontsize=9.5, fontweight='bold',
                transform=ax.get_yaxis_transform())
        ax.text(r['enrichment_ratio'] + 0.04, i,
                f' K={r["pathway_size_K"]}, p={r["p_perm"]:.3g}',
                ha='left', va='center', fontsize=8.8)

    ax.axvline(1.0, color='gray', lw=1.0, linestyle='--', alpha=0.7)
    ax.set_yticks([])
    ax.set_xlabel('Enrichment ratio = Q_obs / median(Q_perm)   (>1 → 富集)', fontsize=10.5)
    ax.set_title(f'MSEA QEA (local GlobalTest, 10k perms){title_suffix}\n'
                 f'bubble size ∝ pathway size; color = -log10 p',
                 fontsize=11.5, fontweight='bold', pad=10)
    ax.grid(axis='x', alpha=0.3, linestyle='--')

    xmax = max(df['enrichment_ratio'].max() * 1.4, 2.0)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.7, n - 0.3)

    cb = plt.colorbar(sc, ax=ax, pad=0.02, shrink=0.85)
    cb.set_label('-log10 (permutation p)', fontsize=9.5)

    plt.subplots_adjust(left=0.3, right=0.92, top=0.88, bottom=0.13)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    print('=== 20. 本地 MSEA QEA (GlobalTest, 10000 perms) ===')
    print(f'通路库: {len(KEGG_PATHWAYS)} 条 KEGG 手工通路 (来自脚本 17)')

    # 主结果: residualized
    resid_path = MA_DIR / 'concentration_residualized_80.csv'
    raw_path = MA_DIR / 'concentration_raw_80.csv'

    out_resid = run_qea(resid_path, 'residualized 80 (协变量已洗, 主结果)')
    out_raw = run_qea(raw_path, 'raw ComBat 80 (未残差化, 对照)')

    # 输出 CSV
    out_resid_path = TABLES / 'msea_qea_local_80.csv'
    out_raw_path = TABLES / 'msea_qea_local_80_raw.csv'
    out_resid.to_csv(out_resid_path, index=False, encoding='utf-8-sig')
    out_raw.to_csv(out_raw_path, index=False, encoding='utf-8-sig')
    print(f'\n✓ {out_resid_path.relative_to(ROOT)}')
    print(f'✓ {out_raw_path.relative_to(ROOT)}')

    # 打印主结果
    print('\n=== 主结果 (residualized): MSEA QEA 5 通路富集 ===')
    print(out_resid[['pathway_id', 'pathway_name', 'pathway_size_K',
                     'Q_obs', 'p_perm', 'q_bh', 'frac_up']].to_string(index=False))

    # Bubble plot (residualized)
    fig_path = FIG_DIR / 'msea_qea_bubble_local.png'
    bubble_plot(out_resid, fig_path, title_suffix=' — residualized (covariate-adjusted)')
    print(f'\n✓ {fig_path.relative_to(ROOT)}')

    # 对比小结
    md_lines = [
        '# MSEA QEA 本地结果对比: residualized vs raw\n',
        '本地 GlobalTest QEA, n=165, K-pathways=5 (KEGG hsa 手工策划), 10000 perms.\n',
        '## residualized (协变量 age + GA + year 已洗)\n',
        out_resid[['pathway_id', 'pathway_name', 'pathway_size_K',
                   'Q_obs', 'p_perm', 'q_bh', 'frac_up']].to_markdown(index=False),
        '\n## raw (未残差化 ComBat 后矩阵)\n',
        out_raw[['pathway_id', 'pathway_name', 'pathway_size_K',
                 'Q_obs', 'p_perm', 'q_bh', 'frac_up']].to_markdown(index=False),
        '\n## 解读',
        '- residualized 是主结果 (与脚本 09 ANCOVA / 脚本 16 家族 GSEA 框架一致)',
        '- raw 用作对照, 看协变量调整是否改变富集 p 值显著性顺序',
        '- 与脚本 17 KEGG GSEA-rank (rank-based) 互校: 若两者方向一致, 信号更稳健',
    ]
    md_path = TABLES / 'msea_qea_compare_residvsraw.md'
    md_path.write_text('\n'.join(md_lines), encoding='utf-8')
    print(f'✓ {md_path.relative_to(ROOT)}')

    print('\n=== 完成 ===')


if __name__ == '__main__':
    main()
