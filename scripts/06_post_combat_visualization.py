"""06. ComBat 校正后可视化与诊断: PCA + PLS-DA + OPLS-DA.

每条数据轨生成一张 2x3 大图:
  (0,0) 校正前 PCA 按年份   - 应该看到三色抱团
  (0,1) 校正后 PCA 按年份   - 应该混匀  ← ComBat 是否生效的判定
  (0,2) 校正后 PCA 按 BMI   - 论文 Figure 1 候选 (BMI 在全局可能仍不强分离, 正常)
  (1,0) PLS-DA t1 vs t2 按 BMI                - 监督学习视角
  (1,1) OPLS-DA t_pred vs t_ortho 按 BMI      - 监督学习 + 旋转
  (1,2) PLS-DA permutation test (R²Y 真实 vs 200 次随机) - 验证是否过拟合

也定量重算 PERMANOVA 年份效应, 看是否从 p=0.001 降到不显著.

输入:
  data/02_preprocessed/ori_n165_filtered{80,50}_log2_pareto.xlsx              (校正前)
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat_pareto.xlsx   (校正后)
  data/02_preprocessed/sample_alignment_n165.csv

输出:
  results/figures/01_QC_and_Batch_PCA/PostCombat_diagnosis_{80,50}.png
  results/tables/post_combat_audit.csv
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Ellipse
from pathlib import Path
from sklearn.decomposition import PCA
from scipy.stats import f as f_dist

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
from utils.opls import opls_da, pls_da, permutation_test_opls

PRE_DIR = ROOT / 'data' / '02_preprocessed'
POST_DIR = ROOT / 'data' / '03_batch_corrected'
FIG_DIR = ROOT / 'results' / 'figures' / '01_QC_and_Batch_PCA'
TBL_DIR = ROOT / 'results' / 'tables'
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'

mpl.rcParams['font.family'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
mpl.rcParams['axes.unicode_minus'] = False

YEAR_COLORS = {'A19 (2019)': '#1f77b4', 'B20 (2020)': '#ff7f0e', 'C21 (2021)': '#2ca02c'}
BMI_COLORS = {'正常': '#4C72B0', '超重肥胖': '#DD8452'}


def load_matrix(path):
    df = pd.read_excel(path)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    return df[sample_cols].values.T, sample_cols  # X (n_samp × n_feat), sample_ids


def permanova(X, group, n_perm=999, random_state=42):
    rng = np.random.default_rng(random_state)
    g = pd.Series(group).reset_index(drop=True)
    levels = g.unique()
    n = len(g)
    a = len(levels)
    diff = X[:, None, :] - X[None, :, :]
    D2 = np.sum(diff ** 2, axis=2)

    def F(perm_g):
        ss_t = D2.sum() / (2 * n)
        ss_w = 0.0
        for lv in levels:
            idx = np.where(perm_g == lv)[0]
            if len(idx) > 1:
                ss_w += D2[np.ix_(idx, idx)].sum() / (2 * len(idx))
        ss_a = ss_t - ss_w
        if ss_w <= 0:
            return np.inf
        return (ss_a / (a - 1)) / (ss_w / (n - a))

    F_obs = F(g.values)
    cnt = sum(F(rng.permutation(g.values)) >= F_obs for _ in range(n_perm))
    return float(F_obs), float((cnt + 1) / (n_perm + 1))


def hotelling_ellipse(x, y, conf=0.95):
    """Hotelling's T² 95% 置信椭圆参数 (二维).
    返回 (center, width, height, angle_deg) 或 None (样本数 <5).
    SIMCA-P 风格: T²_crit = p(n-1)/(n-p) × F_{p, n-p, conf}.
    """
    pts = np.column_stack([x, y])
    n = len(pts)
    if n < 5:
        return None
    mu = pts.mean(axis=0)
    cov = np.cov(pts.T, ddof=1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # 按特征值降序
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    p = 2
    t2_crit = (p * (n - 1) / (n - p)) * f_dist.ppf(conf, p, n - p)
    semi_axes = np.sqrt(np.maximum(eigvals, 0) * t2_crit)
    width = 2 * semi_axes[0]
    height = 2 * semi_axes[1]
    angle = float(np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0])))
    return mu, width, height, angle


def plot_combined_permutation(ax, perm, n_perm, title='OPLS-DA Permutation'):
    """MetaboAnalyst 风格 R²Y + Q² 合并置换直方图.
    向下箭头标实测值, 文字标注 R²Y/Q²/P."""
    all_vals = np.concatenate([perm['r2y_perm'], perm['q2_perm'],
                               [perm['r2y_obs'], perm['q2_obs']]])
    x_min = float(min(all_vals.min(), -0.1)) - 0.03
    x_max = float(max(all_vals.max(), 0.3)) + 0.03
    bins = np.linspace(x_min, x_max, 60)

    counts_r2y, _, _ = ax.hist(perm['r2y_perm'], bins=bins,
                                color='#a8c8e0', alpha=0.7, edgecolor='none',
                                label=f'Perm R²Y')
    counts_q2, _, _ = ax.hist(perm['q2_perm'], bins=bins,
                               color='#f4b6a5', alpha=0.7, edgecolor='none',
                               label=f'Perm Q²')

    max_freq = float(max(counts_r2y.max(), counts_q2.max())) if max(counts_r2y.max(), counts_q2.max()) > 0 else 10
    # 箭头高度: R²Y 较低 (0.55), Q² 较高 (0.85) 错开避免重叠
    arrow_r = max_freq * 0.55
    arrow_q = max_freq * 0.85

    # R²Y 黑色向下箭头
    ax.annotate('', xy=(perm['r2y_obs'], 0),
                xytext=(perm['r2y_obs'], arrow_r),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.5, mutation_scale=15))
    ax.text(perm['r2y_obs'], arrow_r + max_freq * 0.02,
            f'R²Y: {perm["r2y_obs"]:.3f}\np = {perm["p_r2y"]:.3f}',
            fontsize=9, ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.85, edgecolor='gray', lw=0.5))

    # Q² 黑色向下箭头
    ax.annotate('', xy=(perm['q2_obs'], 0),
                xytext=(perm['q2_obs'], arrow_q),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.5, mutation_scale=15))
    ax.text(perm['q2_obs'], arrow_q + max_freq * 0.02,
            f'Q²: {perm["q2_obs"]:.3f}\np = {perm["p_q2"]:.3f}',
            fontsize=9, ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.85, edgecolor='gray', lw=0.5))

    ax.set_xlabel('Permutations')
    ax.set_ylabel('Frequency')
    ax.legend(loc='center right', fontsize=8)
    # 判定: Q² > 0 且 p_Q² < 0.05 是真正的通过标准
    passed = (perm['p_r2y'] < 0.05 and perm['p_q2'] < 0.05 and perm['q2_obs'] > 0)
    verdict = '通过验证' if passed else '过拟合 / 模型无效'
    title_color = 'green' if passed else 'red'
    ax.set_title(f'{title} ({n_perm} 次置换)  →  {verdict}',
                 fontsize=10, color=title_color)
    ax.grid(axis='y', alpha=0.3)


def scatter_group(ax, x, y, groups, color_map, x_label, y_label, title, add_ellipse=True):
    for lv, color in color_map.items():
        m = (groups == lv)
        if m.sum() == 0:
            continue
        # Hotelling 95% 置信椭圆 (虚线 + 半透明填充)
        if add_ellipse:
            res = hotelling_ellipse(x[m], y[m])
            if res is not None:
                mu, w, h, ang = res
                ax.add_patch(Ellipse(xy=mu, width=w, height=h, angle=ang,
                                     facecolor=color, alpha=0.12,
                                     edgecolor=color, linestyle='--', linewidth=1.2))
        ax.scatter(x[m], y[m], c=color, alpha=0.7, edgecolor='white',
                   linewidth=0.5, s=45, label=f'{lv} (n={int(m.sum())})')
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title, fontsize=10)
    ax.legend(loc='best', fontsize=8)


def diagnose_one(tier, tag, tag_short, pre_path, post_path):
    print(f'\n{"="*60}')
    print(f'轨道 {tier}')
    X_pre, sample_ids = load_matrix(pre_path)
    X_post, _ = load_matrix(post_path)
    mp = pd.read_csv(PRE_DIR / 'sample_alignment_n165.csv')
    meta = pd.DataFrame({'omx_id': sample_ids}).merge(
        mp[['omx_id', 'BMI_group']], on='omx_id', how='left'
    )
    meta['year'] = meta['omx_id'].str[:3].map({
        'A19': 'A19 (2019)', 'B20': 'B20 (2020)', 'C21': 'C21 (2021)'
    })
    n_samp = X_pre.shape[0]

    # PCA pre/post
    pca_pre = PCA(n_components=5, random_state=42).fit(X_pre)
    pcs_pre = pca_pre.transform(X_pre)
    var_pre = pca_pre.explained_variance_ratio_ * 100

    pca_post = PCA(n_components=5, random_state=42).fit(X_post)
    pcs_post = pca_post.transform(X_post)
    var_post = pca_post.explained_variance_ratio_ * 100

    print(f'  PCA pre  PC1/PC2: {var_pre[0]:.1f}% / {var_pre[1]:.1f}%')
    print(f'  PCA post PC1/PC2: {var_post[0]:.1f}% / {var_post[1]:.1f}%')

    # PERMANOVA pre vs post (年份)
    print(f'  --- PERMANOVA (年份) ---')
    F_pre, p_pre = permanova(X_pre, meta['year'].values, n_perm=999)
    F_post, p_post = permanova(X_post, meta['year'].values, n_perm=999)
    print(f'    校正前: F = {F_pre:.3f}  p = {p_pre:.3f}')
    print(f'    校正后: F = {F_post:.3f}  p = {p_post:.3f}')
    print(f'    PERMANOVA (BMI) 校正后:', end=' ')
    F_bmi, p_bmi = permanova(X_post, meta['BMI_group'].values, n_perm=999)
    print(f'F = {F_bmi:.3f}  p = {p_bmi:.3f}')

    # PLS-DA (仅 scores 图用)
    y_bin = (meta['BMI_group'] == '超重肥胖').astype(int).values
    print(f'  --- PLS-DA (2 comp) ---')
    pls = pls_da(X_post, y_bin, n_components=2)
    print(f'    R²Y = {pls["r2y"]:.3f}   R²X per comp = {[round(r, 3) for r in pls["r2x_per_comp"]]}')

    # OPLS-DA (用于 scores 图)
    print(f'  --- OPLS-DA (1 predictive + 1 orthogonal) ---')
    op = opls_da(X_post, y_bin, n_ortho=1)
    print(f'    R²Y = {op["r2y"]:.3f}  R²X_pred = {op["r2x_pred"]:.3f}  R²X_ortho = {op["r2x_ortho"]:.3f}')

    # OPLS-DA 置换检验: R²Y + Q² 双指标 (MetaboAnalyst 风格合并图)
    print(f'  --- OPLS-DA 置换检验 (200 次, 7-fold CV) ---')
    opls_perm = permutation_test_opls(X_post, y_bin, n_ortho=1, n_perm=200,
                                       n_folds=7, random_state=42)
    print(f'    R²Y_obs = {opls_perm["r2y_obs"]:.3f}   p_R²Y = {opls_perm["p_r2y"]:.3f}')
    print(f'    Q²_obs  = {opls_perm["q2_obs"]:.3f}   p_Q²  = {opls_perm["p_q2"]:.3f}')

    # === 2x3 大图 ===
    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    fig.suptitle(f'ComBat 校正后诊断与可视化 — {tier}', fontsize=14, fontweight='bold')

    # (0,0) 校正前 PCA 按年份
    scatter_group(ax[0, 0], pcs_pre[:, 0], pcs_pre[:, 1], meta['year'].values, YEAR_COLORS,
                  f'PC1 ({var_pre[0]:.1f}%)', f'PC2 ({var_pre[1]:.1f}%)',
                  f'校正前 PCA · 年份  (PERMANOVA p={p_pre:.3f})')

    # (0,1) 校正后 PCA 按年份
    scatter_group(ax[0, 1], pcs_post[:, 0], pcs_post[:, 1], meta['year'].values, YEAR_COLORS,
                  f'PC1 ({var_post[0]:.1f}%)', f'PC2 ({var_post[1]:.1f}%)',
                  f'校正后 PCA · 年份  (PERMANOVA p={p_post:.3f})  ← 应混匀')

    # (0,2) 校正后 PCA 按 BMI
    scatter_group(ax[0, 2], pcs_post[:, 0], pcs_post[:, 1], meta['BMI_group'].values, BMI_COLORS,
                  f'PC1 ({var_post[0]:.1f}%)', f'PC2 ({var_post[1]:.1f}%)',
                  f'校正后 PCA · BMI  (PERMANOVA p={p_bmi:.3f})')

    # (1,0) PLS-DA
    scatter_group(ax[1, 0], pls['t'][:, 0], pls['t'][:, 1], meta['BMI_group'].values, BMI_COLORS,
                  f't1 (R²X={pls["r2x_per_comp"][0]:.2f})',
                  f't2 (R²X={pls["r2x_per_comp"][1]:.2f})',
                  f'PLS-DA scores  R²Y={pls["r2y"]:.3f}')

    # (1,1) OPLS-DA
    scatter_group(ax[1, 1], op['t_pred'], op['t_ortho'][:, 0], meta['BMI_group'].values, BMI_COLORS,
                  f't_predictive (R²X={op["r2x_pred"]:.2f})',
                  f't_orthogonal (R²X={op["r2x_ortho"]:.2f})',
                  f'OPLS-DA scores  R²Y={op["r2y"]:.3f}, Q²={opls_perm["q2_obs"]:.3f}')

    # (1,2) MetaboAnalyst 风格 OPLS-DA 合并置换直方图 (R²Y + Q²)
    plot_combined_permutation(ax[1, 2], opls_perm, n_perm=200,
                              title='OPLS-DA R²Y + Q² 置换分布')

    plt.tight_layout()
    out_fig = FIG_DIR / f'PostCombat_diagnosis_{tag_short}.png'
    plt.savefig(out_fig, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'  图已保存: {out_fig.relative_to(ROOT)}')

    return {
        'tier': tier,
        'PERMANOVA年份_校正前_F': round(F_pre, 3),
        'PERMANOVA年份_校正前_p': round(p_pre, 4),
        'PERMANOVA年份_校正后_F': round(F_post, 3),
        'PERMANOVA年份_校正后_p': round(p_post, 4),
        'PERMANOVA_BMI_校正后_F': round(F_bmi, 3),
        'PERMANOVA_BMI_校正后_p': round(p_bmi, 4),
        'PLS-DA_R2Y': round(pls['r2y'], 3),
        'OPLS-DA_R2Y': round(op['r2y'], 3),
        'OPLS-DA_R2X_pred': round(op['r2x_pred'], 3),
        'OPLS-DA_R2X_ortho': round(op['r2x_ortho'], 3),
        'OPLS-DA_Q2_CV': round(opls_perm['q2_obs'], 3),
        'OPLS-DA_p_R2Y': round(opls_perm['p_r2y'], 4),
        'OPLS-DA_p_Q2': round(opls_perm['p_q2'], 4),
    }


def main():
    print('=== 06. ComBat 校正后可视化与诊断 (双轨并行) ===')
    records = []
    records.append(diagnose_one(
        '80% 主分析', '80%主', '80main',
        PRE_DIR / 'ori_n165_filtered80_log2_pareto.xlsx',
        POST_DIR / 'ori_n165_filtered80_log2_combat_pareto.xlsx',
    ))
    records.append(diagnose_one(
        '50% 探索性', '50%探', '50expl',
        PRE_DIR / 'ori_n165_filtered50_log2_pareto.xlsx',
        POST_DIR / 'ori_n165_filtered50_log2_combat_pareto.xlsx',
    ))
    audit = pd.DataFrame(records)
    audit_path = TBL_DIR / 'post_combat_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')
    print()
    print(audit.to_string(index=False))


if __name__ == '__main__':
    main()
