"""04. 批次效应诊断 (按采样年份).

样本 ID 前缀对应年份: A19 -> 2019, B20 -> 2020, C21 -> 2021

诊断流程 (两条数据轨各跑一遍):
  Step 1. PCA (Pareto-scaled log2 数据)
  Step 2. 视觉: PC1-PC2 散点 (按年份/按 BMI_group), Scree plot, PC1 密度
  Step 3. 定量:
    a) 每个 PC 按年份和 BMI 分别做 ANOVA, 取 R² 和 P
    b) PERMANOVA 在欧氏距离矩阵上检验年份效应
    c) 年份 x BMI Cramér's V (共线性复核)
  Step 4. 输出图 + 审计 CSV

输入:
  data/02_preprocessed/ori_n165_filtered{80,50}_log2_pareto.xlsx
  data/02_preprocessed/sample_alignment_n165.csv

输出:
  results/figures/01_QC_and_Batch_PCA/PCA_diagnosis_{80,50}.png
  results/tables/batch_diagnosis_audit.csv
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from sklearn.decomposition import PCA
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
PRE_DIR = ROOT / 'data' / '02_preprocessed'
FIG_DIR = ROOT / 'results' / 'figures' / '01_QC_and_Batch_PCA'
TBL_DIR = ROOT / 'results' / 'tables'
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'

# 字体: Windows 中文
mpl.rcParams['font.family'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
mpl.rcParams['axes.unicode_minus'] = False

YEAR_COLORS = {'A19 (2019)': '#1f77b4', 'B20 (2020)': '#ff7f0e', 'C21 (2021)': '#2ca02c'}
BMI_COLORS = {'正常': '#4C72B0', '超重肥胖': '#DD8452'}


def cramers_v(x, y):
    """Cramér's V for two categorical variables."""
    ct = pd.crosstab(x, y)
    chi2 = stats.chi2_contingency(ct).statistic
    n = ct.values.sum()
    r, k = ct.shape
    return float(np.sqrt(chi2 / (n * (min(r, k) - 1))))


def permanova(X, group, n_perm=999, random_state=42):
    """PERMANOVA (Anderson 2001) on Euclidean distances.
    返回 (pseudo-F, p-value)."""
    rng = np.random.default_rng(random_state)
    g = pd.Series(group).reset_index(drop=True)
    levels = g.unique()
    n = len(g)
    a = len(levels)

    # 距离矩阵 (平方欧氏)
    diff = X[:, None, :] - X[None, :, :]
    D2 = np.sum(diff ** 2, axis=2)

    def pseudo_F(perm_g):
        ss_t = D2.sum() / (2 * n)
        ss_w = 0.0
        for lv in levels:
            idx = np.where(perm_g == lv)[0]
            if len(idx) > 1:
                ss_w += D2[np.ix_(idx, idx)].sum() / (2 * len(idx))
        ss_a = ss_t - ss_w
        denom_a = a - 1
        denom_w = n - a
        if ss_w <= 0 or denom_w == 0:
            return np.inf
        return (ss_a / denom_a) / (ss_w / denom_w)

    F_obs = pseudo_F(g.values)
    cnt = 0
    for _ in range(n_perm):
        perm = rng.permutation(g.values)
        if pseudo_F(perm) >= F_obs:
            cnt += 1
    p = (cnt + 1) / (n_perm + 1)
    return float(F_obs), float(p)


def anova_pc_by_group(pc, group):
    """对单个 PC 用 group 做一元 ANOVA, 返回 R² 和 P."""
    df = pd.DataFrame({'pc': pc, 'g': group})
    grand_mean = df['pc'].mean()
    ss_t = ((df['pc'] - grand_mean) ** 2).sum()
    ss_b = sum(len(sub) * (sub['pc'].mean() - grand_mean) ** 2
               for _, sub in df.groupby('g'))
    ss_w = ss_t - ss_b
    levels = df['g'].nunique()
    n = len(df)
    if ss_w <= 0 or n - levels <= 0 or levels < 2:
        return np.nan, np.nan
    F = (ss_b / (levels - 1)) / (ss_w / (n - levels))
    p = 1 - stats.f.cdf(F, levels - 1, n - levels)
    r2 = ss_b / ss_t if ss_t > 0 else np.nan
    return float(r2), float(p)


def diagnose_one(src_name, tag, tag_short):
    print(f'\n{"="*60}')
    print(f'轨道 {tag}: {src_name}')
    df = pd.read_excel(PRE_DIR / src_name)
    mp = pd.read_csv(PRE_DIR / 'sample_alignment_n165.csv')
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]

    # 矩阵: 行=样本, 列=特征 (PCA 标准形式)
    X = df[sample_cols].values.T  # shape: (n_samp, n_feat)
    n_samp, n_feat = X.shape
    print(f'  矩阵形状 (样本x特征): {X.shape}')

    # 样本元数据
    meta = pd.DataFrame({'omx_id': sample_cols})
    meta = meta.merge(mp[['omx_id', 'BMI_group']], on='omx_id', how='left')
    meta['year_prefix'] = meta['omx_id'].str[:3]
    meta['year'] = meta['year_prefix'].map({'A19': 'A19 (2019)', 'B20': 'B20 (2020)', 'C21': 'C21 (2021)'})

    # PCA
    n_components = min(10, n_samp - 1, n_feat)
    pca = PCA(n_components=n_components, random_state=42)
    pcs = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100
    print(f'  PC1: {var[0]:.1f}%  PC2: {var[1]:.1f}%  PC3: {var[2]:.1f}% '
          f'(累积前 5: {var[:5].sum():.1f}%)')

    # === 定量: 每个 PC 的 R² ===
    quant = []
    for k in range(min(5, n_components)):
        r2_y, p_y = anova_pc_by_group(pcs[:, k], meta['year'])
        r2_b, p_b = anova_pc_by_group(pcs[:, k], meta['BMI_group'])
        quant.append({
            'tier': tag,
            'PC': f'PC{k+1}',
            '解释方差%': round(var[k], 2),
            '年份R²': round(r2_y, 3),
            '年份P': f'{p_y:.3g}' if not np.isnan(p_y) else 'NA',
            'BMIR²': round(r2_b, 3),
            'BMI P': f'{p_b:.3g}' if not np.isnan(p_b) else 'NA',
        })
    quant_df = pd.DataFrame(quant)
    print('  --- PC 方差归因 ---')
    print(quant_df.to_string(index=False))

    # PERMANOVA (年份)
    print('  --- PERMANOVA (Pareto 距离, 999 置换) ---')
    F_year, p_year = permanova(X, meta['year'].values, n_perm=999)
    F_bmi, p_bmi = permanova(X, meta['BMI_group'].values, n_perm=999)
    print(f'    年份: pseudo-F = {F_year:.3f}  p = {p_year:.3f}')
    print(f'    BMI:  pseudo-F = {F_bmi:.3f}  p = {p_bmi:.3f}')

    # 共线性
    cv = cramers_v(meta['year'], meta['BMI_group'])
    print(f"  --- 共线性 ---  年份 x BMI Cramér's V = {cv:.3f}  "
          f"({'低' if cv<0.2 else '中' if cv<0.5 else '高'})")

    # === 可视化 (2x2) ===
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(f'批次效应诊断 — {tag}  ({n_samp} 样本 × {n_feat} 代谢物, log2+Pareto)',
                 fontsize=14, fontweight='bold')

    # (0,0) PCA 按年份
    ax = axes[0, 0]
    for lv, color in YEAR_COLORS.items():
        m = meta['year'] == lv
        ax.scatter(pcs[m, 0], pcs[m, 1], c=color, label=f'{lv} (n={int(m.sum())})',
                   alpha=0.7, edgecolor='white', linewidth=0.5, s=50)
    ax.set_xlabel(f'PC1 ({var[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({var[1]:.1f}%)')
    ax.set_title(f"PC1-PC2 按年份  (年份R² on PC1={quant[0]['年份R²']}, P={quant[0]['年份P']})")
    ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
    ax.legend(loc='best', fontsize=9)

    # (0,1) PCA 按 BMI
    ax = axes[0, 1]
    for lv, color in BMI_COLORS.items():
        m = meta['BMI_group'] == lv
        ax.scatter(pcs[m, 0], pcs[m, 1], c=color, label=f'{lv} (n={int(m.sum())})',
                   alpha=0.7, edgecolor='white', linewidth=0.5, s=50)
    ax.set_xlabel(f'PC1 ({var[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({var[1]:.1f}%)')
    ax.set_title(f"PC1-PC2 按 BMI 分组  (BMI R² on PC1={quant[0]['BMIR²']}, P={quant[0]['BMI P']})")
    ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
    ax.legend(loc='best', fontsize=9)

    # (1,0) Scree plot
    ax = axes[1, 0]
    n_show = min(10, len(var))
    ax.bar(range(1, n_show + 1), var[:n_show], color='#5b8ec2', edgecolor='black')
    ax.plot(range(1, n_show + 1), np.cumsum(var[:n_show]), 'o-', color='red', label='累积 %')
    ax.set_xlabel('主成分')
    ax.set_ylabel('解释方差 %')
    ax.set_title(f'Scree plot  (前 5 累积 {var[:5].sum():.1f}%)')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # (1,1) PC1 密度图
    ax = axes[1, 1]
    for lv, color in YEAR_COLORS.items():
        m = meta['year'] == lv
        if m.sum() > 1:
            ax.hist(pcs[m, 0], bins=20, alpha=0.45, color=color, label=lv, density=True)
    ax.set_xlabel(f'PC1 ({var[0]:.1f}%)')
    ax.set_ylabel('密度')
    ax.set_title(f"PC1 分布按年份  (PERMANOVA p={p_year:.3f}, Cramér's V={cv:.2f})")
    ax.legend(fontsize=9)

    plt.tight_layout()
    out_fig = FIG_DIR / f'PCA_diagnosis_{tag_short}.png'
    plt.savefig(out_fig, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'  图已保存: {out_fig.relative_to(ROOT)}')

    return quant + [{
        'tier': tag, 'PC': 'PERMANOVA-year', '解释方差%': '',
        '年份R²': round(F_year, 3), '年份P': f'{p_year:.3f}',
        'BMIR²': round(F_bmi, 3), 'BMI P': f'{p_bmi:.3f}',
    }, {
        'tier': tag, 'PC': "Cramer's V (年份×BMI)", '解释方差%': '',
        '年份R²': round(cv, 3), '年份P': '', 'BMIR²': '', 'BMI P': '',
    }]


def main():
    print('=== 04. 批次效应诊断 (双轨并行) ===')
    all_records = []
    all_records += diagnose_one('ori_n165_filtered80_log2_pareto.xlsx', '80% 主分析', '80main')
    all_records += diagnose_one('ori_n165_filtered50_log2_pareto.xlsx', '50% 探索性', '50expl')

    audit = pd.DataFrame(all_records)
    audit_path = TBL_DIR / 'batch_diagnosis_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
