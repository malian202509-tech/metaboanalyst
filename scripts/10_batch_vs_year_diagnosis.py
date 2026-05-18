"""10. 上机批次 vs 采样年份 联合诊断 (再次校正前必跑).

背景: 我们之前用 omx_id 前缀 (A19/B20/C21) 当批次做了 ComBat. 现在拿到了实际上机批次表
(data/01_raw/上机顺序和批次表.xlsx, C18 柱 sheet), 必须重新评估:
  1. 上机批次效应是否真显著? 还是之前测到的"年份效应"已经吸收了?
  2. 批次和年份是否互为冗余?
  3. 应该校正哪个? 还是都校正?

输入:
  data/02_preprocessed/ori_n165_filtered{80,50}_log2_pareto.xlsx   (校正前 Pareto)
  data/02_preprocessed/sample_alignment_n165.csv
  data/01_raw/上机顺序和批次表.xlsx                                  (C18 柱 sheet)

输出:
  results/figures/01_QC_and_Batch_PCA/PCA_batch_vs_year_{80main,50expl}.png
  results/tables/batch_vs_year_diagnosis.csv

判读规则:
  - 若 batch p<0.05 且 R² 高于 year → batch 是主驱动, year 是 batch 的代理 → 改校正 batch
  - 若 batch p>0.1 且 year 仍显著 → batch 弱, year 真实, 保持原 ComBat(year) 即可
  - 若两者都显著且独立 (Cramér's V 低) → 双源, 需联合校正方案 (双 batch 模型 / 合并标签 / 分而治之)
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from sklearn.decomposition import PCA
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
PRE_DIR = ROOT / 'data' / '02_preprocessed'
RAW_DIR = ROOT / 'data' / '01_raw'
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

BATCH_COLORS = {1: '#d62728', 2: '#9467bd'}
YEAR_COLORS = {'A19 (2019)': '#1f77b4', 'B20 (2020)': '#ff7f0e', 'C21 (2021)': '#2ca02c'}
BMI_COLORS = {'正常': '#4C72B0', '超重肥胖': '#DD8452'}


def cramers_v(x, y):
    ct = pd.crosstab(x, y)
    chi2 = stats.chi2_contingency(ct).statistic
    n = ct.values.sum()
    r, k = ct.shape
    return float(np.sqrt(chi2 / (n * (min(r, k) - 1))))


def permanova(X, group, n_perm=999, random_state=42):
    """PERMANOVA (Anderson 2001) on Euclidean distances. 返回 (pseudo-F, p)."""
    rng = np.random.default_rng(random_state)
    g = pd.Series(group).reset_index(drop=True)
    levels = g.unique()
    n = len(g)
    a = len(levels)
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
        if ss_w <= 0 or (n - a) == 0:
            return np.inf
        return (ss_a / (a - 1)) / (ss_w / (n - a))

    F_obs = pseudo_F(g.values)
    cnt = sum(1 for _ in range(n_perm) if pseudo_F(rng.permutation(g.values)) >= F_obs)
    return float(F_obs), float((cnt + 1) / (n_perm + 1))


def anova_pc(pc, group):
    """单因素 ANOVA R² + P."""
    df = pd.DataFrame({'pc': pc, 'g': group})
    grand = df['pc'].mean()
    ss_t = ((df['pc'] - grand) ** 2).sum()
    ss_b = sum(len(sub) * (sub['pc'].mean() - grand) ** 2 for _, sub in df.groupby('g'))
    ss_w = ss_t - ss_b
    L, n = df['g'].nunique(), len(df)
    if ss_w <= 0 or n - L <= 0 or L < 2:
        return np.nan, np.nan
    F = (ss_b / (L - 1)) / (ss_w / (n - L))
    p = 1 - stats.f.cdf(F, L - 1, n - L)
    return float(ss_b / ss_t) if ss_t > 0 else np.nan, float(p)


def load_batch_table():
    """读 C18 柱 sheet, 只保留 Subject 行, 返回 omx_id -> batch / 进样顺序."""
    p = RAW_DIR / '上机顺序和批次表.xlsx'
    df = pd.read_excel(p, sheet_name='C18柱-上机顺序和批次表', header=1, usecols=range(4))
    df.columns = ['order', 'class', 'omx_id', 'batch']
    df['order'] = pd.to_numeric(df['order'], errors='coerce').astype('Int64')
    df['batch'] = pd.to_numeric(df['batch'], errors='coerce').astype('Int64')
    return df[df['class'] == 'Subject'][['omx_id', 'order', 'batch']].reset_index(drop=True)


def hotelling_ellipse(ax, xy, color, alpha=0.95, n_pts=120):
    """Hotelling T² 95% 置信椭圆 (与已有脚本风格一致)."""
    if len(xy) < 3:
        return
    mean = xy.mean(axis=0)
    cov = np.cov(xy.T)
    if not np.all(np.isfinite(cov)) or np.linalg.matrix_rank(cov) < 2:
        return
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    n = len(xy)
    F = stats.f.ppf(alpha, 2, n - 2)
    scale = np.sqrt(2 * (n - 1) / (n - 2) * F)
    theta = np.linspace(0, 2 * np.pi, n_pts)
    circ = np.stack([np.cos(theta), np.sin(theta)], axis=0)
    ell = vecs @ np.diag(scale * np.sqrt(vals)) @ circ
    ax.plot(ell[0] + mean[0], ell[1] + mean[1], '--', color=color, alpha=0.85, lw=1.2)
    ax.fill(ell[0] + mean[0], ell[1] + mean[1], color=color, alpha=0.08)


def diagnose_one(src_name, tag, tag_short, batch_tbl):
    print(f'\n{"="*60}')
    print(f'轨道 {tag}: {src_name}')
    df = pd.read_excel(PRE_DIR / src_name)
    mp = pd.read_csv(PRE_DIR / 'sample_alignment_n165.csv')
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]

    X = df[sample_cols].values.T
    n_samp, n_feat = X.shape
    print(f'  矩阵 (样本×特征): {X.shape}')

    meta = (pd.DataFrame({'omx_id': sample_cols})
              .merge(mp[['omx_id', 'BMI_group']], on='omx_id', how='left')
              .merge(batch_tbl, on='omx_id', how='left'))
    meta['year_prefix'] = meta['omx_id'].str[:3]
    meta['year'] = meta['year_prefix'].map(
        {'A19': 'A19 (2019)', 'B20': 'B20 (2020)', 'C21': 'C21 (2021)'})
    miss = int(meta['batch'].isna().sum())
    if miss > 0:
        print(f'  ⚠ {miss} 样本批次缺失')
    print(f'  batch 分布: {meta["batch"].value_counts().sort_index().to_dict()}')

    # PCA
    n_comp = min(10, n_samp - 1, n_feat)
    pca = PCA(n_components=n_comp, random_state=42)
    pcs = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100
    print(f'  PC1 {var[0]:.1f}% / PC2 {var[1]:.1f}% / PC3 {var[2]:.1f}% '
          f'(前 5 累积 {var[:5].sum():.1f}%)')

    # === 每 PC 的三因素 R² ===
    rows = []
    for k in range(min(5, n_comp)):
        r2_b, p_b = anova_pc(pcs[:, k], meta['batch'])
        r2_y, p_y = anova_pc(pcs[:, k], meta['year'])
        r2_m, p_m = anova_pc(pcs[:, k], meta['BMI_group'])
        rows.append({
            'tier': tag, 'PC': f'PC{k+1}', '解释方差%': round(var[k], 2),
            'batchR²': round(r2_b, 3), 'batch P': f'{p_b:.3g}',
            'yearR²': round(r2_y, 3), 'year P': f'{p_y:.3g}',
            'BMIR²': round(r2_m, 3), 'BMI P': f'{p_m:.3g}',
        })
    quant_df = pd.DataFrame(rows)
    print('  --- 每 PC 三因素 R² + P ---')
    print(quant_df.to_string(index=False))

    # === PERMANOVA: 三因素 ===
    print('  --- PERMANOVA (999 置换) ---')
    F_b, p_b_pm = permanova(X, meta['batch'].astype(str).values, n_perm=999)
    F_y, p_y_pm = permanova(X, meta['year'].values, n_perm=999)
    F_m, p_m_pm = permanova(X, meta['BMI_group'].values, n_perm=999)
    print(f'    batch: F={F_b:.3f}  p={p_b_pm:.3f}')
    print(f'    year : F={F_y:.3f}  p={p_y_pm:.3f}')
    print(f'    BMI  : F={F_m:.3f}  p={p_m_pm:.3f}')

    # === Cramér's V (3x3) ===
    print('  --- Cramér\'s V 矩阵 (共线性) ---')
    v_by = cramers_v(meta['batch'], meta['year'])
    v_bm = cramers_v(meta['batch'], meta['BMI_group'])
    v_ym = cramers_v(meta['year'], meta['BMI_group'])
    print(f'    batch × year : {v_by:.3f}  ({"低" if v_by<0.2 else "中" if v_by<0.5 else "高"})')
    print(f'    batch × BMI  : {v_bm:.3f}  ({"低" if v_bm<0.2 else "中" if v_bm<0.5 else "高"})')
    print(f'    year  × BMI  : {v_ym:.3f}  ({"低" if v_ym<0.2 else "中" if v_ym<0.5 else "高"})')

    # === 可视化 (2x3) ===
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f'批次 vs 年份联合诊断 — {tag}  ({n_samp} × {n_feat}, log2+Pareto, 校正前)',
                 fontsize=14, fontweight='bold')

    def scatter_by(ax, color_map, label_col, title):
        for lv, color in color_map.items():
            m = meta[label_col] == lv
            if m.sum() == 0:
                continue
            xy = pcs[m.values][:, :2]
            ax.scatter(xy[:, 0], xy[:, 1], c=color, label=f'{lv} (n={int(m.sum())})',
                       alpha=0.7, edgecolor='white', linewidth=0.5, s=50)
            hotelling_ellipse(ax, xy, color)
        ax.set_xlabel(f'PC1 ({var[0]:.1f}%)')
        ax.set_ylabel(f'PC2 ({var[1]:.1f}%)')
        ax.set_title(title)
        ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
        ax.legend(loc='best', fontsize=9)

    scatter_by(axes[0, 0], BATCH_COLORS, 'batch',
               f"PC1-PC2 按 batch  (PC1 R²={rows[0]['batchR²']}, P={rows[0]['batch P']})")
    scatter_by(axes[0, 1], YEAR_COLORS, 'year',
               f"PC1-PC2 按 year  (PC1 R²={rows[0]['yearR²']}, P={rows[0]['year P']})")
    scatter_by(axes[0, 2], BMI_COLORS, 'BMI_group',
               f"PC1-PC2 按 BMI  (PC1 R²={rows[0]['BMIR²']}, P={rows[0]['BMI P']})")

    # (1,0) Scree
    ax = axes[1, 0]
    n_show = min(10, len(var))
    ax.bar(range(1, n_show + 1), var[:n_show], color='#5b8ec2', edgecolor='black')
    ax.plot(range(1, n_show + 1), np.cumsum(var[:n_show]), 'o-', color='red', label='累积 %')
    ax.set_xlabel('主成分'); ax.set_ylabel('解释方差 %')
    ax.set_title(f'Scree (前 5 累积 {var[:5].sum():.1f}%)')
    ax.legend(); ax.grid(axis='y', alpha=0.3)

    # (1,1) batch 信号最强那个 PC 的密度
    best_pc_batch = int(np.argmax([float(r['batchR²']) for r in rows]))
    ax = axes[1, 1]
    for lv, color in BATCH_COLORS.items():
        m = meta['batch'] == lv
        if m.sum() > 1:
            ax.hist(pcs[m.values, best_pc_batch], bins=20, alpha=0.45,
                    color=color, label=f'batch {lv}', density=True)
    ax.set_xlabel(f'PC{best_pc_batch+1} ({var[best_pc_batch]:.1f}%)')
    ax.set_ylabel('密度')
    ax.set_title(f"batch 最强 PC{best_pc_batch+1} 密度  "
                 f"(PERMANOVA p={p_b_pm:.3f}, V={v_bm:.2f})")
    ax.legend(fontsize=9)

    # (1,2) year 信号最强那个 PC 的密度 (对照)
    best_pc_year = int(np.argmax([float(r['yearR²']) for r in rows]))
    ax = axes[1, 2]
    for lv, color in YEAR_COLORS.items():
        m = meta['year'] == lv
        if m.sum() > 1:
            ax.hist(pcs[m.values, best_pc_year], bins=20, alpha=0.45,
                    color=color, label=lv, density=True)
    ax.set_xlabel(f'PC{best_pc_year+1} ({var[best_pc_year]:.1f}%)')
    ax.set_ylabel('密度')
    ax.set_title(f"year 最强 PC{best_pc_year+1} 密度  "
                 f"(PERMANOVA p={p_y_pm:.3f}, V={v_by:.2f})")
    ax.legend(fontsize=9)

    plt.tight_layout()
    out_fig = FIG_DIR / f'PCA_batch_vs_year_{tag_short}.png'
    plt.savefig(out_fig, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'  图已保存: {out_fig.relative_to(ROOT)}')

    perm_rows = [
        {'tier': tag, 'PC': 'PERMANOVA-batch', '解释方差%': '',
         'batchR²': round(F_b, 3), 'batch P': f'{p_b_pm:.3f}',
         'yearR²': '', 'year P': '', 'BMIR²': '', 'BMI P': ''},
        {'tier': tag, 'PC': 'PERMANOVA-year', '解释方差%': '',
         'batchR²': '', 'batch P': '',
         'yearR²': round(F_y, 3), 'year P': f'{p_y_pm:.3f}',
         'BMIR²': '', 'BMI P': ''},
        {'tier': tag, 'PC': 'PERMANOVA-BMI', '解释方差%': '',
         'batchR²': '', 'batch P': '', 'yearR²': '', 'year P': '',
         'BMIR²': round(F_m, 3), 'BMI P': f'{p_m_pm:.3f}'},
        {'tier': tag, 'PC': "Cramer's V (batch×year)", '解释方差%': '',
         'batchR²': round(v_by, 3), 'batch P': '',
         'yearR²': round(v_by, 3), 'year P': '', 'BMIR²': '', 'BMI P': ''},
        {'tier': tag, 'PC': "Cramer's V (batch×BMI)", '解释方差%': '',
         'batchR²': round(v_bm, 3), 'batch P': '', 'yearR²': '', 'year P': '',
         'BMIR²': round(v_bm, 3), 'BMI P': ''},
        {'tier': tag, 'PC': "Cramer's V (year×BMI)", '解释方差%': '',
         'batchR²': '', 'batch P': '',
         'yearR²': round(v_ym, 3), 'year P': '',
         'BMIR²': round(v_ym, 3), 'BMI P': ''},
    ]
    return rows + perm_rows


def main():
    print('=== 10. 批次 vs 年份 联合诊断 (双轨) ===')
    batch_tbl = load_batch_table()
    print(f'C18 柱 Subject 批次表: {len(batch_tbl)} 行')
    print(f'批次分布 (全 198 Subject): {batch_tbl["batch"].value_counts().sort_index().to_dict()}')

    all_rows = []
    all_rows += diagnose_one('ori_n165_filtered80_log2_pareto.xlsx', '80% 主分析', '80main', batch_tbl)
    all_rows += diagnose_one('ori_n165_filtered50_log2_pareto.xlsx', '50% 探索性', '50expl', batch_tbl)

    audit = pd.DataFrame(all_rows)
    audit_path = TBL_DIR / 'batch_vs_year_diagnosis.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
