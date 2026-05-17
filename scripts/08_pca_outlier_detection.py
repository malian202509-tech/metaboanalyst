"""08. PCA 离群值诊断 (Hotelling T² + DModX 双重判定).

依据:
  - Wiklund et al. Anal Chem 2008: T² + DModX 是 SIMCA-P 标准做法
  - Brereton, Chemometrics 2003: T² 抓 strong outlier; DModX 抓 moderate outlier
  - 保守阈值: T² > 99% AND DModX > 1.5× 均值 (双命中) -> 严格离群候选
            T² > 95% AND DModX > 1.5× 均值 (双命中) -> 中等离群候选

策略:
  - 不自动剔除, 只输出候选清单
  - 候选用户人工核对原始病历/上机报告:
      * 有技术问题 -> 剔除
      * 真实生物学极端 -> 保留, 在 Limitations 说明
  - 后续敏感性分析: 主分析 with vs without outlier

输入:
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat_pareto.xlsx
  data/02_preprocessed/sample_alignment_n165.csv

输出:
  results/figures/01_QC_and_Batch_PCA/Outlier_diagnosis_{80main,50expl}.png
  results/tables/outlier_candidates.csv
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
POST_DIR = ROOT / 'data' / '03_batch_corrected'
PRE_DIR = ROOT / 'data' / '02_preprocessed'
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
BMI_COLORS = {'正常': '#4C72B0', '超重肥胖': '#DD8452'}

K_COMPONENTS = 5            # 用前 5 个 PC 计算 T² / DModX (约解释 50% 方差)
DMODX_RATIO_THRESH = 1.5    # DModX 临界: > 1.5 × 均值
T2_ALPHA_STRICT = 0.01      # T² 严格临界 99%
T2_ALPHA_LOOSE = 0.05       # T² 宽松临界 95%


def hotelling_t2(scores, eigvals):
    """每样本 T² = sum_k (t_ki² / λ_k)."""
    return np.sum((scores ** 2) / eigvals.reshape(1, -1), axis=1)


def t2_critical(K, n, alpha):
    return K * (n - 1) / (n - K) * f_dist.ppf(1 - alpha, K, n - K)


def dmodx_per_sample(X, scores, loadings, K):
    """DModX_i = sqrt( sum_j (x_ij - x_hat_ij)² / (n_features - K) )."""
    X_hat = scores @ loadings.T
    residuals = X - X_hat
    n_feat = X.shape[1]
    if n_feat <= K:
        return np.full(X.shape[0], np.nan)
    return np.sqrt(np.sum(residuals ** 2, axis=1) / (n_feat - K))


def hotelling_ellipse_2d(scores_2d, conf):
    """2D Hotelling 置信椭圆 (用于 PCA 散点)."""
    n = len(scores_2d)
    if n < 5: return None
    mu = scores_2d.mean(axis=0)
    cov = np.cov(scores_2d.T, ddof=1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]; eigvecs = eigvecs[:, order]
    p = 2
    t2c = (p * (n - 1) / (n - p)) * f_dist.ppf(conf, p, n - p)
    semi = np.sqrt(np.maximum(eigvals, 0) * t2c)
    angle = float(np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0])))
    return mu, 2 * semi[0], 2 * semi[1], angle


def load_data(post_path):
    df = pd.read_excel(post_path)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    return df[sample_cols].values.T, sample_cols  # X (n_samp, n_feat)


def diagnose_one(tier, tag_short, post_path):
    print(f'\n{"="*60}')
    print(f'轨道: {tier}')
    X, sample_ids = load_data(post_path)
    mp = pd.read_csv(PRE_DIR / 'sample_alignment_n165.csv')
    n_samp, n_feat = X.shape
    print(f'  样本数: {n_samp}, 特征数: {n_feat}, 使用 K={K_COMPONENTS} 个 PC')

    meta = pd.DataFrame({'omx_id': sample_ids}).merge(
        mp[['omx_id', '标本号', 'BMI_group']], on='omx_id', how='left'
    )
    meta['year'] = meta['omx_id'].str[:3]

    # PCA
    pca = PCA(n_components=K_COMPONENTS, random_state=42)
    scores = pca.fit_transform(X)
    eigvals = pca.explained_variance_                # σ²(t_k)
    loadings = pca.components_.T                     # (n_feat, K)
    var_pct = pca.explained_variance_ratio_ * 100
    print(f'  PC 解释方差 (前{K_COMPONENTS}): {var_pct.round(1)}  (累积 {var_pct.sum():.1f}%)')

    # T² + 临界
    t2 = hotelling_t2(scores, eigvals)
    t2c_loose = t2_critical(K_COMPONENTS, n_samp, T2_ALPHA_LOOSE)
    t2c_strict = t2_critical(K_COMPONENTS, n_samp, T2_ALPHA_STRICT)
    print(f'  T² 临界 95%={t2c_loose:.2f}  99%={t2c_strict:.2f}')
    print(f'  T² > 95%: {(t2 > t2c_loose).sum()} 样本  | > 99%: {(t2 > t2c_strict).sum()} 样本')

    # DModX + 临界
    # X 已被 Pareto 缩放, 列均值约为 0; 但 PCA 内部还会再中心化, 这里直接传 X 给 PCA
    # DModX 需要 X 减去 mean 后与 X_hat 比较, 但 sklearn PCA 的 transform 是基于中心化后的
    X_centered = X - X.mean(axis=0)
    dmodx = dmodx_per_sample(X_centered, scores, loadings, K_COMPONENTS)
    dmodx_mean = float(np.mean(dmodx))
    dmodx_crit = DMODX_RATIO_THRESH * dmodx_mean
    print(f'  DModX 均值={dmodx_mean:.3f}  临界 (1.5x)={dmodx_crit:.3f}')
    print(f'  DModX > 临界: {(dmodx > dmodx_crit).sum()} 样本')

    # 双重命中
    hit_strict = (t2 > t2c_strict) & (dmodx > dmodx_crit)
    hit_loose = (t2 > t2c_loose) & (dmodx > dmodx_crit)
    print(f'  双重命中 (严格 T²>99% & DModX>1.5x): {hit_strict.sum()} 个')
    print(f'  双重命中 (宽松 T²>95% & DModX>1.5x): {hit_loose.sum()} 个')

    # === 4 子图诊断 ===
    fig, ax = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(f'PCA 离群值诊断 — {tier}  ({n_samp} 样本 × {n_feat} 代谢物)',
                 fontsize=14, fontweight='bold')

    # (0,0) PCA 散点 + 95%/99% Hotelling 椭圆 + 候选标签
    a = ax[0, 0]
    for lv, color in BMI_COLORS.items():
        m = meta['BMI_group'].values == lv
        if m.sum() == 0: continue
        # 95% (实线) + 99% (虚线) 椭圆
        res95 = hotelling_ellipse_2d(scores[m, :2], 0.95)
        res99 = hotelling_ellipse_2d(scores[m, :2], 0.99)
        if res95:
            mu, w, h, ang = res95
            a.add_patch(Ellipse(xy=mu, width=w, height=h, angle=ang,
                                facecolor=color, alpha=0.10,
                                edgecolor=color, linestyle='-', linewidth=1.2,
                                label=None))
        if res99:
            mu, w, h, ang = res99
            a.add_patch(Ellipse(xy=mu, width=w, height=h, angle=ang,
                                facecolor='none',
                                edgecolor=color, linestyle=':', linewidth=1))
        a.scatter(scores[m, 0], scores[m, 1], c=color, alpha=0.7,
                  edgecolor='white', linewidth=0.5, s=45,
                  label=f'{lv} (n={int(m.sum())})')
    # 标候选离群
    for idx in np.where(hit_loose)[0]:
        sn = meta.iloc[idx]['标本号']
        a.annotate(str(sn), (scores[idx, 0], scores[idx, 1]),
                   fontsize=7, color='red', alpha=0.85,
                   xytext=(4, 4), textcoords='offset points')
    a.axhline(0, color='gray', lw=0.5); a.axvline(0, color='gray', lw=0.5)
    a.set_xlabel(f'PC1 ({var_pct[0]:.1f}%)'); a.set_ylabel(f'PC2 ({var_pct[1]:.1f}%)')
    a.set_title(f'PCA scores + Hotelling 95% (实线) / 99% (虚线) 椭圆')
    a.legend(fontsize=8)

    # (0,1) T² vs DModX 散点
    a = ax[0, 1]
    a.scatter(t2, dmodx, c='#4C72B0', alpha=0.6, s=40, edgecolor='white', linewidth=0.5)
    a.axvline(t2c_loose, color='orange', ls='--', lw=1, label=f'T² 95% = {t2c_loose:.1f}')
    a.axvline(t2c_strict, color='red', ls='--', lw=1, label=f'T² 99% = {t2c_strict:.1f}')
    a.axhline(dmodx_crit, color='green', ls='--', lw=1, label=f'DModX 1.5× = {dmodx_crit:.2f}')
    # 候选标签
    for idx in np.where(hit_loose)[0]:
        sn = meta.iloc[idx]['标本号']
        marker_color = 'red' if hit_strict[idx] else 'orange'
        a.scatter(t2[idx], dmodx[idx], c=marker_color, s=80, marker='o',
                  edgecolor='black', linewidth=1, zorder=3)
        a.annotate(str(sn), (t2[idx], dmodx[idx]),
                   fontsize=7, color=marker_color, fontweight='bold',
                   xytext=(5, 5), textcoords='offset points')
    a.set_xlabel('Hotelling T²'); a.set_ylabel('DModX')
    a.set_title(f"双重诊断  红=严格({hit_strict.sum()}例)  橙=宽松({hit_loose.sum() - hit_strict.sum()}例)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (1,0) T² 排序 (top 15)
    a = ax[1, 0]
    order = np.argsort(t2)[::-1][:15]
    labels = [str(meta.iloc[i]['标本号']) for i in order]
    colors = ['red' if t2[i] > t2c_strict else ('orange' if t2[i] > t2c_loose else '#4C72B0') for i in order]
    a.barh(range(len(order)), t2[order], color=colors, edgecolor='black', linewidth=0.5)
    a.set_yticks(range(len(order))); a.set_yticklabels(labels, fontsize=8)
    a.axvline(t2c_loose, color='orange', ls='--', lw=1)
    a.axvline(t2c_strict, color='red', ls='--', lw=1)
    a.invert_yaxis()
    a.set_xlabel('Hotelling T²'); a.set_title('T² Top 15')
    a.grid(axis='x', alpha=0.3)

    # (1,1) DModX 排序 (top 15)
    a = ax[1, 1]
    order = np.argsort(dmodx)[::-1][:15]
    labels = [str(meta.iloc[i]['标本号']) for i in order]
    colors = ['red' if (dmodx[i] > dmodx_crit and t2[i] > t2c_strict)
              else ('orange' if (dmodx[i] > dmodx_crit and t2[i] > t2c_loose) else '#7cc77c')
              for i in order]
    a.barh(range(len(order)), dmodx[order], color=colors, edgecolor='black', linewidth=0.5)
    a.set_yticks(range(len(order))); a.set_yticklabels(labels, fontsize=8)
    a.axvline(dmodx_crit, color='green', ls='--', lw=1)
    a.invert_yaxis()
    a.set_xlabel('DModX'); a.set_title('DModX Top 15')
    a.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    out_fig = FIG_DIR / f'Outlier_diagnosis_{tag_short}.png'
    plt.savefig(out_fig, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'  图已保存: {out_fig.relative_to(ROOT)}')

    # 候选清单 (所有"宽松命中"以上)
    cand = meta.copy()
    cand['tier'] = tier
    cand['T2'] = t2.round(3)
    cand['T2_临界95%'] = round(t2c_loose, 3)
    cand['T2_临界99%'] = round(t2c_strict, 3)
    cand['DModX'] = dmodx.round(3)
    cand['DModX_临界'] = round(dmodx_crit, 3)
    cand['hit_strict (T²>99% & DModX>1.5x)'] = hit_strict
    cand['hit_loose (T²>95% & DModX>1.5x)'] = hit_loose
    # 只输出至少宽松命中的样本
    cand_out = cand[hit_loose].copy()
    cand_out['优先级'] = cand_out['hit_strict (T²>99% & DModX>1.5x)'].map({True: '🔴 严格', False: '🟠 宽松'})
    return cand_out


def main():
    print('=== 08. PCA 离群值诊断 (Hotelling T² + DModX 双重判定) ===')
    print(f'    参数: K_components={K_COMPONENTS}, T² alpha 严格={T2_ALPHA_STRICT}/宽松={T2_ALPHA_LOOSE}, DModX ratio={DMODX_RATIO_THRESH}')

    all_cands = []
    all_cands.append(diagnose_one(
        '80% 主分析', '80main',
        POST_DIR / 'ori_n165_filtered80_log2_combat_pareto.xlsx',
    ))
    all_cands.append(diagnose_one(
        '50% 探索性', '50expl',
        POST_DIR / 'ori_n165_filtered50_log2_combat_pareto.xlsx',
    ))

    if any(len(c) > 0 for c in all_cands):
        combined = pd.concat(all_cands, ignore_index=True)
        cols = ['tier', 'omx_id', '标本号', 'BMI_group', 'year',
                'T2', 'T2_临界95%', 'T2_临界99%',
                'DModX', 'DModX_临界',
                'hit_strict (T²>99% & DModX>1.5x)',
                'hit_loose (T²>95% & DModX>1.5x)',
                '优先级']
        combined = combined[cols].sort_values(['tier', 'T2'], ascending=[True, False])
        out_csv = TBL_DIR / 'outlier_candidates.csv'
        combined.to_csv(out_csv, index=False, encoding='utf-8-sig')
        print(f'\n候选清单已写入: {out_csv.relative_to(ROOT)}')
        print()
        print('--- 候选汇总 ---')
        print(combined.to_string(index=False))
    else:
        print('\n两条轨道均无离群候选 (双重命中 0).')


if __name__ == '__main__':
    main()
