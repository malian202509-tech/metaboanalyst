"""07. OPLS-DA SIMCA-P 风格置换检验.

每条数据轨生成一张 2x2 大图:
  (0,0) OPLS-DA scores: t_pred vs t_ortho 按 BMI + Hotelling 95% 椭圆
  (0,1) 置换检验主图: |corr(perm_y, y)| 在 x 轴, R²Y 与 Q² 在 y 轴
        画散点 + 回归线 + 在 corr=1 处标出实测值 + 在 corr=0 处标出截距
  (1,0) R²Y 置换分布直方图 + 实测值红线 (PLS-DA 风格)
  (1,1) Q²  置换分布直方图 + 实测值红线

SIMCA-P 报告标准:
  - R²Y intercept (corr=0) < 0.30~0.40  -> 模型不过拟合
  - Q²  intercept (corr=0) < 0.05        -> 预测能力可靠
  - P_R²Y, P_Q² < 0.05                  -> 显著优于随机

输入:
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat_pareto.xlsx
  data/02_preprocessed/sample_alignment_n165.csv

输出:
  results/figures/01_QC_and_Batch_PCA/OPLS_permutation_{80main,50expl}.png
  results/tables/opls_permutation_audit.csv
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Ellipse
from pathlib import Path
from scipy.stats import f as f_dist

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
from utils.opls import opls_da, opls_da_q2, permutation_test_opls

POST_DIR = ROOT / 'data' / '03_batch_corrected'
PRE_DIR = ROOT / 'data' / '02_preprocessed'
FIG_DIR = ROOT / 'results' / 'figures' / '01_QC_and_Batch_PCA'
TBL_DIR = ROOT / 'results' / 'tables'

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'

mpl.rcParams['font.family'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
mpl.rcParams['axes.unicode_minus'] = False
BMI_COLORS = {'正常': '#4C72B0', '超重肥胖': '#DD8452'}

N_PERM = 200
N_FOLDS = 7
N_ORTHO = 1


def hotelling_ellipse(x, y, conf=0.95):
    pts = np.column_stack([x, y])
    n = len(pts)
    if n < 5:
        return None
    mu = pts.mean(axis=0)
    cov = np.cov(pts.T, ddof=1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]; eigvecs = eigvecs[:, order]
    p = 2
    t2_crit = (p * (n - 1) / (n - p)) * f_dist.ppf(conf, p, n - p)
    semi = np.sqrt(np.maximum(eigvals, 0) * t2_crit)
    angle = float(np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0])))
    return mu, 2 * semi[0], 2 * semi[1], angle


def add_ellipse(ax, x, y, color, conf=0.95):
    res = hotelling_ellipse(x, y, conf)
    if res is None:
        return
    mu, w, h, ang = res
    ax.add_patch(Ellipse(xy=mu, width=w, height=h, angle=ang,
                         facecolor=color, alpha=0.12,
                         edgecolor=color, linestyle='--', linewidth=1.2))


def load_matrix(path):
    df = pd.read_excel(path)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    return df[sample_cols].values.T, sample_cols


def validate_one(tier, tag_short, post_path):
    print(f'\n{"="*60}')
    print(f'轨道: {tier}')
    X, sample_ids = load_matrix(post_path)
    mp = pd.read_csv(PRE_DIR / 'sample_alignment_n165.csv')
    meta = pd.DataFrame({'omx_id': sample_ids}).merge(
        mp[['omx_id', 'BMI_group']], on='omx_id', how='left'
    )
    y = (meta['BMI_group'] == '超重肥胖').astype(int).values
    print(f'  X shape: {X.shape}, 正常={int((y==0).sum())}, 超重肥胖={int((y==1).sum())}')

    # OPLS-DA fit (用于 scores 图)
    print(f'  --- OPLS-DA fit ({N_ORTHO} orthogonal) ---')
    model = opls_da(X, y, n_ortho=N_ORTHO)
    print(f'    R²X_pred = {model["r2x_pred"]:.3f}   R²X_ortho = {model["r2x_ortho"]:.3f}   R²Y = {model["r2y"]:.3f}')

    # Q² CV
    print(f'  --- Q² via {N_FOLDS}-fold CV ---')
    q2 = opls_da_q2(X, y, n_ortho=N_ORTHO, n_folds=N_FOLDS, random_state=42)
    print(f'    Q² = {q2:.3f}')

    # 置换检验
    print(f'  --- 置换检验 ({N_PERM} 次, 每次 {N_FOLDS}-fold CV) ---')
    perm = permutation_test_opls(X, y, n_ortho=N_ORTHO, n_perm=N_PERM,
                                  n_folds=N_FOLDS, random_state=42)
    print(f'    R²Y_obs = {perm["r2y_obs"]:.3f}   p_R²Y = {perm["p_r2y"]:.3f}   R²Y intercept(corr=0) = {perm["r2y_intercept"]:.3f}')
    print(f'    Q²_obs  = {perm["q2_obs"]:.3f}    p_Q²  = {perm["p_q2"]:.3f}   Q²  intercept(corr=0) = {perm["q2_intercept"]:.3f}')

    # 判定
    r2y_pass_intercept = perm['r2y_intercept'] < 0.4
    q2_pass_intercept = perm['q2_intercept'] < 0.05
    r2y_pass_p = perm['p_r2y'] < 0.05
    q2_pass_p = perm['p_q2'] < 0.05
    overall = '通过' if (r2y_pass_intercept and q2_pass_intercept and r2y_pass_p and q2_pass_p) else '过拟合 / 不可靠'
    print(f'  判定: R²Y截距<0.4={r2y_pass_intercept}, Q²截距<0.05={q2_pass_intercept}, '
          f'p_R²Y<0.05={r2y_pass_p}, p_Q²<0.05={q2_pass_p}  →  {overall}')

    # === 2x2 大图 ===
    fig, ax = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(f'OPLS-DA 置换检验 SIMCA-P 风格 — {tier}  ({N_PERM} 次置换, {N_FOLDS}-fold CV)',
                 fontsize=14, fontweight='bold')

    # (0,0) OPLS-DA scores 图
    a = ax[0, 0]
    for lv, color in BMI_COLORS.items():
        m = meta['BMI_group'].values == lv
        if m.sum() == 0: continue
        add_ellipse(a, model['t_pred'][m], model['t_ortho'][m, 0], color)
        a.scatter(model['t_pred'][m], model['t_ortho'][m, 0], c=color, alpha=0.7,
                  edgecolor='white', linewidth=0.5, s=45,
                  label=f'{lv} (n={int(m.sum())})')
    a.axhline(0, color='gray', lw=0.5)
    a.axvline(0, color='gray', lw=0.5)
    a.set_xlabel(f't_predictive (R²X={model["r2x_pred"]:.3f})')
    a.set_ylabel(f't_orthogonal (R²X={model["r2x_ortho"]:.3f})')
    a.set_title(f'OPLS-DA scores  R²Y={model["r2y"]:.3f}, Q²={q2:.3f}')
    a.legend(fontsize=9)

    # (0,1) SIMCA-P 主图: corr 上的 R²Y / Q²
    a = ax[0, 1]
    x_all_r = np.concatenate([perm['correlations'], [1.0]])
    y_all_r = np.concatenate([perm['r2y_perm'], [perm['r2y_obs']]])
    y_all_q = np.concatenate([perm['q2_perm'], [perm['q2_obs']]])
    a.scatter(perm['correlations'], perm['r2y_perm'], c='#1f77b4', alpha=0.5, s=20, label='R²Y 置换')
    a.scatter(perm['correlations'], perm['q2_perm'], c='#2ca02c', alpha=0.5, s=20, label='Q² 置换')
    a.scatter([1.0], [perm['r2y_obs']], c='#1f77b4', s=120, marker='*', edgecolor='black',
              linewidth=1, label=f'R²Y 实测 = {perm["r2y_obs"]:.3f}')
    a.scatter([1.0], [perm['q2_obs']], c='#2ca02c', s=120, marker='*', edgecolor='black',
              linewidth=1, label=f'Q² 实测 = {perm["q2_obs"]:.3f}')
    # 回归线
    xx = np.array([0, 1.05])
    a.plot(xx, np.polyval(np.polyfit(x_all_r, y_all_r, 1), xx), '--', color='#1f77b4', lw=1.2)
    a.plot(xx, np.polyval(np.polyfit(x_all_r, y_all_q, 1), xx), '--', color='#2ca02c', lw=1.2)
    # 截距标注
    a.axhline(0, color='gray', lw=0.5)
    a.axvline(0, color='gray', lw=0.5)
    a.annotate(f'R²Y 截距 = {perm["r2y_intercept"]:.3f}  ({"OK" if perm["r2y_intercept"]<0.4 else "FAIL (<0.4)"})',
               xy=(0.02, perm['r2y_intercept']), fontsize=9, color='#1f77b4',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    a.annotate(f'Q² 截距 = {perm["q2_intercept"]:.3f}  ({"OK" if perm["q2_intercept"]<0.05 else "FAIL (<0.05)"})',
               xy=(0.02, perm['q2_intercept'] - 0.05), fontsize=9, color='#2ca02c',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    a.set_xlabel('|corr(perm y, original y)|')
    a.set_ylabel('R²Y / Q²')
    a.set_title('SIMCA-P 置换主图  (回归线截距是判定关键)')
    a.legend(fontsize=8, loc='best')
    a.grid(alpha=0.3)

    # (1,0) R²Y 分布
    a = ax[1, 0]
    a.hist(perm['r2y_perm'], bins=30, color='#5b8ec2', edgecolor='black', alpha=0.7)
    a.axvline(perm['r2y_obs'], color='red', lw=2, label=f'实测 = {perm["r2y_obs"]:.3f}')
    title_c = 'green' if perm['p_r2y'] < 0.05 else 'red'
    a.set_title(f'R²Y 置换分布  p = {perm["p_r2y"]:.3f}  '
                f'({"通过" if perm["p_r2y"] < 0.05 else "未通过"})',
                color=title_c, fontsize=10)
    a.set_xlabel('R²Y'); a.set_ylabel('频数')
    a.legend(); a.grid(axis='y', alpha=0.3)

    # (1,1) Q² 分布
    a = ax[1, 1]
    a.hist(perm['q2_perm'], bins=30, color='#7cc77c', edgecolor='black', alpha=0.7)
    a.axvline(perm['q2_obs'], color='red', lw=2, label=f'实测 = {perm["q2_obs"]:.3f}')
    title_c = 'green' if perm['p_q2'] < 0.05 else 'red'
    a.set_title(f'Q² 置换分布  p = {perm["p_q2"]:.3f}  '
                f'({"通过" if perm["p_q2"] < 0.05 else "未通过"})',
                color=title_c, fontsize=10)
    a.set_xlabel('Q²'); a.set_ylabel('频数')
    a.legend(); a.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    out_fig = FIG_DIR / f'OPLS_permutation_{tag_short}.png'
    plt.savefig(out_fig, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'  图已保存: {out_fig.relative_to(ROOT)}')

    return {
        'tier': tier,
        'R2X_pred': round(model['r2x_pred'], 3),
        'R2X_ortho': round(model['r2x_ortho'], 3),
        'R2Y_obs': round(perm['r2y_obs'], 3),
        'Q2_obs': round(perm['q2_obs'], 3),
        'p_R2Y': round(perm['p_r2y'], 4),
        'p_Q2': round(perm['p_q2'], 4),
        'R2Y_intercept(corr=0)': round(perm['r2y_intercept'], 3),
        'Q2_intercept(corr=0)': round(perm['q2_intercept'], 3),
        'R2Y截距<0.4': perm['r2y_intercept'] < 0.4,
        'Q2截距<0.05': perm['q2_intercept'] < 0.05,
        '判定': overall,
    }


def main():
    print('=== 07. OPLS-DA SIMCA-P 风格置换检验 (双轨并行) ===')
    print(f'    参数: n_orthogonal={N_ORTHO}, n_permutations={N_PERM}, n_folds={N_FOLDS}')
    records = []
    records.append(validate_one(
        '80% 主分析', '80main',
        POST_DIR / 'ori_n165_filtered80_log2_combat_pareto.xlsx',
    ))
    records.append(validate_one(
        '50% 探索性', '50expl',
        POST_DIR / 'ori_n165_filtered50_log2_combat_pareto.xlsx',
    ))
    audit = pd.DataFrame(records)
    audit_path = TBL_DIR / 'opls_permutation_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')
    print()
    cols_print = ['tier', 'R2Y_obs', 'Q2_obs', 'p_R2Y', 'p_Q2',
                  'R2Y_intercept(corr=0)', 'Q2_intercept(corr=0)', '判定']
    print(audit[cols_print].to_string(index=False))


if __name__ == '__main__':
    main()
