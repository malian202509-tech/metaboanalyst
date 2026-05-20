"""14. 候选差异代谢物 boxplot + partial residual plot (SOP §6.2).

为 4 个候选 (50 探索轨, 含 80 主轨的 2 个) 各画 2 个子图:
  - 左: 原始 log2 浓度 boxplot + 散点 jitter (每个样本一点)
  - 右: partial residual plot (剔除 age/GA/year 协变量后, BMI 的纯效应)
        partial_resid = y - (Xβ - β_BMI · BMI) = β_BMI · BMI + residual

候选 (50 探索轨, 按 p_limma 排序):
  1. 5,6-DiHETrE          (AA-CYP/sEH, 50 only)
  2. 12-HETrE             (AA-LOX, 在 80 主轨)
  3. 20-HDoHE             (ω-3 PUFA oxylipins, 在 80 主轨)
  4. 6,15-Diketo-PGF1α    (AA-COX, 50 only)

输出:
  results/figures/04_Metabolite_boxplots/candidate_boxplots_grid.png (4×2 大图, 论文 Figure 候选)
  results/figures/04_Metabolite_boxplots/candidate_{short}.png (每个候选单独一张)
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
COMBAT_DIR = ROOT / 'data' / '03_batch_corrected'
PREP_DIR = ROOT / 'data' / '02_preprocessed'
RAW_DIR = ROOT / 'data' / '01_raw'
TABLES = ROOT / 'results' / 'tables'
FIG_DIR = ROOT / 'results' / 'figures' / '04_Metabolite_boxplots'
FIG_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'
BMI_COLORS = {'正常': '#67A9CF', '超重肥胖': '#EF8A62'}

# 候选清单 (50 探索轨 4 个), 按 p_limma 升序
CANDIDATES = [
    ('5,6-DiHydroxy-8Z,11Z,14Z-eicosatrienoic acid',          '5,6-DiHETrE',         'AA-CYP/sEH'),
    ('12-Hydroxy-8Z,10E,14Z-eicosatrienoic acid',             '12-HETrE',            'AA-LOX'),
    ('20-Hydroxy-4Z,7Z,10Z,13Z,16Z,18E-docosahexaenoic acid', '20-HDoHE',            'ω-3 PUFA oxylipins'),
    ('6,15-Diketo-13,14-dihydro-prostaglandin F1α',           '6,15-Diketo-PGF1α',  'AA-COX'),
]


def load_data():
    df = pd.read_excel(COMBAT_DIR / 'ori_n165_filtered50_log2_combat.xlsx')
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]

    align = pd.read_csv(PREP_DIR / 'sample_alignment_n165.csv')
    clin = pd.read_excel(RAW_DIR / '分娩特征分析_2019-2021_主分析队列n165.xlsx')

    cov = (clin[['标本号', '孕妇年龄', 'GA_decimal', 'BMI_group']]
           .rename(columns={'孕妇年龄': 'age'})
           .merge(align[['标本号', 'omx_id']], on='标本号', how='left'))
    cov['BMI_overweight'] = (cov['BMI_group'] == '超重肥胖').astype(float)
    yp = cov['omx_id'].str[:3]
    cov['year_A19'] = (yp == 'A19').astype(float)
    cov['year_C21'] = (yp == 'C21').astype(float)

    cov_idx = cov.set_index('omx_id').loc[sample_cols]

    anc = pd.read_csv(TABLES / 'ancova_main_50.csv', encoding='utf-8-sig').set_index('Metabolite Name')

    # 真·80 轨候选集合 (区别于 also_in_80_main 表达的"在 80 轨特征过滤集中")
    cand80 = pd.read_csv(TABLES / 'diff_candidates_80.csv', encoding='utf-8-sig')
    main_cand_set = set(cand80['Metabolite Name'])

    return df, sample_cols, cov_idx, anc, main_cand_set


def compute_partial_residual(y, X, contrast_col):
    """剔除非 BMI 协变量效应后, 返回 partial residual = β_BMI · BMI + ε"""
    model = sm.OLS(y, X).fit()
    beta_bmi = model.params[contrast_col]
    bmi_vals = X[contrast_col].values
    # full residual + β_BMI 项 = partial residual for BMI
    partial = model.resid + beta_bmi * bmi_vals
    return partial


def draw_box(ax, y, groups, palette, ylabel, title_lines):
    """单子图: boxplot + 散点 jitter."""
    grp_order = ['正常', '超重肥胖']
    data = [y[groups == g] for g in grp_order]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.55, patch_artist=True,
                    boxprops=dict(linewidth=1.0),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=1.0),
                    capprops=dict(linewidth=1.0),
                    flierprops=dict(marker='', markersize=0))
    for b, g in zip(bp['boxes'], grp_order):
        b.set_facecolor(palette[g])
        b.set_alpha(0.5)
        b.set_edgecolor(palette[g])
    # jitter 散点
    rng = np.random.default_rng(42)
    for x_pos, g in enumerate(grp_order):
        yvals = y[groups == g]
        x_jit = rng.normal(x_pos, 0.07, size=len(yvals))
        ax.scatter(x_jit, yvals, s=18, c=palette[g], edgecolors='black',
                   linewidths=0.4, alpha=0.65, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'{g}\n(n={int((groups == g).sum())})' for g in grp_order], fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title('\n'.join(title_lines), fontsize=10, pad=8)
    ax.grid(axis='y', alpha=0.3, linestyle='--')


def main():
    print('=== 14. 候选 boxplot + partial residual ===\n')
    df, sample_cols, cov, anc, main_cand_set = load_data()
    print(f'  样本: {len(sample_cols)};  协变量列: age, GA_decimal, year_A19, year_C21')

    groups = cov['BMI_group'].values
    X_design = sm.add_constant(cov[['BMI_overweight', 'age', 'GA_decimal',
                                     'year_A19', 'year_C21']].astype(float))

    # 总 figure: 4 行 × 2 列
    fig, axes = plt.subplots(4, 2, figsize=(11, 16))
    palette = BMI_COLORS

    for row_idx, (name, short, family) in enumerate(CANDIDATES):
        y = df.loc[df['Metabolite Name'] == name, sample_cols].iloc[0].values.astype(float)
        a = anc.loc[name]
        fc = a['log2FC']
        ci_lo, ci_hi = a['log2FC_CI_lo'], a['log2FC_CI_hi']
        p_ols = a['p_ols_hc3']
        p_lim = a['p_limma']
        p_wil = a['p_wilcoxon']
        q_lim = a['q_limma_bh']
        is_in_main_candidates = name in main_cand_set

        # 子图 1: 原始 log2 boxplot
        title_left = [
            f'{short}  ({family})',
            f'Raw log2 expression  (n={len(y)})',
        ]
        draw_box(axes[row_idx, 0], y, groups, palette,
                 ylabel='log2 (concentration, ComBat-adjusted)',
                 title_lines=title_left)

        # 子图 2: partial residual plot
        partial = compute_partial_residual(y, X_design, 'BMI_overweight')
        # 把组均值加回来让 y 轴更直观
        partial_shift = partial + np.mean(y)
        title_right = [
            f'log2FC = {fc:+.3f}  [{ci_lo:+.3f}, {ci_hi:+.3f}]',
            f'p_OLS={p_ols:.4f}  p_limma={p_lim:.4f}  p_Wilcox={p_wil:.4f}',
        ]
        draw_box(axes[row_idx, 1], partial_shift, groups, palette,
                 ylabel='Partial residual\n(adjusted for age, GA, year)',
                 title_lines=title_right)

        # 右上角标签: 是真·80 主轨候选 (核心一致) 还是 50 only (低丰度抢救)
        if is_in_main_candidates:
            tag_text = '★★ Core candidate\n(80% main + 50% explor)'
            tag_color = '#D62728'
        else:
            tag_text = '★ 50% exploratory only\n(low-abundance rescue)'
            tag_color = '#FF7F0E'
        axes[row_idx, 1].text(0.98, 0.97, tag_text, transform=axes[row_idx, 1].transAxes,
                              ha='right', va='top', fontsize=8.5,
                              color=tag_color, fontweight='bold',
                              bbox=dict(facecolor='white', alpha=0.85, edgecolor=tag_color))

    fig.suptitle('Candidate differentially-abundant oxylipins (50% exploratory track, n=4)\n'
                 'Raw log2 vs partial residual after ANCOVA adjustment',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.04, hspace=0.55, wspace=0.30)
    grid_path = FIG_DIR / 'candidate_boxplots_grid.png'
    plt.savefig(grid_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {grid_path.relative_to(ROOT)}')

    # 每个候选单独一张 (论文 Figure 拆分用)
    for name, short, family in CANDIDATES:
        y = df.loc[df['Metabolite Name'] == name, sample_cols].iloc[0].values.astype(float)
        a = anc.loc[name]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        title_left = [f'{short}  ({family})', f'Raw log2  (n={len(y)})']
        draw_box(axes[0], y, groups, palette,
                 ylabel='log2 (concentration)', title_lines=title_left)
        partial = compute_partial_residual(y, X_design, 'BMI_overweight')
        partial_shift = partial + np.mean(y)
        title_right = [
            f'log2FC={a["log2FC"]:+.3f}  [{a["log2FC_CI_lo"]:+.3f}, {a["log2FC_CI_hi"]:+.3f}]',
            f'p_OLS={a["p_ols_hc3"]:.4f}  p_limma={a["p_limma"]:.4f}  p_Wilcox={a["p_wilcoxon"]:.4f}',
        ]
        draw_box(axes[1], partial_shift, groups, palette,
                 ylabel='Partial residual\n(adj. age/GA/year)',
                 title_lines=title_right)
        plt.tight_layout()
        # 文件名安全处理
        safe = short.replace('/', '_').replace('α', 'alpha').replace(',', '_')
        single_path = FIG_DIR / f'candidate_{safe}.png'
        plt.savefig(single_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'  ✓ {single_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
