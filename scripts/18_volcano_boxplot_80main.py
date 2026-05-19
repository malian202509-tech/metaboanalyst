"""18. 80 主分析差异代谢物火山图 + 候选箱线图 (论文版).

依据 SOP §6.2 + 论文 Figure 标准:
  火山图 (Volcano plot):
    - x = log2FC, y = -log10(p_limma)
    - 候选 (★★ Core, p_limma<0.05 + |FC|≥log2(1.2) + 离群稳健) 用红/蓝实心 + 名字标注
    - 阈值参考线: x=±log2(1.2), y=-log10(0.05)
    - 灰色 = 非显著, 浅红/蓝 = 仅满足一个阈值
    - 右下角图例 + 全集 n / 候选 n 注释

  候选箱线图 (Boxplot):
    - 80 主轨 2 个候选: 12-HETrE, 20-HDoHE
    - 每候选 1 列, 共 1×2 网格
    - 原始 log2 浓度 + 散点 jitter overlay
    - 显著性 bracket (基于 p_limma): *** p<0.001, ** p<0.01, * p<0.05, NS p≥0.05
    - 副标题: log2FC + 95% CI + 三检验 p 值

输入:
  results/tables/ancova_main_80.csv
  results/tables/diff_candidates_80.csv
  data/03_batch_corrected/ori_n165_filtered80_log2_combat.xlsx
  data/02_preprocessed/sample_alignment_n165.csv
  data/02_preprocessed/metabolite_family_map.csv (取 short_label)

输出:
  results/figures/02_Volcano_plots/volcano_80main.png (论文 Figure 2 候选)
  results/figures/04_Metabolite_boxplots/boxplot_80main_candidates.png (论文 Figure 3 候选)
  results/figures/04_Metabolite_boxplots/boxplot_80main_{short}.png (单张)
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'  # 数学符号用普通字体

ROOT = Path(__file__).resolve().parent.parent
COMBAT_DIR = ROOT / 'data' / '03_batch_corrected'
PREP_DIR = ROOT / 'data' / '02_preprocessed'
TABLES = ROOT / 'results' / 'tables'
FIG_VOLC = ROOT / 'results' / 'figures' / '02_Volcano_plots'
FIG_BOX = ROOT / 'results' / 'figures' / '04_Metabolite_boxplots'
FIG_VOLC.mkdir(parents=True, exist_ok=True)
FIG_BOX.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'
BMI_COLORS = {'正常': '#67A9CF', '超重肥胖': '#EF8A62'}

FC_THRESH = np.log2(1.2)
P_THRESH = 0.05


def stars(p):
    if pd.isna(p): return 'NS'
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'NS'


def draw_volcano(anc, cand_set, short_map, out_path):
    """火山图 — 论文期刊标准 (三色简洁版).
    分类: down-regulated (蓝) / not-significant (灰) / up-regulated (红)
    阈值: |log2FC| ≥ log2(1.2)  AND  p_limma < 0.05
    """
    x = anc['log2FC'].values
    y = -np.log10(anc['p_limma'].clip(1e-10).values)
    fc = anc['log2FC'].values
    p = anc['p_limma'].values

    sig_up = (fc > FC_THRESH) & (p < P_THRESH)
    sig_dn = (fc < -FC_THRESH) & (p < P_THRESH)
    ns = ~(sig_up | sig_dn)

    fig, ax = plt.subplots(figsize=(9.0, 7.0))

    # 三色 (期刊标准: 蓝/灰/红)
    ax.scatter(x[ns], y[ns], s=30, c='#BFBFBF', alpha=0.75,
               edgecolors='none', label='not-significant')
    ax.scatter(x[sig_dn], y[sig_dn], s=34, c='#377EB8', alpha=0.85,
               edgecolors='none', label='down-regulated')
    ax.scatter(x[sig_up], y[sig_up], s=34, c='#E41A1C', alpha=0.85,
               edgecolors='none', label='up-regulated')

    # 仅标注 up-regulated / down-regulated 的显著点 (期刊风, 不标灰色 NS 点)
    up_idx = np.where(sig_up)[0].tolist()
    dn_idx = np.where(sig_dn)[0].tolist()

    # 按 -log10 p 降序排; 同象限多点垂直错开
    up_idx.sort(key=lambda i: -y[i])
    dn_idx.sort(key=lambda i: -y[i])

    def _place(idx_list, side):
        """side: 'right' / 'left' 决定文字朝向."""
        # 同象限多点用阶梯式 dy 错开
        for k, i in enumerate(idx_list):
            name = anc['Metabolite Name'].iloc[i]
            short = short_map.get(name, name[:25])
            col = '#E41A1C' if x[i] > 0 else '#377EB8'
            if side == 'right':
                dx = 22
                dy = 18 - k * 18   # 第 1 个向上, 第 2 个向下, 依次错开
                ha = 'left'
            else:
                dx = -22
                dy = 18 - k * 18
                ha = 'right'
            ax.annotate(short, (x[i], y[i]),
                        xytext=(dx, dy), textcoords='offset points',
                        fontsize=10, color=col,
                        ha=ha, va='center',
                        arrowprops=dict(arrowstyle='-', color=col,
                                        lw=0.6, alpha=0.85,
                                        shrinkA=0, shrinkB=3))

    _place(up_idx, 'right')
    _place(dn_idx, 'left')

    # 阈值参考线 (黑色虚线, 期刊风)
    ax.axhline(-np.log10(P_THRESH), color='black', linestyle='--', lw=0.8, alpha=0.85)
    ax.axvline(FC_THRESH, color='black', linestyle='--', lw=0.8, alpha=0.85)
    ax.axvline(-FC_THRESH, color='black', linestyle='--', lw=0.8, alpha=0.85)

    # 阈值文本 (淡灰色小字, 仅作标记)
    xmax_pre = max(abs(x.min()), abs(x.max())) * 1.15
    ymax_pre = float(y.max()) * 1.30
    ax.text(xmax_pre * 0.98, -np.log10(P_THRESH) + ymax_pre * 0.012,
            r'$p$ = 0.05', fontsize=8.5, color='#666666',
            ha='right', va='bottom')
    ax.text(FC_THRESH + xmax_pre * 0.008, ymax_pre * 0.02,
            r'$\mathrm{log_2}$(1.2)', fontsize=8.5, color='#666666',
            ha='left', va='bottom')
    ax.text(-FC_THRESH - xmax_pre * 0.008, ymax_pre * 0.02,
            r'$-\mathrm{log_2}$(1.2)', fontsize=8.5, color='#666666',
            ha='right', va='bottom')

    # 坐标轴 (mathtext 渲染下标, 不依赖 unicode 字符字形)
    ax.set_xlabel(r'$\mathrm{log_2 Fold\ Change}$', fontsize=13)
    ax.set_ylabel(r'$-\mathrm{log_{10}}\,(P\,\mathrm{value})$', fontsize=13)

    # 图例 (右侧, 标题 'Group')
    leg = ax.legend(title='Group', loc='center left', bbox_to_anchor=(1.01, 0.5),
                    fontsize=11, title_fontsize=12, frameon=False)

    # 坐标轴范围: 对称, 留 25% 顶部给标签
    xmax = max(abs(x.min()), abs(x.max())) * 1.15
    ymax = float(y.max()) * 1.30
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-0.05, ymax)

    # 封闭式四边框
    for s in ['top', 'bottom', 'left', 'right']:
        ax.spines[s].set_visible(True)
        ax.spines[s].set_linewidth(1.0)
        ax.spines[s].set_color('black')
    ax.tick_params(direction='out', length=4, width=1.0, color='black')

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {out_path.relative_to(ROOT)}')


def draw_boxplot_with_bracket(ax, y_vals, groups, p_anc, p_lim, p_wil,
                               fc, ci_lo, ci_hi, short, family):
    """单候选 boxplot + 散点 + 显著性 bracket."""
    palette = BMI_COLORS
    grp_order = ['正常', '超重肥胖']
    data = [y_vals[groups == g] for g in grp_order]

    bp = ax.boxplot(data, positions=[0, 1], widths=0.55, patch_artist=True,
                    boxprops=dict(linewidth=1.2),
                    medianprops=dict(color='black', linewidth=1.8),
                    whiskerprops=dict(linewidth=1.2),
                    capprops=dict(linewidth=1.2),
                    flierprops=dict(marker='', markersize=0))
    for b, g in zip(bp['boxes'], grp_order):
        b.set_facecolor(palette[g])
        b.set_alpha(0.55)
        b.set_edgecolor(palette[g])

    rng = np.random.default_rng(42)
    for x_pos, g in enumerate(grp_order):
        yvals = y_vals[groups == g]
        x_jit = rng.normal(x_pos, 0.075, size=len(yvals))
        ax.scatter(x_jit, yvals, s=22, c=palette[g], edgecolors='black',
                   linewidths=0.45, alpha=0.7, zorder=3)

    # 显著性 bracket (基于 p_limma)
    y_max = float(np.max(y_vals))
    y_min = float(np.min(y_vals))
    span = y_max - y_min
    bar_y = y_max + span * 0.08
    tick = span * 0.025
    ax.plot([0, 0, 1, 1], [bar_y, bar_y + tick, bar_y + tick, bar_y],
            color='black', lw=1.2)
    sig_label = stars(p_lim)
    label_text = sig_label
    ax.text(0.5, bar_y + tick * 1.4, label_text,
            ha='center', va='bottom', fontsize=11,
            fontweight='bold' if sig_label != 'NS' else 'normal',
            color='#D62728' if sig_label != 'NS' else 'gray')

    # 上扩 y 留 bracket 空间
    ax.set_ylim(y_min - span * 0.05, y_max + span * 0.28)

    ax.set_xticks([0, 1])
    x_labels = [
        f'Normal\n(n={int((groups == grp_order[0]).sum())})',
        f'Overweight/Obese\n(n={int((groups == grp_order[1]).sum())})',
    ]
    ax.set_xticklabels(x_labels, fontsize=11)
    ax.set_ylabel(r'$\mathrm{log_2}$ (concentration)', fontsize=11)
    ax.set_title(short, fontsize=13, fontweight='bold', pad=8)
    ax.text(0.5, -0.18, r'* $P_{limma}$ < 0.05; ** $P_{limma}$ < 0.01',
            transform=ax.transAxes, ha='center', va='top', fontsize=9, color='#555555')
    ax.grid(axis='y', alpha=0.3, linestyle=':')
    # 期刊风: 隐藏顶/右框
    for s in ['top', 'bottom', 'left', 'right']:
        ax.spines[s].set_visible(True)
        ax.spines[s].set_linewidth(1.0)
        ax.spines[s].set_color('black')


def main():
    print('=== 18. 80 主分析火山图 + 候选箱线图 ===\n')

    anc = pd.read_csv(TABLES / 'ancova_main_80.csv', encoding='utf-8-sig')
    cand = pd.read_csv(TABLES / 'diff_candidates_80.csv', encoding='utf-8-sig')
    cand_set = set(cand['Metabolite Name'])
    fmap = pd.read_csv(PREP_DIR / 'metabolite_family_map.csv', encoding='utf-8-sig')
    short_map = dict(zip(fmap['Metabolite Name'], fmap['short_label']))
    family_map = dict(zip(fmap['Metabolite Name'], fmap['family_main']))
    print(f'  80 主轨: {len(anc)} 特征, {len(cand_set)} 候选 ★★')

    # === 火山图 ===
    print('\n[1] Volcano plot')
    draw_volcano(anc, cand_set, short_map, FIG_VOLC / 'volcano_80main.png')

    # === 候选 boxplot ===
    print('\n[2] Candidate boxplots')
    df = pd.read_excel(COMBAT_DIR / 'ori_n165_filtered80_log2_combat.xlsx')
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    align = pd.read_csv(PREP_DIR / 'sample_alignment_n165.csv')
    s2g = dict(zip(align['omx_id'], align['BMI_group']))
    groups = np.array([s2g.get(c, '?') for c in sample_cols])

    # 1×2 网格 + 单独 2 张
    cand_sorted = cand.sort_values('p_limma').reset_index(drop=True)
    n_cand = len(cand_sorted)
    fig, axes = plt.subplots(1, n_cand, figsize=(5.5 * n_cand, 6))
    if n_cand == 1:
        axes = [axes]

    for idx, r in cand_sorted.iterrows():
        name = r['Metabolite Name']
        y = df.loc[df['Metabolite Name'] == name, sample_cols].iloc[0].values.astype(float)
        short = short_map.get(name, name[:25])
        family = family_map.get(name, '-')
        draw_boxplot_with_bracket(
            axes[idx], y, groups,
            r['p_ols_hc3'], r['p_limma'], r['p_wilcoxon'],
            r['log2FC'], r['log2FC_CI_lo'], r['log2FC_CI_hi'],
            short, family,
        )

    fig.suptitle('Core differential oxylipins — 80% main analysis track (n=165)',
                 fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    grid_path = FIG_BOX / 'boxplot_80main_candidates.png'
    plt.savefig(grid_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {grid_path.relative_to(ROOT)}')

    # 每候选单独一张
    for idx, r in cand_sorted.iterrows():
        name = r['Metabolite Name']
        y = df.loc[df['Metabolite Name'] == name, sample_cols].iloc[0].values.astype(float)
        short = short_map.get(name, name[:25])
        family = family_map.get(name, '-')
        fig, ax = plt.subplots(figsize=(5.5, 6))
        draw_boxplot_with_bracket(
            ax, y, groups,
            r['p_ols_hc3'], r['p_limma'], r['p_wilcoxon'],
            r['log2FC'], r['log2FC_CI_lo'], r['log2FC_CI_hi'],
            short, family,
        )
        plt.tight_layout()
        safe = short.replace('/', '_').replace('α', 'alpha').replace(',', '_')
        single_path = FIG_BOX / f'boxplot_80main_{safe}.png'
        plt.savefig(single_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'  ✓ {single_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
