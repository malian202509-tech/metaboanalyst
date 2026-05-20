"""15f. AA 池通路全景示意图 — 参考 06_Pathway_Schematics 风格.

把 15e 的 AA 池 ANCOVA 全景结果绘制成发表级 schematic:
  顶部 AA 母体
  4 大酶系分支 (横向并列): COX / LOX / CYP-Epo→sEH / CYP4 ω
  每分支下垂直列出 panel 内核心成员, 标 log2FC + p
  底部独立分支: 非酶 auto-ox HETE + 内源大麻素 AEA

颜色: 通路色 (COX 橙红 / LOX 绿 / CYP-Epo 蓝 / sEH 紫 / CYP4 ω 浅蓝 / 非酶 橙 / AEA 黄绿)
方向: ↑/↓ 文字标注 + 显著性 (★ raw p<0.05, ⚠ 0.05≤p<0.10)

输出:
  results/figures/06_Pathway_Schematics/AA_pool_pathway_schematic.png
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / 'results/figures/06_Pathway_Schematics'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 通路色 (与 06 目录里 AA_LOX_pathway_schematic 风格一致, 浅色填充 + 深色边)
PATH_STYLE = {
    'AA':       {'fill': '#cfd8dc', 'edge': '#37474f'},
    'COX':      {'fill': '#fbb4ae', 'edge': '#b71c1c'},
    'LOX':      {'fill': '#b3de69', 'edge': '#33691e'},
    'CYP-Epo':  {'fill': '#80b1d3', 'edge': '#0d47a1'},
    'sEH':      {'fill': '#bebada', 'edge': '#4527a0'},
    'CYP4':     {'fill': '#bee3f8', 'edge': '#01579b'},
    'Auto-ox':  {'fill': '#fdbf6f', 'edge': '#e65100'},
    'AEA':      {'fill': '#ffffb3', 'edge': '#9e9d24'},
}


def sig_mark(p):
    if pd.isna(p):
        return ''
    if p < 0.05:
        return ' ★'
    if p < 0.10:
        return ' ⚠'
    return ''


def arrow_mark(log2fc):
    if pd.isna(log2fc):
        return '?'
    if log2fc > 0.05:
        return '↑'
    if log2fc < -0.05:
        return '↓'
    return '~'


def draw_box(ax, x, y, w, h, text, style, fontsize=9.5, fontweight='normal',
             text_color='black'):
    box = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                         boxstyle='round,pad=0.02,rounding_size=0.08',
                         linewidth=1.3, facecolor=style['fill'],
                         edgecolor=style['edge'])
    ax.add_patch(box)
    ax.text(x, y, text, ha='center', va='center',
            fontsize=fontsize, fontweight=fontweight, color=text_color,
            wrap=True, linespacing=1.15)


def draw_arrow(ax, x1, y1, x2, y2, color='#37474f', lw=1.2,
               style='->,head_width=4,head_length=5'):
    arr = FancyArrowPatch((x1, y1), (x2, y2),
                          arrowstyle=style,
                          color=color, linewidth=lw,
                          mutation_scale=1, shrinkA=2, shrinkB=2)
    ax.add_patch(arr)


def load_data():
    df = pd.read_csv(ROOT / 'results/tables/aa_pool_pathway_landscape.csv')
    out = {}
    for _, r in df.iterrows():
        out[r['short_label']] = {
            'log2FC': r['log2FC'], 'p': r['p_ols_hc3'],
            'p_limma': r['p_limma'], 'p_wil': r['p_wilcoxon'],
            'pathway': r['pathway'], 'in80': r['in_main_80'],
            'src': r['source_track'],
        }
    return out


def label(d, name, show_p='ols'):
    if name not in d:
        return f"{name}\n(无数据)"
    r = d[name]
    lfc = r['log2FC']
    p = r['p'] if show_p == 'ols' else r['p_limma']
    mark = sig_mark(p)
    arr = arrow_mark(lfc)
    track = ' [50]' if r['src'] == '50_explor' else ''
    return f"{name}{track}\n{arr} log2FC={lfc:+.2f}  p={p:.3f}{mark}"


def main():
    print('=== 15f. AA 池通路 schematic ===\n')
    d = load_data()

    fig, ax = plt.subplots(figsize=(17, 12))
    ax.set_xlim(0, 22)
    ax.set_ylim(0, 16.5)
    ax.axis('off')

    # ---- 标题 ----
    ax.text(11, 16.0,
            'AA 池通路全景：花生四烯酸经 COX / LOX / CYP-Epo / sEH / CYP4-ω '
            '4 大酶系生成的主要代谢物',
            ha='center', fontsize=14, fontweight='bold')
    ax.text(11, 15.5,
            '(超重肥胖 vs 正常, n=165, ANCOVA 调整 age/GA/year/batch; ★ raw p<0.05; ⚠ 0.05≤p<0.10)',
            ha='center', fontsize=10.5, color='#444444')

    # ---- AA 母体 ----
    aa_text = label(d, 'AA').replace('AA\n', 'AA (花生四烯酸)\n')
    draw_box(ax, 11, 14.2, 4.0, 1.2, aa_text, PATH_STYLE['AA'],
             fontsize=10.5, fontweight='bold')

    # ---- 4 大分支酶系顶层标题框 ----
    branch_xs = {'COX': 3.0, 'LOX': 8.5, 'CYP-Epo': 14.0, 'CYP4': 19.0}
    branch_titles = {
        'COX':     'COX 通路\n(PG / TX / PGI2 / HHT)',
        'LOX':     'LOX 通路\n(HETE / HETrE / LT)',
        'CYP-Epo': 'CYP-Epo → sEH 串联\n(EpETrE → DiHETrE)',
        'CYP4':    'CYP4 ω-羟化\n(ω-HETE)',
    }
    for k, x in branch_xs.items():
        draw_box(ax, x, 12.4, 3.6, 1.0, branch_titles[k], PATH_STYLE[k],
                 fontsize=10, fontweight='bold')
        # AA → 各分支
        draw_arrow(ax, 11, 13.55, x, 12.95, color=PATH_STYLE[k]['edge'], lw=1.5)

    # =================== COX 分支 (左) ===================
    cox_x = branch_xs['COX']
    cox_items = [
        ('TXB2',                    11.4),
        ('PGE2',                    10.5),
        ('PGD2',                     9.6),
        ('PGF2α',                    8.7),
        ('PGA2',                     7.8),
        ('PGJ2',                     6.9),
        ('6,15-diketo-dh-PGF1α',     5.7),  # PGI2 metabolite, 占两行
        ('11β-PGE2',                 4.7),
        ('15d-PGJ2-like',            3.8),
        ('15-keto-PGE2',             2.9),
        ('13,14-dh-15k-PGE2',        2.0),
        ('12-HHT',                   1.1),
    ]
    for name, y in cox_items:
        # 6,15-diketo / 13,14-dh 名字太长用 narrow box
        w = 4.0 if len(name) <= 12 else 4.4
        draw_box(ax, cox_x, y, w, 0.78, label(d, name), PATH_STYLE['COX'],
                 fontsize=8.5)

    # =================== LOX 分支 (中左) ===================
    lox_x = branch_xs['LOX']
    # 12-LOX 子簇
    ax.text(lox_x - 1.6, 11.4, '12-LOX', fontsize=9, fontweight='bold',
            ha='right', color=PATH_STYLE['LOX']['edge'])
    draw_box(ax, lox_x, 11.4, 3.6, 0.78, label(d, '12-HETE'),
             PATH_STYLE['LOX'], fontsize=8.5)
    draw_box(ax, lox_x, 10.4, 3.6, 0.78, label(d, '12-HETrE'),
             PATH_STYLE['LOX'], fontsize=8.5, fontweight='bold')
    # 15-LOX 子簇
    ax.text(lox_x - 1.6, 9.0, '15-LOX', fontsize=9, fontweight='bold',
            ha='right', color=PATH_STYLE['LOX']['edge'])
    draw_box(ax, lox_x, 9.0, 3.6, 0.78, label(d, '15-HETE'),
             PATH_STYLE['LOX'], fontsize=8.5)
    # 5-LOX 子簇 (50 轨)
    ax.text(lox_x - 1.6, 7.6, '5-LOX [50]', fontsize=9, fontweight='bold',
            ha='right', color=PATH_STYLE['LOX']['edge'])
    draw_box(ax, lox_x, 7.6, 3.6, 0.78, label(d, 'LTB4'),
             PATH_STYLE['LOX'], fontsize=8.5)
    draw_box(ax, lox_x, 6.6, 3.6, 0.78, label(d, 'LTE4'),
             PATH_STYLE['LOX'], fontsize=8.5)
    # 非酶 HETE 子簇 (放在 LOX 列底部，但用 Auto-ox 颜色区分)
    ax.text(lox_x - 1.6, 5.0, '非酶氧化\n(11R/8-LOX\n人组织罕见)',
            fontsize=8.5, fontweight='bold', ha='right',
            color=PATH_STYLE['Auto-ox']['edge'])
    draw_box(ax, lox_x, 5.0, 3.6, 0.78, label(d, '11-HETE'),
             PATH_STYLE['Auto-ox'], fontsize=8.5)
    draw_box(ax, lox_x, 4.0, 3.6, 0.78, label(d, '8-HETE'),
             PATH_STYLE['Auto-ox'], fontsize=8.5)
    # 8-LOX HETrE 副产物
    ax.text(lox_x - 1.6, 2.6, '8-LOX/auto\nvia HpETrE', fontsize=8.5,
            fontweight='bold', ha='right',
            color=PATH_STYLE['Auto-ox']['edge'])
    draw_box(ax, lox_x, 2.6, 3.6, 0.78, label(d, '8-HETrE'),
             PATH_STYLE['Auto-ox'], fontsize=8.5)

    # =================== CYP-Epo + sEH 分支 (中右) ===================
    cyp_x = branch_xs['CYP-Epo']
    # CYP-Epo 顶部
    ax.text(cyp_x - 1.9, 11.4, 'CYP-Epo\n(CYP2C/2J)', fontsize=9,
            fontweight='bold', ha='right',
            color=PATH_STYLE['CYP-Epo']['edge'])
    draw_box(ax, cyp_x, 11.4, 3.6, 0.78, label(d, '8,9-EpETrE'),
             PATH_STYLE['CYP-Epo'], fontsize=8.5)
    draw_box(ax, cyp_x, 10.4, 3.6, 0.78, label(d, '11,12-EpETrE'),
             PATH_STYLE['CYP-Epo'], fontsize=8.5)
    # panel 缺失提示
    ax.text(cyp_x, 9.4, '(panel 缺 5,6- / 14,15-EpETrE)',
            ha='center', fontsize=8.2, style='italic', color='#555555')
    # sEH 标签 + 箭头
    ax.annotate('', xy=(cyp_x, 8.3), xytext=(cyp_x, 9.1),
                arrowprops=dict(arrowstyle='->', lw=1.5,
                                color=PATH_STYLE['sEH']['edge']))
    ax.text(cyp_x + 1.9, 8.7, 'sEH 水解', fontsize=9, fontweight='bold',
            color=PATH_STYLE['sEH']['edge'])
    # sEH 产物 4 个 DiHETrE
    draw_box(ax, cyp_x, 7.8, 3.6, 0.78, label(d, '5,6-DiHETrE'),
             PATH_STYLE['sEH'], fontsize=8.5)
    draw_box(ax, cyp_x, 6.8, 3.6, 0.78, label(d, '8,9-DiHETrE'),
             PATH_STYLE['sEH'], fontsize=8.5)
    draw_box(ax, cyp_x, 5.8, 3.6, 0.78, label(d, '11,12-DiHETrE'),
             PATH_STYLE['sEH'], fontsize=8.5)
    draw_box(ax, cyp_x, 4.8, 3.6, 0.78, label(d, '14,15-DiHETrE'),
             PATH_STYLE['sEH'], fontsize=8.5)

    # =================== CYP4 ω-羟化 (右) ===================
    cyp4_x = branch_xs['CYP4']
    ax.text(cyp4_x - 1.9, 11.4, 'CYP4F/4A', fontsize=9,
            fontweight='bold', ha='right',
            color=PATH_STYLE['CYP4']['edge'])
    draw_box(ax, cyp4_x, 11.4, 3.6, 0.78, label(d, '16-HETE'),
             PATH_STYLE['CYP4'], fontsize=8.5)
    draw_box(ax, cyp4_x, 10.4, 3.6, 0.78, label(d, '18-HETE'),
             PATH_STYLE['CYP4'], fontsize=8.5)
    # 提示：CYP4 在 AA 池反向, 在 DHA 池 ↑↑
    note_box = FancyBboxPatch((cyp4_x - 2.0, 7.4), 4.0, 2.4,
                              boxstyle='round,pad=0.05,rounding_size=0.08',
                              linewidth=1.0, facecolor='#fffde7',
                              edgecolor='#f57f17', linestyle='--')
    ax.add_patch(note_box)
    ax.text(cyp4_x, 8.6,
            '★ 关键对比\nAA→ω-HETE: 0/2 ↑ (反向 ↓)\n'
            'DHA→20-HDoHE: ↑↑ p=0.022\n→ CYP4 偏好 DHA 而非 AA',
            ha='center', va='center', fontsize=9, color='#bf360c',
            fontweight='bold')

    # =================== AEA (底部右) ===================
    aea_y = 1.5
    aea_x = 19.5
    ax.text(aea_x - 1.9, aea_y + 0.5, 'NAPE-PLD', fontsize=9, fontweight='bold',
            ha='right', color=PATH_STYLE['AEA']['edge'])
    draw_box(ax, aea_x, aea_y, 3.6, 1.0, label(d, 'AEA'),
             PATH_STYLE['AEA'], fontsize=8.8)
    # AA → AEA 分支 (从 AA 母体单独发一根虚线箭头)
    aa_to_aea = FancyArrowPatch((11, 13.7), (aea_x, aea_y + 0.6),
                                arrowstyle='->,head_width=3,head_length=4',
                                color=PATH_STYLE['AEA']['edge'],
                                linewidth=1.0, linestyle=':',
                                connectionstyle='arc3,rad=-0.3')
    ax.add_patch(aa_to_aea)
    ax.text(15.5, 7.5, 'AA + 乙醇胺\n→ AEA (内源大麻素)',
            fontsize=8.5, color=PATH_STYLE['AEA']['edge'],
            style='italic', ha='center')

    # =================== 底部图例 / 关键发现汇总 ===================
    summary_box = FancyBboxPatch((0.5, 0.05), 17.5, 0.95,
                                  boxstyle='round,pad=0.05,rounding_size=0.1',
                                  linewidth=1.0, facecolor='#f5f5f5',
                                  edgecolor='#666666')
    ax.add_patch(summary_box)
    ax.text(9.25, 0.52,
            '关键发现 (AA 池):  AA 母体 ↓10%  |  COX 15/16↑ (TXB2 ★ p=0.045, PGI2-met +0.64)  |  '
            'LOX 5/6↑ (12-HETrE ★★ p=0.021)  |  CYP-Epo flat  |  sEH 异构体特异 (5,6↑ vs 8,9↓ p=0.055)  |  '
            'CYP4 ω 反向 ↓ (偏好 DHA)  |  非酶 11/8-HETE ↑ (11-HETE ⚠ p=0.067)  |  AEA ↓',
            ha='center', va='center', fontsize=8.8)

    out_fig = OUT_DIR / 'AA_pool_pathway_schematic.png'
    plt.savefig(out_fig, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  ✓ {out_fig.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
