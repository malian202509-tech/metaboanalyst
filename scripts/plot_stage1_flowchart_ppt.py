"""
阶段 1 预处理流程图 —— PPT 紧凑版，2×3 布局
配色参考：浅蓝 #B8D4E3 / 浅珊瑚 #E8A598
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(7.5, 3.6))
ax.set_xlim(0, 7.5)
ax.set_ylim(0, 3.6)
ax.axis('off')

# 配色
colors_top = dict(face='#B8D4E3', edge='#6B9BC2', text='#2F4F4F')
colors_bot = dict(face='#E8A598', edge='#C77D6B', text='#2F4F4F')

def draw_box(ax, x, y, w, h, line1, line2, colors):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.01,rounding_size=0.18",
                         facecolor=colors['face'], edgecolor=colors['edge'], linewidth=1.2)
    ax.add_patch(box)
    ax.text(x, y + 0.1, line1, ha='center', va='center', fontsize=10,
            color=colors['text'], weight='bold')
    if line2:
        ax.text(x, y - 0.14, line2, ha='center', va='center', fontsize=9,
                color=colors['text'])

def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='#444444', lw=1.3,
                               connectionstyle="arc3,rad=0"))

# 上排：QC / 过滤 / 插补
draw_box(ax, 1.5, 2.6, 1.6, 0.7, 'QC 质控', 'CV% <= 30%', colors_top)
draw_box(ax, 3.75, 2.6, 1.6, 0.7, '特征过滤', '80% 规则', colors_top)
draw_box(ax, 6.0, 2.6, 1.6, 0.7, '缺失值插补', 'QRILC', colors_top)

# 下排：归一化 / log2 / Pareto
draw_box(ax, 1.5, 1.1, 1.6, 0.7, '归一化', '已 ng/g', colors_bot)
draw_box(ax, 3.75, 1.1, 1.6, 0.7, '对数转换', 'log2', colors_bot)
draw_box(ax, 6.0, 1.1, 1.6, 0.7, '缩放', 'Pareto', colors_bot)

# 上排横向箭头
draw_arrow(ax, 2.3, 2.6, 2.95, 2.6)
draw_arrow(ax, 4.55, 2.6, 5.2, 2.6)

# 从插补向下到归一化（折线：向下 → 向左 → 向下）
verts = [(6.0, 2.25), (6.0, 1.8), (1.5, 1.8), (1.5, 1.45)]
codes = [mpatches.Path.MOVETO, mpatches.Path.LINETO,
         mpatches.Path.LINETO, mpatches.Path.LINETO]
path = mpatches.Path(verts, codes)
arrow_patch = FancyArrowPatch(path=path, arrowstyle='->', mutation_scale=15,
                              color='#444444', linewidth=1.3)
ax.add_patch(arrow_patch)

# 下排横向箭头
draw_arrow(ax, 2.3, 1.1, 2.95, 1.1)
draw_arrow(ax, 4.55, 1.1, 5.2, 1.1)

plt.tight_layout()
out_path = 'results/figures/stage1_flowchart_ppt.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"PPT 流程图已保存至: {out_path}")
plt.close()
