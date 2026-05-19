"""
阶段 1 预处理流程图 —— 紧凑版，适合 PPT 右下角
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(10, 3.2))
ax.set_xlim(0, 10)
ax.set_ylim(0, 3.2)
ax.axis('off')

# 颜色
c_main = '#4472C4'
c_sub  = '#5B9BD5'
c_out  = '#7030A0'

def box(ax, x, y, w, h, line1, line2='', color=c_main, fs1=8, fs2=7):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.01,rounding_size=0.15",
                         facecolor=color, edgecolor='black', linewidth=0.8)
    ax.add_patch(box)
    ax.text(x, y + 0.08, line1, ha='center', va='center', fontsize=fs1,
            color='white', weight='bold')
    if line2:
        ax.text(x, y - 0.12, line2, ha='center', va='center', fontsize=fs2,
                color='#E7E6E6', style='italic')

def arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='#333333', lw=1.0,
                               connectionstyle="arc3,rad=0"))

# 标题
ax.text(5, 2.9, '阶段 1：原始数据质量控制与预处理', ha='center', va='center',
        fontsize=11, weight='bold', color='#1F4E79')

# 7 个步骤横向排列
steps = [
    (0.9, 1.9, 1.0, 0.7, '1.1 QC 质控', 'CV% <= 30%'),
    (2.2, 1.9, 1.0, 0.7, '1.2 过滤', '80% 规则'),
    (3.5, 1.9, 1.0, 0.7, '1.3 插补', 'QRILC'),
    (4.8, 1.9, 1.0, 0.7, '1.4 归一化', '已 ng/g'),
    (6.1, 1.9, 1.0, 0.7, '1.5 转换', 'log2'),
    (7.4, 1.9, 1.0, 0.7, '1.6 缩放', 'Pareto'),
    (8.9, 1.9, 1.1, 0.7, '1.7 产物', '原始 / 分析', c_out),
]
for x, y, w, h, l1, l2, *c in steps:
    box(ax, x, y, w, h, l1, l2, color=c[0] if c else c_main)

# 箭头
for i in range(len(steps)-1):
    x1 = steps[i][0] + steps[i][2]/2 + 0.02
    x2 = steps[i+1][0] - steps[i+1][2]/2 - 0.02
    arrow(ax, x1, 1.9, x2, 1.9)

# 下方小字：双轨说明
ax.text(5, 0.9, '主轨 80% (n=67)  +  探索轨 50% (n=73)   |   单位: ng/g 组织湿重',
        ha='center', va='center', fontsize=8, color='#555555')

# 底部推荐组合
ax.text(5, 0.4, '推荐流程:  QRILC  →  log2  →  Pareto',
        ha='center', va='center', fontsize=8, color='#7F6000',
        bbox=dict(boxstyle='round,pad=0.25', facecolor='#FFF2CC', edgecolor='#D6B656', lw=0.8))

plt.tight_layout()
out_path = 'results/figures/stage1_flowchart_compact.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"紧凑流程图已保存至: {out_path}")
plt.close()
