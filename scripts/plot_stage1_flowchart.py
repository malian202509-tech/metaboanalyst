"""
绘制阶段 1：原始数据质量控制与预处理 流程图
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(12, 16))
ax.set_xlim(0, 12)
ax.set_ylim(0, 16)
ax.axis('off')

# 颜色定义
color_main = '#4472C4'      # 主步骤
color_decision = '#ED7D31'  # 决策/过滤
color_note = '#70AD47'      # 注释
color_output = '#7030A0'    # 产物

def draw_box(ax, x, y, w, h, text, color, fontsize=10, text_color='white'):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.02,rounding_size=0.2",
                         facecolor=color, edgecolor='black', linewidth=1.2)
    ax.add_patch(box)
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            color=text_color, weight='bold', wrap=True)
    return box

def draw_arrow(ax, x1, y1, x2, y2, color='black'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

def draw_note(ax, x, y, text, fontsize=9, color='#555555'):
    ax.text(x, y, text, ha='left', va='center', fontsize=fontsize,
            color=color, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#F2F2F2', edgecolor='#CCCCCC', lw=0.8))

# 标题
ax.text(6, 15.5, '阶段 1：原始数据质量控制与预处理', ha='center', va='center',
        fontsize=16, weight='bold', color='#1F4E79')
ax.text(6, 15.0, '目标：在进入统计模型前，将数据矩阵处理到"残差近似正态、特征尺度可比、低质量特征已剔除"的状态',
        ha='center', va='center', fontsize=10, color='#555555')

# 各步骤位置
steps = [
    (6, 13.8, 4.0, 0.7, '开始：原始数据矩阵', '#2E75B6', 11),
    (6, 12.6, 5.0, 0.9, '1.1  QC 样本质控', color_main, 11),
    (6, 11.2, 5.0, 0.9, '1.2  特征级缺失与零值过滤', color_main, 11),
    (6, 9.8, 5.0, 0.9, '1.3  缺失值插补', color_main, 11),
    (6, 8.4, 5.0, 0.9, '1.4  归一化（Normalization）', color_main, 11),
    (6, 7.0, 5.0, 0.9, '1.5  对数转换（Log transformation）', color_main, 11),
    (6, 5.6, 5.0, 0.9, '1.6  缩放（Scaling）', color_main, 11),
    (6, 4.2, 5.0, 0.9, '1.7  阶段产物', color_output, 11),
]

for x, y, w, h, text, color, fs in steps:
    draw_box(ax, x, y, w, h, text, color, fs)

# 连接箭头
arrow_pairs = [
    (6, 13.45, 6, 13.05),
    (6, 12.15, 6, 11.65),
    (6, 10.75, 6, 10.25),
    (6, 9.35, 6, 8.85),
    (6, 7.95, 6, 7.45),
    (6, 6.55, 6, 6.05),
    (6, 5.15, 6, 4.65),
]
for x1, y1, x2, y2 in arrow_pairs:
    draw_arrow(ax, x1, y1, x2, y2)

# 1.1 注释
draw_note(ax, 9.0, 12.9, '- QC 样本 PCA 聚集度评估\n- 特征 CV% <= 30% 保留\n- TIC 随进样顺序稳定性', fontsize=9)
# 1.1 剔除分支
ax.annotate('', xy=(3.2, 12.6), xytext=(3.5, 12.6),
            arrowprops=dict(arrowstyle='->', color=color_decision, lw=1.5))
ax.text(2.4, 12.6, 'CV% > 30%\n剔除特征', ha='center', va='center', fontsize=9, color=color_decision, weight='bold')

# 1.2 注释
draw_note(ax, 9.0, 11.5, '- 80% 规则 (modified 80% rule)\n- 任一组非零检出率 >= 80% 保留', fontsize=9)
# 1.2 剔除分支
ax.annotate('', xy=(3.2, 11.2), xytext=(3.5, 11.2),
            arrowprops=dict(arrowstyle='->', color=color_decision, lw=1.5))
ax.text(2.4, 11.2, '两组均 < 80%\n剔除特征', ha='center', va='center', fontsize=9, color=color_decision, weight='bold')

# 1.3 注释
draw_note(ax, 9.0, 10.1, '- 首选：半最小值 (1/2 min)\n- 备选：KNN / missForest', fontsize=9)

# 1.4 注释
draw_note(ax, 9.0, 8.7, '- 首选：PQN\n- 组织样本：归一化到蛋白浓度', fontsize=9)

# 1.5 注释
draw_note(ax, 9.0, 7.3, '- 推荐：log2 或 log10\n- 含零值：glog (generalized log)', fontsize=9)

# 1.6 注释
draw_note(ax, 9.0, 5.9, '- 推荐：Pareto scaling\n- 兼顾大/小丰度特征', fontsize=9)

# 1.7 产物分支
ax.text(3.5, 4.6, '原始矩阵\n(仅 1.1–1.3)', ha='center', va='center', fontsize=9, color=color_output, weight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#E2D4F0', edgecolor=color_output, lw=1.2))
ax.annotate('', xy=(3.5, 4.2), xytext=(3.5, 4.5),
            arrowprops=dict(arrowstyle='->', color=color_output, lw=1.5))

ax.text(3.5, 3.7, '分析矩阵\n(完整 1.1–1.6)', ha='center', va='center', fontsize=9, color=color_output, weight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#E2D4F0', edgecolor=color_output, lw=1.2))

# 典型组合注释（底部）
combo_box = FancyBboxPatch((1.0, 0.8), 10.0, 1.8,
                           boxstyle="round,pad=0.02,rounding_size=0.2",
                           facecolor='#FFF2CC', edgecolor='#D6B656', linewidth=1.5, linestyle='--')
ax.add_patch(combo_box)
ax.text(6, 2.1, '典型推荐组合', ha='center', va='center', fontsize=11, weight='bold', color='#7F6000')
ax.text(6, 1.5, '半最小值插补  →  PQN 归一化  →  log2 转换  →  Pareto scaling',
        ha='center', va='center', fontsize=10, color='#7F6000')
ax.text(6, 1.0, '注：文章 Methods 必须写明每一步的顺序与参数；原始矩阵与分析矩阵需分别归档',
        ha='center', va='center', fontsize=9, color='#595959', style='italic')

plt.tight_layout()
out_path = 'results/figures/stage1_qc_preprocessing_flowchart.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"流程图已保存至: {out_path}")
plt.close()
