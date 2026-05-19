import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

OUTPUT_PATH = (
    "results/figures/06_Pathway_Schematics/"
    "DHA_derivatives_pathway_schematic.png"
)

def draw_box(ax, x, y, width, height, text, color, fontsize=10):
    """基于中心坐标(x,y)绘制一个圆角矩形框"""
    box = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width, height,
        boxstyle="round,pad=0.05,rounding_size=0.1",
        linewidth=1.2, edgecolor="#555555", facecolor=color, zorder=3
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color="#111111", zorder=4)
    return box

def draw_ortho_arrow(ax, x1, y1, x2, y2, y_turn):
    """绘制带转折的正交箭头 (下->平->下的路径)"""
    ax.plot([x1, x1], [y1, y_turn], color="#555555", lw=1.5, zorder=1)
    ax.plot([x1, x2], [y_turn, y_turn], color="#555555", lw=1.5, zorder=1)
    ax.annotate("", xy=(x2, y2), xytext=(x2, y_turn),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)

def draw_vert_arrow(ax, x, y1, y2):
    """绘制垂直向下的箭头"""
    ax.annotate("", xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)

def draw_horiz_arrow(ax, x1, y1, x2):
    """绘制水平方向的箭头 (从x1指向x2)"""
    ax.annotate("", xy=(x2, y1), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)

def main():
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")

    # 1. 标题
    ax.text(6, 7.7, "DHA 氧化代谢主要途径及衍生物上游倾向判别", ha="center", va="center", fontsize=15, fontweight="bold")

    # 2. 顶层中心节点：DHA
    draw_box(ax, 6, 6.8, 2.5, 0.7, "DHA\n(二十二碳六烯酸)", "#C9D7F8", 12)

    # 3. 中间层（三个主要途径）
    y_header = 4.8
    draw_box(ax, 2.5, y_header, 2.6, 0.6, "12/15-LOX 酶促倾向", "#F7C1C1", 11)
    draw_box(ax, 6.0, y_header, 2.6, 0.6, "CYP ω-羟化酶 倾向", "#C5E6C5", 11)
    draw_box(ax, 9.5, y_header, 3.4, 0.6, "自氧化 (非酶促) 及其它混合倾向", "#E2E2E2", 11)

    # 4. 从 DHA 到三个主要途径的正交连线
    y_turn = 5.6
    draw_ortho_arrow(ax, 6, 6.45, 2.5, 5.1, y_turn)
    draw_ortho_arrow(ax, 6, 6.45, 6.0, 5.1, y_turn)
    draw_ortho_arrow(ax, 6, 6.45, 9.5, 5.1, y_turn)

    # 5. 生成代谢物框级结构

    # === 左分支 (LOX) ===
    # 将其垂直排列在中心下方
    draw_box(ax, 2.5, 3.6, 1.8, 0.5, "17-HDHA", "#FCE2E2", 10)
    draw_box(ax, 2.5, 2.6, 1.8, 0.5, "14-HDHA", "#FCE2E2", 10)
    
    # 箭头连接
    draw_vert_arrow(ax, 2.5, 4.5, 3.85)  # LOX 到 17-HDHA
    draw_vert_arrow(ax, 2.5, 3.35, 2.85) # 17-HDHA 到 14-HDHA

    # === 中分支 (CYP) ===
    draw_box(ax, 6.0, 3.6, 1.8, 0.5, "20-HDHA", "#E1F4E1", 10)
    draw_vert_arrow(ax, 6.0, 4.5, 3.85)  # CYP 到 20-HDHA

    # === 右分支 (自氧化/混合) ===
    mx_items = ["4-HDHA", "7-HDHA", "8-HDHA", "10-HDHA", "11-HDHA", "13-HDHA", "16-HDHA"]
    # 采用左右对称两列的“树状枝干”布局
    spine_x = 9.5
    top_y = 4.5
    bottom_y = 1.6
    
    # 画一条贯穿中央的树干
    ax.plot([spine_x, spine_x], [top_y, bottom_y], color="#555555", lw=1.5, zorder=1)

    # 每一行的Y轴位置
    y_starts = [3.7, 2.9, 2.1, 1.4]
    box_w = 1.3
    
    for i, name in enumerate(mx_items):
        row = i // 2
        col = i % 2
        
        # 计算框的中心X坐标：左列在8.5，右列在10.5
        cx = 8.6 if col == 0 else 10.4
        cy = y_starts[row]
        
        draw_box(ax, cx, cy, box_w, 0.45, name, "#F2F2F2", 10)
        
        # 画横向的树枝及箭头 (从树干指向框的边缘)
        edge_x = cx + box_w/2 if col == 0 else cx - box_w/2
        draw_horiz_arrow(ax, spine_x, cy, edge_x)

    # 调整并保存输出
    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    print(f"Saved: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()
