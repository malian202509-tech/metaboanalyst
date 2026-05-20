import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

OUTPUT_PATH = (
    "results/figures/06_Pathway_Schematics/"
    "AA_pool_pathway_schematic.png"
)

def draw_box(ax, x, y, width, height, text, color, fontsize=9):
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

def draw_tree_column(ax, spine_x, top_y, bottom_y, items, color, box_w, box_h, y_spacing):
    """双排树状排列的函数"""
    ax.plot([spine_x, spine_x], [top_y, bottom_y], color="#555555", lw=1.5, zorder=1)
    for i, name in enumerate(items):
        row = i // 2
        col = i % 2
        # 计算框的中心X坐标 (增加偏移量，避免拥挤)
        cx = spine_x - box_w/2 - 0.35 if col == 0 else spine_x + box_w/2 + 0.35
        cy = top_y - 0.6 - row * y_spacing
        draw_box(ax, cx, cy, box_w, box_h, name, color, 9)
        # 画横向的树枝及箭头 (从树干指向框的边缘)
        edge_x = cx + box_w/2 if col == 0 else cx - box_w/2
        draw_horiz_arrow(ax, spine_x, cy, edge_x)

def draw_fishbone_column(ax, header_x, header_y_bottom, spine_offset, box_offset, groups, box_w, box_h, y_spacing):
    """单边鱼骨图的形式"""
    spine_x = header_x + spine_offset
    box_x = header_x + box_offset
    
    y_curr = header_y_bottom - 0.4
    
    # 模拟高度计算以确定底端
    spine_bottom = y_curr
    for grp in groups:
        if grp.get("label"): spine_bottom -= 0.5
        spine_bottom -= len(grp["items"]) * y_spacing
        spine_bottom -= 0.2
        
    spine_bottom += 0.2 + y_spacing
    
    # 直接绘制由左侧垂直向下主轴 (无起步折线)
    ax.plot([spine_x, spine_x], [header_y_bottom, spine_bottom], color="#555555", lw=1.5, zorder=1)

    real_y = y_curr
    for grp in groups:
        if grp.get("label"):
            ax.text(box_x, real_y, grp["label"], ha="center", va="center", fontsize=10, fontweight='bold', color=grp.get("label_color", "#333"))
            if grp.get("label_desc"):
                ax.text(box_x, real_y - 0.2, grp["label_desc"], ha="center", va="center", fontsize=8, color=grp.get("label_color", "#333"))
                real_y -= 0.2
            real_y -= 0.45
            
        for item in grp["items"]:
            draw_box(ax, box_x, real_y, box_w, box_h, item, grp["color"], 9)
            edge_x = box_x - box_w/2 if box_x > spine_x else box_x + box_w/2
            draw_horiz_arrow(ax, spine_x, real_y, edge_x)
            real_y -= y_spacing
            
        real_y -= 0.2

def draw_single_column(ax, cx, top_y, items, color, box_w, box_h, y_spacing):
    """垂直单列排布并带有向下连线"""
    y_current = top_y
    for i, name in enumerate(items):
        draw_box(ax, cx, y_current, box_w, box_h, name, color, 9)
        if i < len(items) - 1:
            draw_vert_arrow(ax, cx, y_current - box_h/2, y_current - y_spacing + box_h/2)
        y_current -= y_spacing

def draw_group_label(ax, x, y, text, color):
    ax.text(x, y, text, ha="center", va="center", fontsize=10, fontweight='bold', color=color)

def main():
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # 1. 标题
    ax.text(9, 9.5, "AA 氧化代谢的主要途径", ha="center", va="center", fontsize=16, fontweight="bold")

    # 2. 顶层中心节点：AA
    draw_box(ax, 9, 8.5, 3.0, 0.7, "AA\n(花生四烯酸)", "#C9D7F8", 12)

    # 3. 中间层分支（五大途径）
    y_header = 6.8
    branch_xs = [2.2, 5.8, 9.4, 13.0, 16.2]
    
    # 画框
    draw_box(ax, branch_xs[0], y_header, 2.8, 0.6, "COX 通路", "#Fbb4ae", 11)
    draw_box(ax, branch_xs[1], y_header, 2.8, 0.6, "LOX 通路", "#B3de69", 11)
    draw_box(ax, branch_xs[2], y_header, 3.2, 0.6, "CYP-Epo / sEH 通路", "#80b1d3", 11)
    draw_box(ax, branch_xs[3], y_header, 2.8, 0.6, "CYP ω-羟化酶 通路", "#Bee3f8", 11)
    draw_box(ax, branch_xs[4], y_header, 2.8, 0.6, "非酶促/内源大麻素", "#Fdbf6f", 11)

    # 连接母体与各个分支头
    y_turn = 7.6
    for x in branch_xs:
        draw_ortho_arrow(ax, 9, 8.15, x, 7.1, y_turn)

    # ================== 4. 具体代谢物排布 ==================

    box_h = 0.45
    
    # 4.1 COX 分支 (数量最多用双行树状)
    cox_items = ["TXB2", "PGE2", "PGD2", "PGF2α", "PGA2", "PGJ2", 
                 "6,15-diketo-dh-PGF1α", "11β-PGE2", "15d-PGJ2-like", 
                 "15-keto-PGE2", "13,14-dh-15k-PGE2", "12-HHT"]
    draw_tree_column(ax, branch_xs[0], 6.5, 1.2, cox_items, "#Fce2e2", 1.5, box_h, 0.8)

    # 4.2 LOX 分支 (分 12-, 15-, 5- 三组串联)
    grp_lox = [
        {"label": "12-LOX", "items": ["12-HETE", "12-HETrE"], "color": "#E1f4e1", "label_color": "#33691e"},
        {"label": "15-LOX", "items": ["15-HETE"], "color": "#E1f4e1", "label_color": "#33691e"},
        {"label": "5-LOX", "items": ["LTB4", "LTE4"], "color": "#E1f4e1", "label_color": "#33691e"}
    ]
    draw_fishbone_column(ax, branch_xs[1], 6.5, -1.1, 0.3, grp_lox, 1.8, box_h, 0.6)

    # 4.3 CYP-Epo -> sEH
    grp_cyp = [
        {"label": "CYP-Epo", "items": ["8,9-EpETrE", "11,12-EpETrE"], "color": "#D5e8f4", "label_color": "#0d47a1"},
        {"label": "↓ sEH 水解", "items": ["5,6-DiHETrE", "8,9-DiHETrE", "11,12-DiHETrE", "14,15-DiHETrE"], "color": "#E8e6f4", "label_color": "#4527a0"}
    ]
    draw_fishbone_column(ax, branch_xs[2], 6.5, -1.2, 0.4, grp_cyp, 2.2, box_h, 0.55)

    # 4.4 CYP ω-羟化酶
    grp_cyp4 = [
        {"items": ["16-HETE", "18-HETE"], "color": "#Ddf1fa"}
    ]
    draw_fishbone_column(ax, branch_xs[3], 6.5, -1.1, 0.3, grp_cyp4, 1.8, box_h, 0.8)

    # 4.5 自氧化及其他
    grp_auto = [
        {"label": "非酶促自氧化", "items": ["11-HETE", "8-HETE", "8-HETrE"], "color": "#Ffe6cc", "label_color": "#e65100"},
        {"label": "内源大麻素", "items": ["AEA"], "color": "#Ffffd9", "label_color": "#9e9d24"}
    ]
    draw_fishbone_column(ax, branch_xs[4], 6.5, -1.1, 0.3, grp_auto, 1.8, box_h, 0.65)


    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()