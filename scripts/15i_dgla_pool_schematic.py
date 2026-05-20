import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

OUTPUT_PATH = (
    "results/figures/06_Pathway_Schematics/"
    "DGLA_pool_pathway_schematic.png"
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

def draw_horiz_arrow(ax, x1, y1, x2, style="-|>"):
    """绘制水平方向的箭头(从x1指向x2)"""
    ax.annotate("", xy=(x2, y1), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color="#555555", lw=1.5, linestyle='--' if style=='->' else '-'), zorder=2)

def draw_tree_column(ax, spine_x, top_y, bottom_y, items, color, box_w, box_h, y_spacing):
    """双排树状排列的函数"""
    ax.plot([spine_x, spine_x], [top_y, bottom_y], color="#555555", lw=1.5, zorder=1)
    for i, name in enumerate(items):
        row = i // 2
        col = i % 2
        cx = spine_x - box_w/2 - 0.35 if col == 0 else spine_x + box_w/2 + 0.35
        cy = top_y - 0.6 - row * y_spacing
        draw_box(ax, cx, cy, box_w, box_h, name, color, 8.5)
        edge_x = cx + box_w/2 if col == 0 else cx - box_w/2
        draw_horiz_arrow(ax, spine_x, cy, edge_x)

def draw_fishbone_column(ax, header_x, header_y_bottom, spine_offset, box_offset, groups, box_w, box_h, y_spacing):
    """单边鱼骨图的形式"""
    spine_x = header_x + spine_offset
    box_x = header_x + box_offset
    y_curr = header_y_bottom - 0.4
    spine_bottom = y_curr
    for grp in groups:
        if grp.get("label"): spine_bottom -= 0.5
        spine_bottom -= len(grp["items"]) * y_spacing
        spine_bottom -= 0.2
    spine_bottom += 0.2 + y_spacing
    if spine_bottom > y_curr: spine_bottom = y_curr - 0.5
    ax.plot([spine_x, spine_x], [header_y_bottom, spine_bottom], color="#555555", lw=1.5, zorder=1)
    real_y = y_curr
    for grp in groups:
        if grp.get("label"):
            ax.text(box_x, real_y, grp["label"], ha="center", va="center", fontsize=10, fontweight="bold", color=grp.get("label_color", "#333"))
            real_y -= 0.45
        for item in grp["items"]:
            draw_box(ax, box_x, real_y, box_w, box_h, item, grp["color"], 9)
            edge_x = box_x - box_w/2 if box_x > spine_x else box_x + box_w/2
            draw_horiz_arrow(ax, spine_x, real_y, edge_x)
            real_y -= y_spacing
        real_y -= 0.2

def main():
    fig, ax = plt.subplots(figsize=(13, 9))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 10)
    ax.axis("off")

    ax.text(8, 9.5, "DGLA 氧化代谢的主要途径", ha="center", va="center", fontsize=16, fontweight="bold")

    dgla_x, dgla_y = 6.5, 8.5
    aa_x, aa_y = 11.5, 8.5
    draw_box(ax, dgla_x, dgla_y, 3.2, 0.7, "DGLA\n(双高-γ-亚麻酸)", "#F0D88C", 11)
    draw_box(ax, aa_x, aa_y, 2.6, 0.6, "AA (花生四烯酸)", "#D6E4F4", 10)
    
    ax.annotate("", xy=(aa_x - 1.3, aa_y), xytext=(dgla_x + 1.6, dgla_y),
                arrowprops=dict(arrowstyle="-|>", color="#888888", lw=1.2, linestyle="--"), zorder=2)
    ax.text((dgla_x + aa_x)/2, dgla_y + 0.35, "FADS1\n(Δ5-去饱和酶)", ha="center", va="bottom", fontsize=9, style="italic", color="#666666")

    y_header = 7.0
    branch_xs = [3.5, 8.0, 12.5]
    
    draw_box(ax, branch_xs[0], y_header, 3.2, 0.6, "COX 通路", "#F2B5B5", 11)
    draw_box(ax, branch_xs[1], y_header, 3.0, 0.6, "LOX 通路", "#C5E0C5", 11)
    draw_box(ax, branch_xs[2], y_header, 3.6, 0.6, "非酶促 (8-LOX/Auto)", "#F4C798", 11)

    y_turn = 7.6
    for x in branch_xs:
        draw_ortho_arrow(ax, dgla_x, 8.15, x, 7.3, y_turn)

    box_h = 0.45
    
    # 4.1 COX 分支
    cox_items = ["PGE1", "15-keto-PGE1", "PGD1", "13,14-dh-15k-PGE1", "PGF1α", "13,14-dh-15k-PGF1α"]
    draw_tree_column(ax, branch_xs[0], 6.7, 4.0, cox_items, "#FBE3E3", 2.1, box_h, 0.8)

    # 4.2 LOX 分支
    grp_lox = [
        {"label": "12-LOX", "items": ["12-HETrE"], "color": "#E8F3E8", "label_color": "#2F7A2F"},
        {"label": "15-LOX", "items": ["15-HETrE"], "color": "#E8F3E8", "label_color": "#2F7A2F"}
    ]
    draw_fishbone_column(ax, branch_xs[1], 6.7, -1.0, 0.4, grp_lox, 2.0, box_h, 0.7)

    # 4.3 Auto 分支
    grp_auto = [
        {"label": "8-LOX / Auto-ox", "items": ["8-HETrE"], "color": "#FBE5CC", "label_color": "#B85F1A"}
    ]
    draw_fishbone_column(ax, branch_xs[2], 6.7, -1.0, 0.4, grp_auto, 2.0, box_h, 0.7)

    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    print("Saved:", OUTPUT_PATH)

if __name__ == "__main__":
    main()