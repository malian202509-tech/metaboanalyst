import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

OUTPUT_PATH = "results/figures/06_Pathway_Schematics/AA_pool_pathway_schematic_horizontal.png"

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

def draw_ortho_arrow_horiz(ax, x1, y1, x2, y2, x_turn):
    """横向出发，折转到垂直，再横向指向目标"""
    ax.plot([x1, x_turn], [y1, y1], color="#555555", lw=1.5, zorder=1)
    ax.plot([x_turn, x_turn], [y1, y2], color="#555555", lw=1.5, zorder=1)
    ax.annotate("", xy=(x2, y2), xytext=(x_turn, y2),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)

def draw_tree_row(ax, spine_y, header_x_right, items, color, box_w, box_h, x_spacing):
    cols = (len(items) + 1) // 2
    start_x = header_x_right + 0.5
    end_x = start_x + (cols - 1) * x_spacing + 0.5

    ax.plot([header_x_right, end_x], [spine_y, spine_y], color="#555555", lw=1.5, zorder=1)

    for i, name in enumerate(items):
        col = i // 2
        row = i % 2 
        cx = start_x + col * x_spacing
        cy = spine_y + 0.85 if row == 0 else spine_y - 0.85
        draw_box(ax, cx, cy, box_w, box_h, name, color, 9)
        
        # 箭头：主轴单向指向框
        edge_y = cy - box_h/2 if row == 0 else cy + box_h/2
        ax.annotate("", xy=(cx, edge_y), xytext=(cx, spine_y),
                    arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)

def draw_fishbone_row(ax, header_y, header_x_right, groups, box_w, box_h, x_spacing):
    spine_y = header_y
    box_y = header_y - 1.2
    
    start_x = header_x_right + 0.5
    
    # 计算总延展长度
    x_curr = start_x
    for grp in groups:
        x_curr += len(grp["items"]) * x_spacing + 0.5
    spine_end = x_curr - 0.5 + 0.3
    
    # 绘制主结构直干 (不再折线，直接从上一个框的右边发出平线)
    ax.plot([header_x_right, spine_end], [spine_y, spine_y], color="#555555", lw=1.5, zorder=1)
    
    real_x = start_x
    for grp in groups:
        grp_start_x = real_x
        
        for item in grp["items"]:
            draw_box(ax, real_x, box_y, box_w, box_h, item, grp["color"], 9)
            
            # 垂直连接线从干向下挂框
            edge_y = box_y + box_h/2
            ax.annotate("", xy=(real_x, edge_y), xytext=(real_x, spine_y),
                        arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)
            real_x += x_spacing
            
        grp_end_x = real_x - x_spacing
        grp_center_x = (grp_start_x + grp_end_x) / 2
        
        # 标签放在方框的正上方 (脊柱和方框之间)
        if grp.get("label"):
            label_y = box_y + box_h/2 + 0.35 
            ax.text(grp_center_x, label_y, grp["label"], ha="center", va="center", fontsize=10, fontweight='bold', color=grp.get("label_color", "#333"))
        
        real_x += 0.5

def main():
    fig, ax = plt.subplots(figsize=(24, 14))
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 15)
    ax.axis("off")

    # 1. 标题
    ax.text(12, 14.5, "AA 氧化代谢的主要途径 (水平排布)", ha="center", va="center", fontsize=16, fontweight="bold")

    # 2. 顶级中心节点：AA
    draw_box(ax, 1.8, 7.5, 2.8, 0.8, "AA\n(花生四烯酸)", "#C9D7F8", 12)

    # 3. Y轴层次分布
    branch_ys = [12.5, 9.8, 7.1, 4.4, 1.7]
    branch_x = 5.6
    
    draw_box(ax, branch_x, branch_ys[0], 2.8, 0.6, "COX 通路", "#Fbb4ae", 11)
    draw_box(ax, branch_x, branch_ys[1], 2.8, 0.6, "LOX 通路", "#B3de69", 11)
    draw_box(ax, branch_x, branch_ys[2], 3.2, 0.6, "CYP-Epo / sEH 通路", "#80b1d3", 11)
    draw_box(ax, branch_x, branch_ys[3], 2.8, 0.6, "CYP ω-羟化酶 通路", "#Bee3f8", 11)
    draw_box(ax, branch_x, branch_ys[4], 2.8, 0.6, "非酶促/内源大麻素", "#Fdbf6f", 11)

    # 画引出折线
    x_turn = 3.6
    for y in branch_ys:
        draw_ortho_arrow_horiz(ax, 3.2, 7.5, branch_x - 1.6, y, x_turn)

    box_h = 0.5
    
    # ================= 4.1 COX 分支 =================
    cox_items = ["TXB2", "PGE2", "PGD2", "PGF2α", "PGA2", "PGJ2", 
                 "6,15-diketo-dh-\nPGF1α", "11β-PGE2", "15d-PGJ2-like", 
                 "15-keto-PGE2", "13,14-dh-15k-\nPGE2", "12-HHT"]
    draw_tree_row(ax, branch_ys[0], branch_x + 1.4, cox_items, "#Fce2e2", 2.2, box_h, 2.4)

    # ================= 4.2 LOX 分支 =================
    grp_lox = [
        {"label": "12-LOX", "items": ["12-HETE", "12-HETrE"], "color": "#E1f4e1", "label_color": "#33691e"},
        {"label": "15-LOX", "items": ["15-HETE"], "color": "#E1f4e1", "label_color": "#33691e"},
        {"label": "5-LOX", "items": ["LTB4", "LTE4"], "color": "#E1f4e1", "label_color": "#33691e"}
    ]
    draw_fishbone_row(ax, branch_ys[1], branch_x + 1.4, grp_lox, 2.0, box_h, 2.2)

    # ================= 4.3 CYP-Epo -> sEH =================
    grp_cyp = [
        {"label": "CYP-Epo", "items": ["8,9-EpETrE", "11,12-EpETrE"], "color": "#D5e8f4", "label_color": "#0d47a1"},
        {"label": "↓ sEH 水解", "items": ["5,6-DiHETrE", "8,9-DiHETrE", "11,12-DiHETrE", "14,15-DiHETrE"], "color": "#E8e6f4", "label_color": "#4527a0"}
    ]
    draw_fishbone_row(ax, branch_ys[2], branch_x + 1.6, grp_cyp, 2.1, box_h, 2.3)

    # ================= 4.4 CYP ω-羟化酶 =================
    grp_cyp4 = [
        {"items": ["16-HETE", "18-HETE"], "color": "#Ddf1fa"}
    ]
    draw_fishbone_row(ax, branch_ys[3], branch_x + 1.4, grp_cyp4, 2.0, box_h, 2.2)

    # ================= 4.5 自氧化及其他 =================
    grp_auto = [
        {"label": "非酶促自氧化", "items": ["11-HETE", "8-HETE", "8-HETrE"], "color": "#Ffe6cc", "label_color": "#e65100"},
        {"label": "内源大麻素", "items": ["AEA"], "color": "#Ffffd9", "label_color": "#9e9d24"}
    ]
    draw_fishbone_row(ax, branch_ys[4], branch_x + 1.4, grp_auto, 2.0, box_h, 2.2)
    
    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()