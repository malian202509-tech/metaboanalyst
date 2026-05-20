"""AA-LOX 通路示意图 (与 06_dha_pathway_schematic.py 同风格).

布局:
  顶层      AA (花生四烯酸)
  中层      5-LOX 倾向   |   12-LOX 倾向   |   15-LOX 倾向
  底层      5-HETE        |   12-HETE      |   15-HETE
            ↓             |    ↓           |    ↓
            LTA4          |   Hepoxilin    |   Lipoxin A4
            ↓             |                |   Lipoxin B4
            LTB4 / Cys-LTs|                |
配色逻辑 (按生物学倾向):
  5-LOX  (促炎 — 白三烯)         浅橙 #F8C9A0
  12-LOX (混合 — HETE/hepoxilin) 浅黄 #F8E5A8
  15-LOX (促消退 — 脂氧素)        浅绿 #C5E6C5  (与 DHA-SPM 同色, 强调"促消退"角色)
"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

OUTPUT_PATH = (
    "results/figures/06_Pathway_Schematics/"
    "AA_LOX_pathway_schematic.png"
)


def draw_box(ax, x, y, width, height, text, color, fontsize=10):
    box = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width, height,
        boxstyle="round,pad=0.05,rounding_size=0.1",
        linewidth=1.2, edgecolor="#555555", facecolor=color, zorder=3,
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fontsize, color="#111111", zorder=4)


def draw_ortho_arrow(ax, x1, y1, x2, y2, y_turn):
    ax.plot([x1, x1], [y1, y_turn], color="#555555", lw=1.5, zorder=1)
    ax.plot([x1, x2], [y_turn, y_turn], color="#555555", lw=1.5, zorder=1)
    ax.annotate("", xy=(x2, y2), xytext=(x2, y_turn),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)


def draw_vert_arrow(ax, x, y1, y2):
    ax.annotate("", xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.5), zorder=2)


def main():
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(6, 7.7, "AA-LOX 通路：花生四烯酸经脂氧合酶生成的主要代谢物",
            ha="center", va="center", fontsize=15, fontweight="bold")

    draw_box(ax, 6, 6.8, 2.5, 0.7, "AA\n(花生四烯酸)", "#C9D7F8", 12)

    y_header = 5.4
    col_x = {"5-LOX": 2.5, "12-LOX": 6.0, "15-LOX": 9.5}
    headers = [
        ("5-LOX 酶促倾向\n(促炎 — 白三烯)",   "#F8C9A0", col_x["5-LOX"]),
        ("12-LOX 酶促倾向\n(混合 — HETE/hepoxilin)", "#F8E5A8", col_x["12-LOX"]),
        ("15-LOX 酶促倾向\n(促消退 — 脂氧素)", "#C5E6C5", col_x["15-LOX"]),
    ]
    for text, color, x in headers:
        draw_box(ax, x, y_header, 3.0, 0.85, text, color, 10.5)

    y_turn = 6.2
    for x in col_x.values():
        draw_ortho_arrow(ax, 6, 6.45, x, 5.85, y_turn)

    products = {
        "5-LOX":  [("5-HETE",    "#FBE2C8"),
                   ("LTA4",      "#FBE2C8"),
                   ("LTB4 / Cys-LTs", "#FBE2C8")],
        "12-LOX": [("12-HETE",     "#FCF2C8"),
                   ("Hepoxilin A3", "#FCF2C8"),
                   ("Trioxilin A3", "#FCF2C8")],
        "15-LOX": [("15-HETE",    "#E1F4E1"),
                   ("Lipoxin A4", "#E1F4E1"),
                   ("Lipoxin B4", "#E1F4E1")],
    }

    y_prod = [3.9, 2.9, 1.9]
    for branch, items in products.items():
        x = col_x[branch]
        prev_y = 4.95
        for (name, color), cy in zip(items, y_prod):
            draw_box(ax, x, cy, 2.0, 0.55, name, color, 10.5)
            draw_vert_arrow(ax, x, prev_y, cy + 0.30)
            prev_y = cy - 0.30

    ax.text(6, 0.8,
            "注：5-LOX 与 15-LOX 经胞外协同 (transcellular) 可生成脂氧素 (LXA4/LXB4)，作为内源性促消退介质",
            ha="center", va="center", fontsize=9, color="#666666", style="italic")

    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
