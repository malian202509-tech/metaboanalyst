"""Redraw participant flow diagram as clean vector graphics (SVG + PDF), Times New Roman."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FIG_W, FIG_H = 10.5, 8.5
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_axis_off()

MAIN_X, MAIN_W = 12, 38
EXCL_X, EXCL_W = 54, 44
LW = 1.1
FS_MAIN = 13
FS_EXCL = 12
FS_SUB = 11


def draw_box(x, y, w, h, lines, fontsize, align="center"):
    ax.add_patch(Rectangle((x, y), w, h, fill=False, edgecolor="black", linewidth=LW))
    if align == "center":
        ax.text(x + w / 2, y + h / 2, "\n".join(lines),
                ha="center", va="center", fontsize=fontsize)
    else:
        ax.text(x + 1.5, y + h - 1.4, "\n".join(lines),
                ha="left", va="top", fontsize=fontsize)


def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                 arrowstyle="-|>", mutation_scale=14,
                                 color="black", linewidth=LW,
                                 shrinkA=0, shrinkB=0))


main_boxes = [
    {"y": 88, "h": 7, "text": ["172,385 enrolled"]},
    {"y": 64, "h": 7, "text": ["151,359 children aged 3 to 5 years"]},
    {"y": 40, "h": 7, "text": ["151,267 singleton live births"]},
    {"y": 6,  "h": 7, "text": ["137,269 mother-child pairs"]},
]
for b in main_boxes:
    draw_box(MAIN_X, b["y"], MAIN_W, b["h"], b["text"], FS_MAIN)

excl_boxes = [
    {"y": 78, "h": 6,  "text": ["21,026 excluded (child age < 3 or ≥ 6 years)"], "align": "center", "fs": FS_EXCL},
    {"y": 54, "h": 6,  "text": ["92 excluded (multiple gestations)"], "align": "center", "fs": FS_EXCL},
    {"y": 18, "h": 18, "text": [
        "13,998 excluded",
        "   11,536 duplicate questionnaire submissions",
        "    1,665 missing key covariates",
        "      763 missing maternal pre-pregnancy height/weight data",
        "       34 missing LDCDQ scores",
    ], "align": "left", "fs": FS_SUB},
]
for b in excl_boxes:
    draw_box(EXCL_X, b["y"], EXCL_W, b["h"], b["text"], b["fs"], align=b["align"])

mid_x = MAIN_X + MAIN_W / 2
for top, bot in [(88, 71), (64, 47), (40, 13)]:
    arrow(mid_x, top, mid_x, bot + 0.05)

branch_ys = [81, 57, 27]
for y in branch_ys:
    ax.plot([mid_x, EXCL_X - 0.05], [y, y], color="black", linewidth=LW)
    arrow(EXCL_X - 0.05, y, EXCL_X - 0.05 + 0.5, y)

for stem in [(88, 81), (64, 57), (40, 27)]:
    pass

plt.tight_layout()

svg_path = OUT_DIR / "flowchart_consort.svg"
pdf_path = OUT_DIR / "flowchart_consort.pdf"
png_path = OUT_DIR / "flowchart_consort.png"
fig.savefig(svg_path, format="svg", bbox_inches="tight")
fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
fig.savefig(png_path, format="png", dpi=400, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {svg_path}")
print(f"Saved: {pdf_path}")
print(f"Saved: {png_path}")
