"""DHA 通路修订版 — 三条主通路 ratio 的森林图.

数据源: results/tables/lox_dha_and_autoox_corrected.csv (tier == main)
样式: 复用 15c_16hdohe_origin_check.py 的森林图风格
输出: results/figures/04_Metabolite_boxplots/dha_pathway_revised_forest.{png,svg}
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS"]
rcParams["axes.unicode_minus"] = False
rcParams["pdf.fonttype"] = 42
rcParams["svg.fonttype"] = "none"

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "results/tables/lox_dha_and_autoox_corrected.csv"
OUT_DIR = ROOT / "results/figures/04_Metabolite_boxplots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

KEEP = {
    "CYP4-ω 严格: 20-HDoHE / DHA":               ("CYP4-ω",   "#6FB56F"),
    "LOX SPM 前体: (14+17-HDoHE) / DHA":          ("LOX SPM",  "#D88080"),
    "自氧化主指数: (7+10+11+13+16-HDoHE) / DHA":  ("非酶",     "#8C8C8C"),
}

df = pd.read_csv(CSV)
df = df[df["ratio"].isin(KEEP.keys())].copy()
df["label"] = df["ratio"].map(lambda r: KEEP[r][0])
df["color"] = df["ratio"].map(lambda r: KEEP[r][1])
df = df.sort_values("log2FC").reset_index(drop=True)

fig, ax = plt.subplots(figsize=(9, 3.6))
ys = range(len(df))

for i, r in df.iterrows():
    ax.errorbar(r["log2FC"], i,
                xerr=[[r["log2FC"] - r["CI_lo"]], [r["CI_hi"] - r["log2FC"]]],
                fmt="o", color=r["color"], markersize=12, capsize=6, linewidth=2.4,
                markeredgecolor="white", markeredgewidth=1.0)
    ax.text(r["CI_hi"] + 0.03, i, f"p = {r['p_ols_hc3']:.4f}",
            va="center", ha="left", fontsize=11, color=r["color"], fontweight="bold")

ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
ax.set_yticks(list(ys))
ax.set_yticklabels(df["label"].tolist(), fontsize=12.5, fontweight="bold")
for lab, color in zip(ax.get_yticklabels(), df["color"]):
    lab.set_color(color)

ax.set_xlabel("log2FC", fontsize=12)
ax.grid(axis="x", alpha=0.25, linestyle="--")
ax.set_xlim(-0.25, df["CI_hi"].max() + 0.45)
ax.set_ylim(-0.6, len(df) - 0.4)
ax.tick_params(axis="y", length=0, pad=6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

plt.tight_layout()
for ext in ("png", "svg", "pdf"):
    out = OUT_DIR / f"dha_pathway_revised_forest.{ext}"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"  [OK] {out.relative_to(ROOT)}")
plt.close()
