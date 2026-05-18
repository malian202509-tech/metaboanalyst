"""13. 家族 heatmap (SOP §5.1): 按 7 大家族分块展示氧化脂质 z-score 矩阵.

设计:
  - 行: 73 个特征 (50 轨) 或 67 个 (80 轨), 按 family_main → 候选优先 → t_limma 排序
  - 列: 165 样本, 按 BMI 分组 (正常 n=120 | 超重肥胖 n=45) 排列; 组内按样本均值
  - 值: log2 矩阵每行 z-score (跨样本); 限幅 [-3, 3] 防离群拉色阶
  - 行注释 (左侧色条):
      · family_main (7 色块)
      · is_candidate ★ 标记
      · direction ↑↓ (来自 ANCOVA log2FC 符号)
  - 列注释 (顶部色条): BMI_group (2 色)

输入:
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat.xlsx
  data/02_preprocessed/sample_alignment_n165.csv
  data/02_preprocessed/metabolite_family_map.csv
  results/tables/ancova_main_{80,50}.csv (用 t_limma + log2FC + tier)

输出:
  results/figures/03_Heatmaps/family_heatmap_{80main,50expl}.png (300 DPI)
  results/figures/03_Heatmaps/family_heatmap_audit.csv (每特征绘图属性审计)
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec

# Windows 中文字体配置
matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

ROOT = Path(__file__).resolve().parent.parent
COMBAT_DIR = ROOT / 'data' / '03_batch_corrected'
PREP_DIR = ROOT / 'data' / '02_preprocessed'
TABLES = ROOT / 'results' / 'tables'
FIG_DIR = ROOT / 'results' / 'figures' / '03_Heatmaps'
FIG_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'

# 7 家族固定显示顺序 + 配色 (论文 friendly, 区分度高)
FAMILY_ORDER = [
    'Free PUFA',
    'Endocannabinoid',
    'LA-oxylipin',
    'AA-COX',
    'AA-LOX',
    'AA-CYP/sEH',
    'EPA/DHA/DPA-oxylipin',
]
FAMILY_COLORS = {
    'Free PUFA':            '#A0A0A0',
    'Endocannabinoid':      '#9E9AC8',
    'LA-oxylipin':          '#74C476',
    'AA-COX':               '#FB6A4A',
    'AA-LOX':               '#FD8D3C',
    'AA-CYP/sEH':           '#6BAED6',
    'EPA/DHA/DPA-oxylipin': '#41AB5D',
}
BMI_COLORS = {'正常': '#67A9CF', '超重肥胖': '#EF8A62'}

Z_CLIP = 3.0


def load_data(tag):
    combat_xlsx = COMBAT_DIR / f'ori_n165_filtered{tag}_log2_combat.xlsx'
    df = pd.read_excel(combat_xlsx)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]

    align = pd.read_csv(PREP_DIR / 'sample_alignment_n165.csv')
    s2g = dict(zip(align['omx_id'], align['BMI_group']))

    fmap = pd.read_csv(PREP_DIR / 'metabolite_family_map.csv', encoding='utf-8-sig')
    fmap_idx = fmap.set_index('Metabolite Name')

    anc = pd.read_csv(TABLES / f'ancova_main_{tag}.csv', encoding='utf-8-sig')
    anc_idx = anc.set_index('Metabolite Name')

    return df, sample_cols, s2g, fmap_idx, anc_idx


def build_matrix(df, sample_cols, s2g, fmap_idx, anc_idx):
    """组装绘图矩阵 (行 z-scored), 行按 family → 候选优先 → t_limma 排序."""
    metabs = df['Metabolite Name'].tolist()

    # 每行 metadata
    rows_meta = []
    for m in metabs:
        f_row = fmap_idx.loc[m]
        a_row = anc_idx.loc[m]
        rows_meta.append({
            'name': m,
            'short': f_row['short_label'],
            'family_main': f_row['family_main'],
            'family_sub': f_row['family_sub'],
            'substrate': f_row['substrate'],
            'enzyme': f_row['enzyme'],
            'is_cand_80': bool(f_row.get('is_candidate_80', False)),
            'is_cand_50': bool(f_row.get('is_candidate_50', False)),
            't_limma': float(a_row['t_limma']),
            'log2FC': float(a_row['log2FC']),
            'p_limma': float(a_row['p_limma']),
            'tier': str(a_row['tier']),
        })
    meta_df = pd.DataFrame(rows_meta)
    # 排序: family_main (按 FAMILY_ORDER) → 候选优先 → |t_limma| 降序
    meta_df['fam_rank'] = meta_df['family_main'].map({f: i for i, f in enumerate(FAMILY_ORDER)})
    meta_df['cand_rank'] = ~(meta_df['is_cand_80'] | meta_df['is_cand_50'])  # 候选在前
    meta_df['abs_t'] = -meta_df['t_limma'].abs()  # 大 |t| 在前 (descending)
    meta_df = meta_df.sort_values(['fam_rank', 'cand_rank', 'abs_t']).reset_index(drop=True)

    # 取数 + z-score
    matrix = df.set_index('Metabolite Name').loc[meta_df['name'], sample_cols]
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1), axis=0)
    z = z.clip(-Z_CLIP, Z_CLIP)

    # 列按 BMI_group 排列
    col_meta = pd.DataFrame({'omx': sample_cols,
                             'bmi': [s2g.get(c, '?') for c in sample_cols]})
    col_meta['bmi_rank'] = col_meta['bmi'].map({'正常': 0, '超重肥胖': 1, '?': 2})
    # 组内按列均值排 (从高到低)
    col_mean = matrix.mean(axis=0)
    col_meta['col_mean'] = col_meta['omx'].map(col_mean.to_dict())
    col_meta = col_meta.sort_values(['bmi_rank', 'col_mean'], ascending=[True, False]).reset_index(drop=True)
    z = z[col_meta['omx'].tolist()]

    return z, meta_df, col_meta


def draw_heatmap(z, meta_df, col_meta, tag, title, out_path):
    n_feat, n_samp = z.shape
    cand_col = 'is_cand_80' if tag == '80' else 'is_cand_50'
    is_cand = meta_df[cand_col].values
    fam_arr = meta_df['family_main'].values
    direction_arr = np.where(meta_df['log2FC'] > 0, '↑', '↓')

    # 布局: 三个左色条(family/cand/direction) + 主热图 + colorbar; 顶部 BMI 色条
    fig_h = max(8.0, n_feat * 0.16 + 1.5)
    fig_w = 14.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        2, 5, figure=fig,
        height_ratios=[0.25, fig_h - 0.5],
        width_ratios=[0.25, 0.18, 0.18, 8.5, 0.30],
        wspace=0.04, hspace=0.02,
    )

    # --- 顶部 BMI 色条 ---
    ax_top = fig.add_subplot(gs[0, 3])
    bmi_codes = col_meta['bmi'].map({'正常': 0, '超重肥胖': 1}).values
    cmap_bmi = ListedColormap([BMI_COLORS['正常'], BMI_COLORS['超重肥胖']])
    ax_top.imshow(bmi_codes[None, :], aspect='auto', cmap=cmap_bmi, interpolation='nearest')
    ax_top.set_xticks([]); ax_top.set_yticks([])
    # 在两组中点写 label
    n_norm = int((col_meta['bmi'] == '正常').sum())
    ax_top.text(n_norm / 2, 0, f'正常 (n={n_norm})', ha='center', va='center',
                fontsize=11, color='white', fontweight='bold')
    ax_top.text(n_norm + (n_samp - n_norm) / 2, 0, f'超重肥胖 (n={n_samp - n_norm})',
                ha='center', va='center', fontsize=11, color='white', fontweight='bold')
    ax_top.set_title(title, fontsize=13, fontweight='bold', pad=10)

    # --- 左侧 family 色条 ---
    ax_fam = fig.add_subplot(gs[1, 0])
    fam_to_idx = {f: i for i, f in enumerate(FAMILY_ORDER)}
    fam_codes = np.array([fam_to_idx[f] for f in fam_arr])
    cmap_fam = ListedColormap([FAMILY_COLORS[f] for f in FAMILY_ORDER])
    ax_fam.imshow(fam_codes[:, None], aspect='auto', cmap=cmap_fam, interpolation='nearest')
    ax_fam.set_xticks([]); ax_fam.set_yticks([])
    # 在每家族段中心写名字
    for f in FAMILY_ORDER:
        positions = np.where(fam_arr == f)[0]
        if len(positions) == 0: continue
        ctr = positions.mean()
        ax_fam.text(0, ctr, f, ha='center', va='center',
                    rotation=90, fontsize=8.5, color='black', fontweight='bold')

    # --- 候选标记 (★ 用纯色块) ---
    ax_cand = fig.add_subplot(gs[1, 1])
    cand_codes = is_cand.astype(int)
    cmap_cand = ListedColormap(['#F5F5F5', '#1A1A1A'])
    ax_cand.imshow(cand_codes[:, None], aspect='auto', cmap=cmap_cand, interpolation='nearest')
    ax_cand.set_xticks([]); ax_cand.set_yticks([])
    for i, c in enumerate(is_cand):
        if c:
            ax_cand.text(0, i, '★', ha='center', va='center',
                         fontsize=10, color='gold', fontweight='bold')

    # --- 方向 ↑↓ ---
    ax_dir = fig.add_subplot(gs[1, 2])
    dir_codes = (meta_df['log2FC'].values > 0).astype(int)  # 0=down,1=up
    cmap_dir = ListedColormap(['#3182BD', '#DE2D26'])  # blue ↓ / red ↑
    ax_dir.imshow(dir_codes[:, None], aspect='auto', cmap=cmap_dir, interpolation='nearest')
    ax_dir.set_xticks([]); ax_dir.set_yticks([])
    for i, d in enumerate(direction_arr):
        ax_dir.text(0, i, d, ha='center', va='center',
                    fontsize=8, color='white', fontweight='bold')

    # --- 主 heatmap ---
    ax = fig.add_subplot(gs[1, 3])
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-Z_CLIP, vmax=Z_CLIP,
                   interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks(range(n_feat))
    # 行标签 (short_label, 候选加 ★ 前缀)
    ylabels = []
    for i, r in meta_df.iterrows():
        prefix = '★ ' if (r[cand_col]) else '  '
        ylabels.append(f'{prefix}{r["short"]}')
    ax.set_yticklabels(ylabels, fontsize=7.5)
    # 候选 label 加粗加色
    for tick, c in zip(ax.get_yticklabels(), is_cand):
        if c:
            tick.set_fontweight('bold')
            tick.set_color('#D62728')

    # 家族分界线 (横线)
    fam_change = np.where(np.array([fam_arr[i] != fam_arr[i-1] for i in range(1, n_feat)]))[0]
    for boundary in fam_change:
        ax.axhline(boundary + 0.5, color='black', lw=0.8, alpha=0.6)
    # BMI 组分界线 (竖线)
    ax.axvline(n_norm - 0.5, color='black', lw=1.2, alpha=0.8)

    # --- colorbar ---
    ax_cb = fig.add_subplot(gs[1, 4])
    cb = plt.colorbar(im, cax=ax_cb)
    cb.set_label('z-score (per-feature)', fontsize=10)

    # 图例 (家族 + 候选 + 方向) — 放右下角
    legend_handles = [mpatches.Patch(color=FAMILY_COLORS[f], label=f) for f in FAMILY_ORDER]
    legend_handles += [
        mpatches.Patch(color='#DE2D26', label='↑ 超重肥胖 > 正常'),
        mpatches.Patch(color='#3182BD', label='↓ 超重肥胖 < 正常'),
        mpatches.Patch(color='gold', label='★ 差异候选 (p_limma<0.05, |FC|≥1.2, robust)'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=5, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.005), frameon=False)

    plt.subplots_adjust(left=0.04, right=0.96, top=0.96, bottom=0.06)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  ✓ {out_path.relative_to(ROOT)}')


def main():
    print('=== 13. 家族 heatmap ===\n')
    audit_rows = []

    for tag, label in [('80', '80% 主分析轨 (67 特征)'),
                       ('50', '50% 探索轨 (73 特征)')]:
        print(f'轨道 {tag}: {label}')
        df, sample_cols, s2g, fmap_idx, anc_idx = load_data(tag)
        z, meta_df, col_meta = build_matrix(df, sample_cols, s2g, fmap_idx, anc_idx)
        out_path = FIG_DIR / f'family_heatmap_{"80main" if tag == "80" else "50expl"}.png'
        title = f'Family-grouped oxylipin heatmap — {label}\n(z-score per feature, ComBat-by-batch corrected)'
        draw_heatmap(z, meta_df, col_meta, tag, title, out_path)

        # 审计行
        for _, r in meta_df.iterrows():
            audit_rows.append({
                'track': tag,
                'family_main': r['family_main'],
                'short': r['short'],
                'name': r['name'],
                'is_candidate': bool(r['is_cand_80' if tag == '80' else 'is_cand_50']),
                'log2FC': round(r['log2FC'], 4),
                't_limma': round(r['t_limma'], 4),
                'p_limma': round(r['p_limma'], 6),
            })

    audit = pd.DataFrame(audit_rows)
    audit_path = FIG_DIR / 'family_heatmap_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计: {audit_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
