"""11. 差异代谢物候选筛选 (基于 ANCOVA v2 输出).

筛选标准 (2026-05-18 用户锁定):
  p_limma < 0.05  AND  |log2FC| ≥ log2(1.2)

策略层级:
  - 标准 = SOP §4.3 'standard' tier 的宽松版 (用 raw p 替代 q_limma<0.10)
  - 不控制 FDR (BH-q 全不显著, n_feat 下 H0 假阳性期望约 3-4 个)
  - 论文 Methods 应明示: "hypothesis-generating candidates", 不能写 "differentially abundant metabolites"
  - 所有候选另由 Wilcoxon + 离群稳健性 (新版多列判据) 互校

为下游 (家族 heatmap / 酶活性比值 / 富集) 准备:
  - 候选表 (含 ANCOVA 全量列)
  - 候选 metabolite 列表 (供下游脚本 import)
  - Markdown 摘要 (论文 Results 草稿)

输入:
  results/tables/ancova_main_{80,50}.csv

输出:
  results/tables/diff_candidates_{80,50}.csv
  results/tables/diff_candidates_summary.md
"""
import os
import sys
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, OSError):
    pass

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / 'results' / 'tables'

P_THRESH = 0.05
FC_LOG2_THRESH = np.log2(1.2)  # ≈ 0.2630


def atomic_write_csv(df, target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=target.stem + '_', suffix='.csv.tmp', dir=str(target.parent))
    os.close(fd)
    try:
        df.to_csv(tmp, index=False, encoding='utf-8-sig')
        os.replace(tmp, target)
    except PermissionError as e:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise PermissionError(
            f'无法写入 {target.name}: 文件可能正被 Excel/WPS 打开. 请关闭后重跑.'
        ) from e
    except Exception:
        if os.path.exists(tmp):
            try: os.unlink(tmp)
            except OSError: pass
        raise


def screen_track(tag):
    src = TABLES / f'ancova_main_{tag}.csv'
    df = pd.read_csv(src, encoding='utf-8-sig')
    n_feat = len(df)

    mask = (df['p_limma'] < P_THRESH) & (df['log2FC'].abs() >= FC_LOG2_THRESH)
    cand = df.loc[mask].copy().sort_values('p_limma').reset_index(drop=True)

    out = TABLES / f'diff_candidates_{tag}.csv'
    atomic_write_csv(cand, out)

    n_h0_expect = n_feat * P_THRESH
    n_robust = int(cand['is_robust_to_outlier'].sum())
    n_up = int((cand['log2FC'] > 0).sum())
    n_down = int((cand['log2FC'] < 0).sum())

    print(f'\n=== 轨道 {tag}% (n_feat={n_feat}) ===')
    print(f'  筛选标准: p_limma<{P_THRESH} AND |log2FC|≥log2(1.2)={FC_LOG2_THRESH:.4f}')
    print(f'  候选数:           {len(cand)}')
    print(f'  H0 假阳性期望:    {n_h0_expect:.1f}')
    print(f'  离群稳健候选:     {n_robust} / {len(cand)}')
    print(f'  方向: ↑ {n_up}, ↓ {n_down}')
    print(f'  输出: {out.relative_to(ROOT)}')

    return cand, {
        'tag': tag, 'n_feat': n_feat, 'n_candidates': len(cand),
        'h0_expect': round(n_h0_expect, 2),
        'n_robust': n_robust, 'n_up': n_up, 'n_down': n_down,
    }


def make_summary_md(cand80, cand50, audit80, audit50):
    """生成 Markdown 摘要 (论文 Results 草稿)."""
    lines = []
    lines.append('# 差异代谢物候选筛选结果\n')
    lines.append(f'**生成时间**: 2026-05-18  ')
    lines.append(f'**筛选标准**: `p_limma < 0.05` AND `|log2FC| ≥ log2(1.2) ≈ 0.263`  ')
    lines.append(f'**输入**: `results/tables/ancova_main_{{80,50}}.csv` (ANCOVA v2, 方案 D + 3 协变量)\n')

    lines.append('## 双轨结果对比\n')
    lines.append('| 项 | 80% 主轨 | 50% 探索轨 |')
    lines.append('|---|---|---|')
    lines.append(f'| 总特征数 | {audit80["n_feat"]} | {audit50["n_feat"]} |')
    lines.append(f'| 候选数 | **{audit80["n_candidates"]}** | **{audit50["n_candidates"]}** |')
    lines.append(f'| H0 假阳性期望 | {audit80["h0_expect"]} | {audit50["h0_expect"]} |')
    lines.append(f'| 实际命中 vs 期望 | {audit80["n_candidates"]} vs {audit80["h0_expect"]} '
                 f'({"低于期望 ✓" if audit80["n_candidates"] < audit80["h0_expect"] else "略高 ⚠"}) | '
                 f'{audit50["n_candidates"]} vs {audit50["h0_expect"]} '
                 f'({"低于期望 ✓" if audit50["n_candidates"] < audit50["h0_expect"] else "略高 ⚠"}) |')
    lines.append(f'| 离群稳健候选 | {audit80["n_robust"]} / {audit80["n_candidates"]} | '
                 f'{audit50["n_robust"]} / {audit50["n_candidates"]} |')
    lines.append(f'| 方向 (↑ / ↓) | {audit80["n_up"]} / {audit80["n_down"]} | '
                 f'{audit50["n_up"]} / {audit50["n_down"]} |\n')

    def render_table(cand, has_also_col=False):
        cols = ['Metabolite Name', 'Chinese Name', 'Class', 'direction',
                'log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi',
                'p_ols_hc3', 'p_limma', 'p_wilcoxon',
                'n_outliers', 'attenuation_pct', 'is_robust_to_outlier']
        if has_also_col and 'also_in_80_main' in cand.columns:
            cols.append('also_in_80_main')
        sub = cand[cols].copy()
        for c in ['log2FC', 'log2FC_CI_lo', 'log2FC_CI_hi', 'attenuation_pct']:
            sub[c] = sub[c].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '')
        for c in ['p_ols_hc3', 'p_limma', 'p_wilcoxon']:
            sub[c] = sub[c].apply(lambda x: f'{x:.4f}' if pd.notna(x) else '')
        header = '| ' + ' | '.join(cols) + ' |'
        sep = '| ' + ' | '.join(['---'] * len(cols)) + ' |'
        body = ['| ' + ' | '.join(str(v) for v in row) + ' |' for row in sub.values]
        return '\n'.join([header, sep, *body])

    lines.append('## 80% 主轨候选清单 (主分析)\n')
    if len(cand80) > 0:
        lines.append(render_table(cand80) + '\n')
    else:
        lines.append('_无候选_\n')

    lines.append('## 50% 探索轨候选清单 (补救 / 敏感性)\n')
    if len(cand50) > 0:
        lines.append(render_table(cand50, has_also_col=True) + '\n')
    else:
        lines.append('_无候选_\n')

    overlap = set(cand80['Metabolite Name']) & set(cand50['Metabolite Name'])
    only50 = set(cand50['Metabolite Name']) - set(cand80['Metabolite Name'])
    lines.append('## 双轨重叠分析\n')
    lines.append(f'- 80 主轨 ∩ 50 探索轨 = **{len(overlap)}** 个 (核心一致信号)')
    lines.append(f'- 仅 50 轨独有 = **{len(only50)}** 个 (50% 检出阈值放宽抢救的低丰度信号)')
    if overlap:
        lines.append(f'\n**核心一致候选 ({len(overlap)})**:')
        for m in sorted(overlap):
            row = cand80[cand80['Metabolite Name'] == m].iloc[0]
            lines.append(f'  - `{m}` → {row["direction"]} log2FC={row["log2FC"]:.3f}, '
                         f'p_limma={row["p_limma"]:.4f}')
    if only50:
        lines.append(f'\n**50 轨独有 ({len(only50)})**:')
        for m in sorted(only50):
            row = cand50[cand50['Metabolite Name'] == m].iloc[0]
            lines.append(f'  - `{m}` → {row["direction"]} log2FC={row["log2FC"]:.3f}, '
                         f'p_limma={row["p_limma"]:.4f}')

    lines.append('\n## Methods 写作要点\n')
    lines.append('- 筛选标准应明示为 **hypothesis-generating tier**, 不作 FDR 控制声明')
    lines.append('- 推荐措辞: "Candidate metabolites were defined by `p_limma < 0.05` and '
                 '`|log2FC| ≥ log2(1.2)`. As no metabolite reached BH-q < 0.05, this threshold '
                 'serves as a hypothesis-generating tier rather than confirmed differential abundance."')
    lines.append('- 每个候选另由 OLS+HC3、Mann-Whitney U、studentized-residual outlier 稳健性互校 (见表)')
    lines.append('- 与上游 PERMANOVA (BMI p≈0.40) 与 OPLS-DA (Q² 为负) 的全局阴性一致, 进一步支持 '
                 '"per-feature evidence is weak but family-level direction is informative"')

    return '\n'.join(lines)


def main():
    print('=== 11. 差异代谢物候选筛选 ===')
    print(f'标准: p_limma<{P_THRESH} AND |log2FC|≥log2(1.2)={FC_LOG2_THRESH:.4f}')

    cand80, audit80 = screen_track('80')
    cand50, audit50 = screen_track('50')

    md = make_summary_md(cand80, cand50, audit80, audit50)
    md_path = TABLES / 'diff_candidates_summary.md'
    md_path.write_text(md, encoding='utf-8')
    print(f'\n摘要: {md_path.relative_to(ROOT)}')

    overlap = set(cand80['Metabolite Name']) & set(cand50['Metabolite Name'])
    only50 = set(cand50['Metabolite Name']) - set(cand80['Metabolite Name'])
    print(f'\n=== 双轨重叠 ===')
    print(f'  80 ∩ 50 = {len(overlap)} (核心一致)')
    print(f'  50 独有 = {len(only50)} (低丰度抢救)')


if __name__ == '__main__':
    main()
