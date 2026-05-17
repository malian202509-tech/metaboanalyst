"""生成主分析队列 (n=165) 的 Table 1 基线表。
- 连续变量自动检测正态性 (Shapiro-Wilk on each group), 选 t / Welch / Mann-Whitney
- 分类变量按期望频数选 χ² 或 Fisher 精确检验
- SMD 按四档阈值标注: <0.10 / 0.10-0.25 / 0.25-0.50 / >0.50
- 输出 Excel + Markdown 两种格式
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import random_table


def fisher_freeman_halton(ct, n_sim=20000, random_state=42):
    """RxC Fisher 精确检验的 Monte Carlo 实现 (Freeman-Halton 扩展).
    固定行列边际 -> 模拟 n_sim 张表 -> p = Pr(LR统计量 >= 观测).
    用对数似然比 (G²) 作为统计量, 对零格更稳定."""
    obs = np.asarray(ct, dtype=float)
    row_sums = obs.sum(axis=1).astype(int)
    col_sums = obs.sum(axis=0).astype(int)

    def g2(table):
        n = table.sum()
        if n == 0: return 0.0
        rs = table.sum(axis=1, keepdims=True)
        cs = table.sum(axis=0, keepdims=True)
        exp = rs * cs / n
        mask = (table > 0) & (exp > 0)
        return float(2 * np.sum(table[mask] * np.log(table[mask] / exp[mask])))

    g2_obs = g2(obs)
    rt = random_table(row_sums, col_sums, seed=random_state)
    sim = rt.rvs(size=n_sim)
    g2_sim = np.array([g2(s) for s in sim])
    return (np.sum(g2_sim >= g2_obs - 1e-12) + 1) / (n_sim + 1)  # 加1平滑

SRC = '分娩特征分析_2019-2021_主分析队列n165.xlsx'
GROUP_COL = 'BMI_group'
GROUP_REF = '正常'
GROUP_EXP = '超重肥胖'

CONT_VARS = [
    ('孕妇年龄', '孕妇年龄 (岁)'),
    ('GA_decimal', '孕周 (周)'),
    ('孕期增重', '孕期增重 (kg)'),
    ('出生体重', '出生体重 (g)'),
]

CAT_VARS = [
    ('胎儿性别', '胎儿性别 [n (%)]', None),
    ('GA_category', '孕周分类 (ACOG) [n (%)]', ['早期足月','足月','晚期/过期足月']),
    ('出生体重_分类', '出生体重分类 [n (%)]', ['低出生体重','正常体重','巨大儿']),
    ('GDM', 'GDM [n (%)]', [0,1]),
    ('PIH（妊娠期高血压）', 'PIH [n (%)]', [0,1]),
    ('脂代谢异常', '脂代谢异常 [n (%)]', [0,1]),
]


def smd_continuous(x1, x2):
    x1, x2 = np.asarray(x1, float), np.asarray(x2, float)
    x1, x2 = x1[~np.isnan(x1)], x2[~np.isnan(x2)]
    m1, m2 = x1.mean(), x2.mean()
    v1, v2 = x1.var(ddof=1), x2.var(ddof=1)
    denom = np.sqrt((v1 + v2) / 2)
    return abs(m1 - m2) / denom if denom > 0 else np.nan


def smd_categorical(s1, s2, levels):
    """多分类 SMD (Haukoos & Lewis 2015): sqrt(T' S^-1 T), T 与 S 用 K-1 维。"""
    s1, s2 = s1.dropna(), s2.dropna()
    K = len(levels)
    if K < 2:
        return np.nan
    p1 = np.array([(s1 == lv).mean() for lv in levels[:-1]])
    p2 = np.array([(s2 == lv).mean() for lv in levels[:-1]])
    T = p1 - p2

    def cov(p):
        S = -np.outer(p, p)
        np.fill_diagonal(S, p * (1 - p))
        return S
    S = (cov(p1) + cov(p2)) / 2
    try:
        return float(np.sqrt(T @ np.linalg.pinv(S) @ T))
    except np.linalg.LinAlgError:
        return np.nan


def smd_flag(s):
    if pd.isna(s): return ''
    if s < 0.10: return '√'
    if s < 0.25: return '⚠'
    if s < 0.50: return '!'
    return '!!'


def fmt_cont(x):
    """根据正态性返回 (描述字符串, P值, 用的检验名)。"""
    x = x.dropna()
    return x


def describe_cont(s, normal):
    s = s.dropna()
    if normal:
        return f'{s.mean():.2f} ± {s.std(ddof=1):.2f}'
    q1, q2, q3 = s.quantile([0.25, 0.5, 0.75])
    return f'{q2:.2f} ({q1:.2f}, {q3:.2f})'


def main():
    df = pd.read_excel(SRC)
    g1 = df[df[GROUP_COL] == GROUP_REF]
    g2 = df[df[GROUP_COL] == GROUP_EXP]
    n1, n2, n_all = len(g1), len(g2), len(df)
    print(f'读入: n={n_all}  正常={n1}  超重肥胖={n2}')

    rows = []

    # --- 连续变量 ---
    for col, label in CONT_VARS:
        x1, x2 = g1[col].dropna(), g2[col].dropna()
        x_all = df[col].dropna()
        # 正态性: 两组分别 Shapiro
        p_sw1 = stats.shapiro(x1).pvalue if len(x1) >= 3 else 0
        p_sw2 = stats.shapiro(x2).pvalue if len(x2) >= 3 else 0
        normal = (p_sw1 > 0.05) and (p_sw2 > 0.05)
        if normal:
            # Levene 检验方差齐
            p_lev = stats.levene(x1, x2).pvalue
            if p_lev > 0.05:
                stat, pval = stats.ttest_ind(x1, x2, equal_var=True)
                test = 't'
            else:
                stat, pval = stats.ttest_ind(x1, x2, equal_var=False)
                test = "Welch's t"
        else:
            stat, pval = stats.mannwhitneyu(x1, x2, alternative='two-sided')
            test = 'Mann-Whitney U'
        smd = smd_continuous(x1, x2)
        miss = df[col].isna().sum()
        rows.append({
            '变量': label + (f'  [缺失 {miss}]' if miss > 0 else ''),
            f'正常 (n={n1})': describe_cont(x1, normal),
            f'超重肥胖 (n={n2})': describe_cont(x2, normal),
            f'Overall (n={n_all})': describe_cont(x_all, normal),
            'P': ('<0.001' if pval < 0.001 else f'{pval:.3f}'),
            '|SMD|': f'{smd:.3f}',
            'SMD 提示': smd_flag(smd),
            '检验': test,
        })

    # --- 分类变量 ---
    for col, label, levels in CAT_VARS:
        s1, s2, s_all = g1[col].dropna(), g2[col].dropna(), df[col].dropna()
        if levels is None:
            levels = sorted(s_all.unique().tolist())
        # 列联表
        ct = pd.crosstab(df[col], df[GROUP_COL]).reindex(index=levels, columns=[GROUP_REF, GROUP_EXP], fill_value=0)
        # 检验: 期望频数 <5 用 Fisher
        try:
            chi2, pval, dof, expected = stats.chi2_contingency(ct)
            use_fisher = (expected < 5).any()
        except ValueError:
            use_fisher = True
            pval = np.nan
        if use_fisher:
            if ct.shape == (2, 2):
                _, pval = stats.fisher_exact(ct)
                test = "Fisher 精确"
            else:
                pval = fisher_freeman_halton(ct, n_sim=20000, random_state=42)
                test = "Fisher-Freeman-Halton (20k模拟)"
        else:
            test = 'χ²'

        smd = smd_categorical(s1, s2, levels)

        # 父行
        miss = df[col].isna().sum()
        rows.append({
            '变量': label + (f'  [缺失 {miss}]' if miss > 0 else ''),
            f'正常 (n={n1})': '',
            f'超重肥胖 (n={n2})': '',
            f'Overall (n={n_all})': '',
            'P': ('<0.001' if (not pd.isna(pval) and pval < 0.001) else (f'{pval:.3f}' if not pd.isna(pval) else 'NA')),
            '|SMD|': f'{smd:.3f}' if not pd.isna(smd) else 'NA',
            'SMD 提示': smd_flag(smd),
            '检验': test,
        })
        # 子行 (各 level)
        for lv in levels:
            n1l = int((s1 == lv).sum()); p1l = n1l / max(len(s1), 1) * 100
            n2l = int((s2 == lv).sum()); p2l = n2l / max(len(s2), 1) * 100
            nal = int((s_all == lv).sum()); pal = nal / max(len(s_all), 1) * 100
            disp_lv = '是' if (col in ['GDM','PIH（妊娠期高血压）','脂代谢异常'] and lv == 1) else ('否' if (col in ['GDM','PIH（妊娠期高血压）','脂代谢异常'] and lv == 0) else str(lv))
            rows.append({
                '变量': f'    {disp_lv}',
                f'正常 (n={n1})': f'{n1l} ({p1l:.1f})',
                f'超重肥胖 (n={n2})': f'{n2l} ({p2l:.1f})',
                f'Overall (n={n_all})': f'{nal} ({pal:.1f})',
                'P': '',
                '|SMD|': '',
                'SMD 提示': '',
                '检验': '',
            })

    out = pd.DataFrame(rows)
    out.to_excel('Table1_baseline_n165.xlsx', index=False)
    print('已写入 Table1_baseline_n165.xlsx')

    # Markdown
    md_lines = ['| ' + ' | '.join(out.columns) + ' |']
    md_lines.append('| ' + ' | '.join(['---']*len(out.columns)) + ' |')
    for _, r in out.iterrows():
        md_lines.append('| ' + ' | '.join(str(v) for v in r.values) + ' |')
    md = '\n'.join(md_lines)
    with open('Table1_baseline_n165.md', 'w', encoding='utf-8') as f:
        f.write('# Table 1. Baseline characteristics (n=165)\n\n')
        f.write(md + '\n\n')
        f.write('**SMD 提示**: √ = 可忽略 (<0.10);  ⚠ = 轻度不均衡 (0.10–0.25, 建议入 GLM);  ! = 明显不均衡 (0.25–0.50);  !! = 严重不均衡 (>0.50)\n')
    print('已写入 Table1_baseline_n165.md')

    print()
    print(md)


if __name__ == '__main__':
    main()
