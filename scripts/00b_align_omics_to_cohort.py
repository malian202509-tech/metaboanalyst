"""00b. 组学样本对齐到主分析队列 (n=165)。

输入:
  - data/01_raw/ori.xlsx                                 198 样本的原始峰面积/浓度
  - data/01_raw/分娩特征分析_2019-2021_主分析队列n165.xlsx   主分析临床表

输出:
  - data/02_preprocessed/ori_n165_aligned.xlsx           对齐后的组学矩阵 (114 代谢物 × 165 样本 + 元信息列)
  - data/02_preprocessed/sample_alignment_n165.csv       样本 ID 映射审计表 (标本号 / omx_id / BMI_group)

ID 映射规则:
  19TPT0XXX  <->  A19-XXX
  20TPT0XXX  <->  B20-XXX
  21TPT0XXX  <->  C21-XXX

执行后, 下游 01_qc_and_preprocessing.py 直接读取 ori_n165_aligned.xlsx 即可。
"""
import re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OMX_IN = ROOT / 'data' / '01_raw' / 'ori.xlsx'
COH_IN = ROOT / 'data' / '01_raw' / '分娩特征分析_2019-2021_主分析队列n165.xlsx'
OMX_OUT = ROOT / 'data' / '02_preprocessed' / 'ori_n165_aligned.xlsx'
MAP_OUT = ROOT / 'data' / '02_preprocessed' / 'sample_alignment_n165.csv'

META_COLS = ['Metabolite Name','Chinese Name','Retention Time(min)','KEGG ID',
             'HMDB ID','Cellular Locations','Super Class','Class','Sub Class','Direct Parent']
UNIT_COL = '浓度单位'


def clin_to_omx(clin_id):
    """临床标本号 -> 组学样本ID (按年份前缀对应 A/B/C)."""
    m = re.match(r'^(\d{2})TPT0?(\d+)$', str(clin_id).strip())
    if not m:
        return None
    year, num = m.group(1), int(m.group(2))
    prefix = {'19': 'A19', '20': 'B20', '21': 'C21'}.get(year)
    if not prefix:
        return None
    return f'{prefix}-{num:03d}'


def main():
    print('=== 00b. 组学样本对齐到主分析队列 ===')

    # 读入
    omx = pd.read_excel(OMX_IN)
    coh = pd.read_excel(COH_IN)
    omx_samples = [c for c in omx.columns if c not in META_COLS and c != UNIT_COL]
    print(f'输入: 组学 {omx.shape[0]} 代谢物 × {len(omx_samples)} 样本')
    print(f'      临床主队列 n={len(coh)}')

    # 构建映射
    coh['omx_id'] = coh['标本号'].apply(clin_to_omx)
    if coh['omx_id'].isna().any():
        bad = coh[coh['omx_id'].isna()]
        raise ValueError(f'有 {len(bad)} 个临床标本号无法转换:\n{bad[["标本号"]].to_string()}')

    omx_set = set(omx_samples)
    coh['matched'] = coh['omx_id'].isin(omx_set)
    n_match = int(coh['matched'].sum())
    print(f'匹配: {n_match} / {len(coh)}')

    if n_match != len(coh):
        miss = coh.loc[~coh['matched'], ['标本号', 'omx_id', 'BMI_group']]
        raise ValueError(f'{len(miss)} 个临床样本在组学表中找不到:\n{miss.to_string()}')

    # 筛选: 保留 10 元信息列 + 165 样本列 + 单位列 (按原顺序)
    keep_samples = [s for s in omx_samples if s in set(coh['omx_id'])]
    keep_cols = META_COLS + keep_samples + [UNIT_COL]
    omx_out = omx[keep_cols].copy()
    print(f'输出: {omx_out.shape[0]} 代谢物 × {len(keep_samples)} 样本 (+ {len(META_COLS)} 元信息列 + 单位列)')

    # 写出
    OMX_OUT.parent.mkdir(parents=True, exist_ok=True)
    omx_out.to_excel(OMX_OUT, index=False)
    print(f'已写入: {OMX_OUT.relative_to(ROOT)}')

    # 审计映射表
    audit = coh[['标本号', 'omx_id', 'BMI_group']].copy()
    audit.to_csv(MAP_OUT, index=False, encoding='utf-8-sig')
    print(f'已写入: {MAP_OUT.relative_to(ROOT)}')

    # 报告
    print()
    print('--- 对齐后 BMI 分组分布 ---')
    print(audit['BMI_group'].value_counts().to_string())
    print()
    print('--- 按年份分布 (用于第 2 阶段批次诊断) ---')
    audit['year'] = audit['omx_id'].str[:3]
    print(pd.crosstab(audit['year'], audit['BMI_group'], margins=True).to_string())


if __name__ == '__main__':
    main()
