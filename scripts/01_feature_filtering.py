"""01. 特征过滤 (按分组检出率).

依据:
  - Bijlsma et al. Anal Chem 2006: modified 80% rule (>=80% in at least one group)
  - Do et al. Metabolomics 2018: 50% missing 是插补可靠性上限

策略 (双轨并行):
  - 主分析: 80% 规则 -> 严格, 残余缺失率低, 插补稳健
  - 探索性: 50% 规则 -> 多保留 6 个 LOD 边缘但生物学关键的分子
                       (LTB4, LTE4, 8-HDHA, 11-HEPE 等)

输入:
  data/02_preprocessed/ori_n165_aligned.xlsx       (114 代谢物 x 165 样本)
  data/02_preprocessed/sample_alignment_n165.csv   (含 BMI_group)

输出:
  data/02_preprocessed/ori_n165_filtered80.xlsx    主分析输入 (67 代谢物)
  data/02_preprocessed/ori_n165_filtered50.xlsx    探索性输入 (73 代谢物)
  data/02_preprocessed/feature_filtering_audit.csv 每个代谢物的检出率与去留决定

缺失定义: NaN 或 0 (LC-MS 中 0 通常为 LOD 截断, 本数据仅含 NaN, 见诊断脚本)
"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / 'data' / '02_preprocessed' / 'ori_n165_aligned.xlsx'
MAP = ROOT / 'data' / '02_preprocessed' / 'sample_alignment_n165.csv'
OUT_DIR = ROOT / 'data' / '02_preprocessed'

OUT_MAIN = OUT_DIR / 'ori_n165_filtered80.xlsx'
OUT_EXPL = OUT_DIR / 'ori_n165_filtered50.xlsx'
OUT_AUDIT = OUT_DIR / 'feature_filtering_audit.csv'

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'
THR_MAIN = 0.80
THR_EXPL = 0.50


def main():
    print('=== 01. 特征过滤 (双轨: 80% 主 + 50% 探索) ===')
    df = pd.read_excel(SRC)
    mp = pd.read_csv(MAP)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    print(f'输入: {df.shape[0]} 代谢物 x {len(sample_cols)} 样本')

    # 缺失矩阵: NaN 或 0
    sub = df[sample_cols]
    is_missing = sub.isna() | (sub == 0)

    # 分组样本列
    g1_cols = [c for c in mp.loc[mp['BMI_group'] == '正常', 'omx_id'] if c in sample_cols]
    g2_cols = [c for c in mp.loc[mp['BMI_group'] == '超重肥胖', 'omx_id'] if c in sample_cols]

    # 检出率
    det_g1 = 1 - is_missing[g1_cols].mean(axis=1)
    det_g2 = 1 - is_missing[g2_cols].mean(axis=1)
    max_det = np.maximum(det_g1, det_g2)

    # 两档筛选 (max group detection rate >= threshold)
    keep_main = max_det >= THR_MAIN
    keep_expl = max_det >= THR_EXPL

    print(f'\n80% 规则: 保留 {keep_main.sum()} / {len(df)} 代谢物')
    print(f'50% 规则: 保留 {keep_expl.sum()} / {len(df)} 代谢物')
    print(f'仅 50% 多保留 (探索性专属): {int((keep_expl & ~keep_main).sum())} 个')

    # 落盘 (保留全部元信息列 + 样本列 + 单位列, 行用布尔过滤)
    keep_cols = META_COLS + sample_cols + [UNIT_COL]
    df_main = df.loc[keep_main, keep_cols].reset_index(drop=True)
    df_expl = df.loc[keep_expl, keep_cols].reset_index(drop=True)

    df_main.to_excel(OUT_MAIN, index=False)
    df_expl.to_excel(OUT_EXPL, index=False)
    print(f'\n已写入: {OUT_MAIN.relative_to(ROOT)}  ({df_main.shape})')
    print(f'已写入: {OUT_EXPL.relative_to(ROOT)}  ({df_expl.shape})')

    # 审计表
    audit = pd.DataFrame({
        '代谢物': df['Metabolite Name'],
        '中文名': df['Chinese Name'],
        'Super Class': df['Super Class'],
        'Class': df['Class'],
        '检出率_正常': det_g1.round(4),
        '检出率_超重肥胖': det_g2.round(4),
        '最大组检出率': max_det.round(4),
        '总缺失率': is_missing.mean(axis=1).round(4),
        '主分析_80%': np.where(keep_main, '保留', '剔除'),
        '探索性_50%': np.where(keep_expl, '保留', '剔除'),
        '类别': np.select(
            [keep_main, keep_expl & ~keep_main, ~keep_expl],
            ['两档都保留', '仅探索性保留', '两档都剔除'],
            default='?'
        ),
    })
    audit.to_csv(OUT_AUDIT, index=False, encoding='utf-8-sig')
    print(f'已写入: {OUT_AUDIT.relative_to(ROOT)}')

    # 仅探索性保留的分子明细 (常被审稿人问到)
    explor_only = audit[audit['类别'] == '仅探索性保留']
    print(f'\n--- 仅探索性 (50%) 保留的 {len(explor_only)} 个 LOD 边缘分子 ---')
    print(explor_only[['代谢物', '中文名', '检出率_正常', '检出率_超重肥胖']].to_string(index=False))

    # 剔除后两组的残余缺失率 (插补可靠性预估)
    res_main = is_missing.loc[keep_main, sample_cols].values.mean()
    res_expl = is_missing.loc[keep_expl, sample_cols].values.mean()
    print(f'\n--- 过滤后残余全局缺失率 (用于评估插补难度) ---')
    print(f'  80% 主分析数据: {res_main*100:.2f}%   (Do 2018 要求 <50% 才插补)')
    print(f'  50% 探索数据:   {res_expl*100:.2f}%')


if __name__ == '__main__':
    main()
