"""05. ComBat 批次校正 (按采样年份, 保护 BMI_group).

依据: Johnson, Li & Rabinovic, Biostatistics 2007.
- batch = year (A19 / B20 / C21)
- mod   = BMI_group (one-hot 哑变量) → 防止 ComBat 抹掉生物学方差
- 输入: log2 矩阵 (ComBat 假设近似正态)
- 输出: 校正后的 log2 矩阵 + 同步生成 Pareto 缩放版供 PCA 使用

输入:
  data/02_preprocessed/ori_n165_filtered{80,50}_log2.xlsx
  data/02_preprocessed/sample_alignment_n165.csv

输出:
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat.xlsx          (ANCOVA/limma 输入)
  data/03_batch_corrected/ori_n165_filtered{80,50}_log2_combat_pareto.xlsx   (PCA/可视化输入)
  data/03_batch_corrected/combat_audit.csv
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
from utils.combat import combat

IN_DIR = ROOT / 'data' / '02_preprocessed'
OUT_DIR = ROOT / 'data' / '03_batch_corrected'
OUT_DIR.mkdir(parents=True, exist_ok=True)

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'


def pareto_scale(matrix):
    mean = matrix.mean(axis=1).values.reshape(-1, 1)
    sd = matrix.std(axis=1, ddof=1).values.reshape(-1, 1)
    sd_safe = np.where(sd > 0, sd, 1.0)
    return pd.DataFrame((matrix.values - mean) / np.sqrt(sd_safe),
                        index=matrix.index, columns=matrix.columns)


def process_one(src_name, tag):
    print(f'\n{"="*60}')
    print(f'轨道 {tag}: {src_name}')
    df = pd.read_excel(IN_DIR / src_name)
    mp = pd.read_csv(IN_DIR / 'sample_alignment_n165.csv')
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    matrix = df[sample_cols].copy()
    n_feat, n_samp = matrix.shape
    print(f'  形状: {n_feat} 代谢物 x {n_samp} 样本')

    # 元信息: year (batch) + BMI (mod)
    meta = pd.DataFrame({'omx_id': sample_cols})
    meta = meta.merge(mp[['omx_id', 'BMI_group']], on='omx_id', how='left')
    meta['year'] = meta['omx_id'].str[:3]
    print(f"  批次 (year): {meta['year'].value_counts().to_dict()}")
    print(f"  生物学保护 (BMI): {meta['BMI_group'].value_counts().to_dict()}")

    # mod = BMI 哑变量 (单列即可, 二分类)
    mod = pd.get_dummies(meta['BMI_group'], drop_first=True, dtype=float)

    # === ComBat ===
    # eb=False: 67/73 个特征太少, EB 会过度收缩 (诊断证实保留率 <50%).
    # 用经典 location-scale 校正 (直接用各批次均值/方差估计, 不做 EB 收缩).
    print(f'  --- 跑 ComBat (batch=year, mod=BMI_group, eb=False) ---')
    pre_var_total = float((matrix.values ** 2).sum())
    corrected = combat(matrix, batch=meta['year'].values, mod=mod.values, eb=False)
    post_var_total = float((corrected.values ** 2).sum())
    print(f'  全局方差: {pre_var_total:.2f} -> {post_var_total:.2f} '
          f'(变化 {(post_var_total - pre_var_total) / pre_var_total * 100:+.2f}%)')

    # 落盘: log2_combat
    df_log = df.copy()
    df_log[sample_cols] = corrected.values
    out_log = src_name.replace('_log2.xlsx', '_log2_combat.xlsx')
    df_log.to_excel(OUT_DIR / out_log, index=False)
    print(f'  已写入: {(OUT_DIR / out_log).relative_to(ROOT)}')

    # 落盘: log2_combat + Pareto (PCA 用)
    pareto = pareto_scale(corrected)
    df_par = df.copy()
    df_par[sample_cols] = pareto.values
    out_par = src_name.replace('_log2.xlsx', '_log2_combat_pareto.xlsx')
    df_par.to_excel(OUT_DIR / out_par, index=False)
    print(f'  已写入: {(OUT_DIR / out_par).relative_to(ROOT)}')

    return {
        'tier': tag,
        'src': src_name,
        'n_feat': n_feat,
        'n_samp': n_samp,
        'pre_total_var': round(pre_var_total, 2),
        'post_total_var': round(post_var_total, 2),
        'var_change_%': round((post_var_total - pre_var_total) / pre_var_total * 100, 2),
        'out_log': out_log,
        'out_pareto': out_par,
    }


def main():
    print('=== 05. ComBat 批次校正 (双轨并行) ===')
    records = []
    records.append(process_one('ori_n165_filtered80_log2.xlsx', '80% 主分析'))
    records.append(process_one('ori_n165_filtered50_log2.xlsx', '50% 探索性'))

    audit = pd.DataFrame(records)
    audit_path = OUT_DIR / 'combat_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
