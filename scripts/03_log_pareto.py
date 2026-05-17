"""03. log2 转换 + Pareto 缩放 (双轨并行).

依据:
  - log2: 代谢物浓度跨数个数量级, 残差严重偏态. log2 后近似正态, 可作 ANCOVA 因变量.
    系数解释直接是 log2 fold change.
  - Pareto scaling (van den Berg 2006, BMC Genomics):
      x' = (x - mean) / sqrt(SD)
    比 autoscaling (除以 SD) 温和, 不会让微量信号过度放大.
    适用 PCA / 热图 / 聚类.
  - 跳过 PQN: 单位已是 ng/g (归一化到组织湿重), 样本间归一化已生物学完成.

下游用法:
  - ANCOVA / limma 输入: 用 *_log2.xlsx (保留原始尺度信息, 系数 = log2 FC)
  - PCA / 热图 / 聚类输入: 用 *_log2_pareto.xlsx

输入:
  data/02_preprocessed/ori_n165_filtered80_imputed.xlsx
  data/02_preprocessed/ori_n165_filtered50_imputed.xlsx

输出:
  data/02_preprocessed/ori_n165_filtered80_log2.xlsx          (ANCOVA 输入, 主)
  data/02_preprocessed/ori_n165_filtered80_log2_pareto.xlsx   (PCA 输入, 主)
  data/02_preprocessed/ori_n165_filtered50_log2.xlsx          (ANCOVA 输入, 探索)
  data/02_preprocessed/ori_n165_filtered50_log2_pareto.xlsx   (PCA 输入, 探索)
  data/02_preprocessed/transform_audit.csv
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / 'data' / '02_preprocessed'

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'


def pareto_scale(matrix):
    """Feature-wise Pareto scaling: (x - mean) / sqrt(SD).
    matrix: rows=features, cols=samples.
    返回: 每行 mean=0, 每行 SD=sqrt(原SD).
    """
    mean = matrix.mean(axis=1).values.reshape(-1, 1)
    sd = matrix.std(axis=1, ddof=1).values.reshape(-1, 1)
    sd_safe = np.where(sd > 0, sd, 1.0)
    scaled = (matrix.values - mean) / np.sqrt(sd_safe)
    return pd.DataFrame(scaled, index=matrix.index, columns=matrix.columns)


def shapiro_pass_rate(matrix, alpha=0.05):
    """每个特征做 Shapiro-Wilk 正态性检验, 返回 P>alpha 的比例 (即'通过正态'的特征比例)."""
    passes = 0
    total = 0
    for i in range(matrix.shape[0]):
        x = matrix.iloc[i].dropna().values
        if len(x) < 3:
            continue
        try:
            p = stats.shapiro(x).pvalue
            if p > alpha:
                passes += 1
            total += 1
        except Exception:
            pass
    return passes, total, passes / total if total > 0 else np.nan


def process_one(src_name, tag):
    print(f'\n{"="*60}')
    print(f'轨道 {tag}: {src_name}')
    src = OUT_DIR / src_name
    df = pd.read_excel(src)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    matrix = df[sample_cols].copy()
    n_feat, n_samp = matrix.shape
    print(f'  形状: {n_feat} 代谢物 x {n_samp} 样本')

    # 健全性检查
    n_nan = int(matrix.isna().sum().sum())
    n_nonpos = int((matrix <= 0).sum().sum())
    print(f'  NaN: {n_nan}  | 非正值: {n_nonpos}  (log2 要求 > 0)')
    if n_nan > 0 or n_nonpos > 0:
        raise ValueError(f'  输入矩阵含 NaN 或非正值, 无法 log2. 请先跑插补.')

    # --- log2 ---
    log2_matrix = np.log2(matrix.astype(float))
    print(f'  log2 后: 范围 [{log2_matrix.values.min():.3f}, {log2_matrix.values.max():.3f}]')

    # 正态性: log2 前后比较
    _, _, pre_rate = shapiro_pass_rate(matrix)
    _, _, post_log_rate = shapiro_pass_rate(log2_matrix)
    print(f'  Shapiro-Wilk 通过率 (P>0.05): 原始 {pre_rate*100:.1f}% -> log2 {post_log_rate*100:.1f}%')

    # --- Pareto scaling ---
    pareto_matrix = pareto_scale(log2_matrix)
    row_means = pareto_matrix.mean(axis=1)
    row_sds = pareto_matrix.std(axis=1, ddof=1)
    print(f'  Pareto 后: 行均值 [{row_means.min():.2e}, {row_means.max():.2e}]  (应接近 0)')
    print(f'             行SD   [{row_sds.min():.3f}, {row_sds.max():.3f}]')

    # 落盘 (两个文件)
    df_log = df.copy()
    df_log[sample_cols] = log2_matrix.values
    out_log = src_name.replace('_imputed.xlsx', '_log2.xlsx')
    (OUT_DIR / out_log).parent.mkdir(parents=True, exist_ok=True)
    df_log.to_excel(OUT_DIR / out_log, index=False)
    print(f'  已写入 (log2 only):     {out_log}')

    df_par = df.copy()
    df_par[sample_cols] = pareto_matrix.values
    out_par = src_name.replace('_imputed.xlsx', '_log2_pareto.xlsx')
    df_par.to_excel(OUT_DIR / out_par, index=False)
    print(f'  已写入 (log2 + Pareto): {out_par}')

    return {
        'tier': tag,
        'src': src_name,
        'out_log2': out_log,
        'out_log2_pareto': out_par,
        'n_feat': n_feat,
        'n_samp': n_samp,
        'shapiro_pass_pre%': round(pre_rate * 100, 1),
        'shapiro_pass_log2%': round(post_log_rate * 100, 1),
        'log2_min': round(float(log2_matrix.values.min()), 3),
        'log2_max': round(float(log2_matrix.values.max()), 3),
        'pareto_row_mean_max_abs': float(row_means.abs().max()),
        'pareto_row_sd_min': float(row_sds.min()),
        'pareto_row_sd_max': float(row_sds.max()),
    }


def main():
    print('=== 03. log2 + Pareto (双轨并行) ===')
    records = []
    records.append(process_one('ori_n165_filtered80_imputed.xlsx', '80% 主分析'))
    records.append(process_one('ori_n165_filtered50_imputed.xlsx', '50% 探索性'))

    audit = pd.DataFrame(records)
    audit_path = OUT_DIR / 'transform_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表: {audit_path.relative_to(ROOT)}')
    print()
    print('--- 汇总 ---')
    print(audit[['tier', 'n_feat', 'shapiro_pass_pre%', 'shapiro_pass_log2%',
                 'log2_min', 'log2_max']].to_string(index=False))


if __name__ == '__main__':
    main()
