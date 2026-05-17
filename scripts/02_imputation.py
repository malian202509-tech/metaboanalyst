"""02. 缺失值插补 (双轨: 80% 主分析 + 50% 探索性).

依据:
  - Wei et al. Sci Rep 2018: 左截断 MNAR (LOD truncation) 推荐 QRILC
  - Hrydziuszko & Viant Metabolomics 2012: MCAR 推荐 KNN
  - Do et al. Metabolomics 2018: 缺失率 <50% 时插补功效可靠

流程 (每条数据轨独立):
  Step 1. 缺失机制诊断
    - 计算每个特征的 (缺失率, 观测中位数) 的 Spearman 相关
    - rho < -0.3 -> 高缺失特征对应低丰度 -> MNAR (LOD 截断) 主导
    - 决定主插补方法

  Step 2. 主插补 (per-sample QRILC)
    - log2 转换 -> 估计每个样本的 N(mu, sigma) -> 在 1% 分位左侧抽样填补
    - exp2 还原到原始尺度

  Step 3. 验证
    - 输出: 无 NaN, 插补值范围合理 (不超出观测最小值)

输入:
  data/02_preprocessed/ori_n165_filtered80.xlsx
  data/02_preprocessed/ori_n165_filtered50.xlsx

输出:
  data/02_preprocessed/ori_n165_filtered80_imputed.xlsx
  data/02_preprocessed/ori_n165_filtered50_imputed.xlsx
  data/02_preprocessed/imputation_audit.csv
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import truncnorm, spearmanr

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / 'data' / '02_preprocessed'

META_COLS = ['Metabolite Name', 'Chinese Name', 'Retention Time(min)', 'KEGG ID',
             'HMDB ID', 'Cellular Locations', 'Super Class', 'Class', 'Sub Class',
             'Direct Parent']
UNIT_COL = '浓度单位'
QRILC_Q = 0.01
RANDOM_STATE = 42


def diagnose_mnar(matrix):
    """诊断缺失机制: 特征缺失率 vs 观测中位数的 Spearman 相关.
    强负相关 (rho < -0.3) -> 低丰度特征更容易缺失 -> MNAR (LOD 截断)."""
    miss_rate = matrix.isna().mean(axis=1)
    median_obs = matrix.median(axis=1, skipna=True)
    valid = ~median_obs.isna() & (miss_rate > 0) & (miss_rate < 1)
    if valid.sum() < 5:
        return np.nan, np.nan, '样本太少'
    rho, pval = spearmanr(miss_rate[valid], median_obs[valid])
    verdict = 'MNAR 强证据' if rho < -0.3 else ('MNAR 弱证据' if rho < 0 else '无 MNAR 信号 (近 MCAR)')
    return float(rho), float(pval), verdict


def qrilc_impute_per_sample(matrix, q=QRILC_Q, random_state=RANDOM_STATE):
    """QRILC 逐样本插补 (Wei 2018; R imputeLCMD 包).

    在 log2 尺度上: 对每个样本估计 N(mu, sigma), 在 q 分位左侧的截断正态分布抽样.
    返回原始尺度结果.
    """
    rng = np.random.default_rng(random_state)
    raw = matrix.astype(float).values
    log_data = np.log2(np.where(raw > 0, raw, np.nan))  # 0 或负值视为 NaN
    imputed_log = log_data.copy()
    n_imputed = 0
    sample_diag = []

    for j in range(log_data.shape[1]):
        x = log_data[:, j]
        mask_obs = ~np.isnan(x)
        x_obs = x[mask_obs]
        n_miss = int((~mask_obs).sum())
        if n_miss == 0:
            sample_diag.append({'col_idx': j, 'n_miss': 0, 'method': 'no_missing'})
            continue
        if len(x_obs) < 5:
            # fallback: feature 最小值 / 2 (log 尺度上 -1)
            fill_val = x_obs.min() - 1 if len(x_obs) > 0 else 0
            imputed_log[~mask_obs, j] = fill_val
            n_imputed += n_miss
            sample_diag.append({'col_idx': j, 'n_miss': n_miss, 'method': 'half_min_fallback'})
            continue

        mu = float(np.mean(x_obs))
        sigma = float(np.std(x_obs, ddof=1))
        thresh = float(np.quantile(x_obs, q))
        # 右截断: b = (thresh - mu) / sigma; 左侧无限
        b_std = (thresh - mu) / sigma if sigma > 0 else 0
        if sigma <= 0:
            imputed_log[~mask_obs, j] = mu
            n_imputed += n_miss
            sample_diag.append({'col_idx': j, 'n_miss': n_miss, 'method': 'sigma_zero_mu'})
            continue

        samples = truncnorm.rvs(
            a=-np.inf, b=b_std, loc=mu, scale=sigma,
            size=n_miss, random_state=int(rng.integers(0, 2**31 - 1))
        )
        imputed_log[~mask_obs, j] = samples
        n_imputed += n_miss
        sample_diag.append({'col_idx': j, 'n_miss': n_miss, 'method': 'qrilc',
                            'mu_log2': mu, 'sigma_log2': sigma, 'thresh_log2': thresh})

    imputed_raw = np.power(2.0, imputed_log)
    return pd.DataFrame(imputed_raw, index=matrix.index, columns=matrix.columns), n_imputed, sample_diag


def process_one(src_name, tag):
    """对单一过滤数据集执行完整插补流程."""
    print(f'\n{"="*60}')
    print(f'轨道 {tag}: {src_name}')
    src = OUT_DIR / src_name
    df = pd.read_excel(src)
    sample_cols = [c for c in df.columns if c not in META_COLS and c != UNIT_COL]
    matrix = df[sample_cols].copy()
    n_feat, n_samp = matrix.shape
    print(f'  形状: {n_feat} 代谢物 x {n_samp} 样本')
    print(f'  缺失数: {int(matrix.isna().sum().sum())}  ({matrix.isna().values.mean()*100:.2f}%)')

    # === Step 1: 诊断 ===
    rho, pval, verdict = diagnose_mnar(matrix)
    print(f'  --- MNAR 诊断 (Spearman 缺失率 vs 观测中位数) ---')
    print(f'    rho = {rho:.3f}  p = {pval:.3g}  -> {verdict}')

    # === Step 2: QRILC 插补 ===
    print(f'  --- QRILC 插补 (per-sample, log2 尺度) ---')
    imputed, n_imp, diag = qrilc_impute_per_sample(matrix)
    print(f'    实际插补值数: {n_imp}')

    # === Step 3: 验证 ===
    n_remain_nan = int(imputed.isna().sum().sum())
    n_neg = int((imputed < 0).sum().sum())
    obs_min = matrix.min().min()
    imp_only_mask = matrix.isna() & ~imputed.isna()
    imp_values = imputed.values[imp_only_mask.values]
    print(f'  --- 验证 ---')
    print(f'    残留 NaN: {n_remain_nan}  (应为 0)')
    print(f'    负值: {n_neg}  (应为 0)')
    print(f'    观测最小值: {obs_min:.4g}')
    print(f'    插补值最小: {imp_values.min():.4g}   最大: {imp_values.max():.4g}')
    print(f'    插补值中位数: {np.median(imp_values):.4g}   (理应 << 观测中位数)')
    print(f'    观测中位数: {matrix.median().median():.4g}')

    # 落盘
    df_out = df.copy()
    df_out[sample_cols] = imputed.values
    out_name = src_name.replace('.xlsx', '_imputed.xlsx')
    out_path = OUT_DIR / out_name
    df_out.to_excel(out_path, index=False)
    print(f'  已写入: {out_path.relative_to(ROOT)}')

    return {
        'tier': tag,
        'src': src_name,
        'out': out_name,
        'n_feat': n_feat,
        'n_samp': n_samp,
        'pre_missing_pct': round(matrix.isna().values.mean() * 100, 3),
        'n_imputed': n_imp,
        'mnar_rho': round(rho, 3),
        'mnar_pval': round(pval, 4),
        'mnar_verdict': verdict,
        'obs_min': float(obs_min),
        'imp_min': float(imp_values.min()),
        'imp_max': float(imp_values.max()),
        'imp_median': float(np.median(imp_values)),
        'obs_median': float(matrix.median().median()),
        'method': f'QRILC q={QRILC_Q}, seed={RANDOM_STATE}',
    }


def main():
    print('=== 02. 缺失值插补 (双轨并行) ===')
    records = []
    records.append(process_one('ori_n165_filtered80.xlsx', '80% 主分析'))
    records.append(process_one('ori_n165_filtered50.xlsx', '50% 探索性'))

    audit = pd.DataFrame(records)
    audit_path = OUT_DIR / 'imputation_audit.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    print(f'\n审计表已写入: {audit_path.relative_to(ROOT)}')
    print()
    print('--- 汇总 ---')
    print(audit[['tier', 'n_feat', 'pre_missing_pct', 'n_imputed',
                 'mnar_rho', 'mnar_verdict', 'imp_median', 'obs_median']].to_string(index=False))


if __name__ == '__main__':
    main()
