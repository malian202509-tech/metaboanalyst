"""Empirical Bayes moderated t-statistics (Smyth, Stat Appl Genet Mol Biol 2004).

跨特征收缩残差方差, 实现 limma 风格的差异分析:
  s²_post = (d_0 · s²_0 + d_g · s²_g) / (d_0 + d_g)
  t_mod   = β̂ / (s_post · √(c' (X'X)⁻¹ c))
  df      = d_0 + d_g

参考: Smyth GK (2004). Linear models and empirical Bayes methods for assessing
      differential expression in microarray experiments. Stat Appl Genet Mol Biol 3:3.

适用代谢组学小样本场景: 特征数 67/73 时, 仅靠单特征方差不稳, EB 收缩显著提升功效.
"""
import numpy as np
from scipy import stats, special


def _trigamma_inv(x, max_iter=50, tol=1e-8):
    """trigamma(y) = x 的 Newton 反解 (Smyth 2004 附录算法).

    步长公式与 limma::trigammaInverse 一致 (R 源代码):
        dif = trigamma(y) * (1 - trigamma(y)/x) / psigamma(y, deriv=2)
        y  <- y + dif      # 注意是 + 不是 -
    psigamma(y, deriv=2) < 0, 与分子配合保证 y 朝正确方向迭代.
    """
    if not np.isfinite(x) or x <= 0:
        return np.inf
    y = 0.5 + 1.0 / x if x <= 1e7 else 1.0 / np.sqrt(x)
    for _ in range(max_iter):
        tri = special.polygamma(1, y)
        psipp = special.polygamma(2, y)
        dy = tri * (1.0 - tri / x) / psipp
        y_new = y + dy
        if y_new <= 0:
            y_new = y / 2.0
        if abs(dy / y_new) < tol:
            y = y_new
            break
        y = y_new
    return y


def fit_prior(s2, df):
    """方法估计先验 (d_0, s²_0).

    Args:
        s2: (n_features,) 各特征残差方差 SSE/df
        df: 标量或 (n_features,) 残差 df

    Returns:
        (d0, s2_0): 标量先验
    """
    s2 = np.asarray(s2, dtype=float)
    df_arr = np.broadcast_to(np.asarray(df, dtype=float), s2.shape)
    mask = (s2 > 0) & np.isfinite(s2) & (df_arr > 0)
    s2v, dfv = s2[mask], df_arr[mask]
    if s2v.size < 2:
        return np.inf, float(np.nanmean(s2)) if s2v.size else 1.0

    z = np.log(s2v)
    e = z - special.digamma(dfv / 2) + np.log(dfv / 2)
    e_mean = float(np.mean(e))
    e_var = float(np.var(e, ddof=1))

    target = e_var - float(np.mean(special.polygamma(1, dfv / 2)))
    if target <= 0:
        return np.inf, float(np.exp(e_mean))

    d0_half = _trigamma_inv(target)
    if not np.isfinite(d0_half):
        return np.inf, float(np.exp(e_mean))
    d0 = 2.0 * d0_half
    s2_0 = float(np.exp(e_mean + special.digamma(d0_half) - np.log(d0_half)))
    return d0, s2_0


def moderated_t(beta, se_ord_factor, s2, df):
    """对 β̂ 应用 EB 方差收缩, 返回 ordinary + moderated 统计量.

    Args:
        beta: (n_features,) 单个对比的系数估计
        se_ord_factor: 标量, c' (X'X)⁻¹ c (跨特征常数, 因为同一设计矩阵)
        s2: (n_features,) 残差方差
        df: 标量或 (n_features,) 残差 df

    Returns:
        dict: beta, se_ord, se_mod, t_ord, t_mod, p_ord, p_mod,
              df_prior, s2_prior, df_total (标量或数组), s2_post
    """
    beta = np.asarray(beta, dtype=float)
    s2 = np.asarray(s2, dtype=float)
    df_arr = np.broadcast_to(np.asarray(df, dtype=float), s2.shape).astype(float)

    d0, s2_0 = fit_prior(s2, df_arr)

    if np.isinf(d0):
        s2_post = np.full_like(s2, s2_0)
        df_total = df_arr.copy()
    else:
        s2_post = (d0 * s2_0 + df_arr * s2) / (d0 + df_arr)
        df_total = df_arr + d0

    se_ord = np.sqrt(np.maximum(s2 * se_ord_factor, 0.0))
    se_mod = np.sqrt(np.maximum(s2_post * se_ord_factor, 0.0))

    with np.errstate(divide='ignore', invalid='ignore'):
        t_ord = np.where(se_ord > 0, beta / se_ord, np.nan)
        t_mod = np.where(se_mod > 0, beta / se_mod, np.nan)
        p_ord = 2.0 * stats.t.sf(np.abs(t_ord), df=df_arr)
        p_mod = 2.0 * stats.t.sf(np.abs(t_mod), df=df_total)

    return {
        'beta': beta, 's2': s2, 's2_post': s2_post,
        'se_ord': se_ord, 'se_mod': se_mod,
        't_ord': t_ord, 't_mod': t_mod,
        'p_ord': p_ord, 'p_mod': p_mod,
        'df_prior': d0, 's2_prior': s2_0, 'df_total': df_total,
    }


def bh_qvalues(pvals):
    """Benjamini-Hochberg q-values (与 R p.adjust(..., method='BH') 一致)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    out = np.full(n, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    m = pv.size
    if m == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * m / (np.arange(1, m + 1))
    # 单调强制 (从右向左取累积最小)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out_valid = np.empty(m)
    out_valid[order] = q
    out[valid] = out_valid
    return out
