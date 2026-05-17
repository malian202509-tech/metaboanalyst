"""Parametric Empirical Bayes ComBat (Johnson, Li & Rabinovic 2007).

参考: sva::ComBat R 实现 (Leek lab).
- 行=特征 (代谢物), 列=样本.
- batch: 每个样本的批次标签
- mod:   生物学协变量矩阵 (n_samples × p), 不含截距, 通常含分组哑变量
        传 mod 后, ComBat 会保护这些因子的方差不被抹掉.
- 默认参数先验 (par_prior=True).
"""
import numpy as np
import pandas as pd


def _it_sol(gamma_hat, delta_sq_hat, n_g, gamma_g_bar, tau_g_sq, lambda_g, theta_g,
            conv=1e-4, max_iter=200):
    """EB 后验的固定点迭代 (sva 的 it.sol)."""
    g_old = gamma_hat.copy()
    d_old = delta_sq_hat.copy()
    for _ in range(max_iter):
        g_new = (n_g * tau_g_sq * g_old + d_old * gamma_g_bar) / (n_g * tau_g_sq + d_old)
        ss = (n_g - 1) * delta_sq_hat + n_g * (gamma_hat - g_new) ** 2
        d_new = (theta_g + 0.5 * ss) / (n_g / 2 + lambda_g - 1)
        change = max(
            np.max(np.abs(g_new - g_old) / (np.abs(g_old) + 1e-12)),
            np.max(np.abs(d_new - d_old) / (np.abs(d_old) + 1e-12)),
        )
        g_old, d_old = g_new, d_new
        if change < conv:
            break
    return g_old, d_old


def combat(data, batch, mod=None, mean_only=False, eb=True):
    """参数 EB ComBat.

    Args:
        data: pd.DataFrame, rows=features (代谢物), cols=samples
        batch: array-like, 每个样本的批次, 长度 = n_samples
        mod: pd.DataFrame 或 ndarray, shape (n_samples, p), 生物学协变量(不含截距)
        mean_only: bool, 仅校正均值不校正方差 (默认 False)
        eb: bool, 是否做经验贝叶斯收缩 (默认 True). 小特征数 (<200) 建议设 False,
            因为 EB 假设有大量特征可估稳定先验, 否则容易过度收缩.

    Returns:
        pd.DataFrame, 与输入同形状, 已校正
    """
    Y = data.values.astype(float)          # n_feat × n_samp
    n_feat, n_samp = Y.shape
    batch = pd.Series(batch).reset_index(drop=True)
    batches = batch.unique()
    n_batches = len(batches)
    n_per_batch = np.array([(batch == b).sum() for b in batches])

    # 设计矩阵: batch one-hot (作为截距) + mod
    batch_oh = pd.get_dummies(batch, dtype=float).reindex(columns=batches).values  # n_samp × n_batches
    if mod is None:
        design = batch_oh
        mod_mat = None
    else:
        mod_mat = (mod.values if isinstance(mod, pd.DataFrame) else np.asarray(mod)).astype(float)
        if mod_mat.ndim == 1:
            mod_mat = mod_mat.reshape(-1, 1)
        design = np.hstack([batch_oh, mod_mat])

    # 单特征 OLS: B = (X'X)^-1 X'Y'
    XtX_inv = np.linalg.pinv(design.T @ design)
    B = XtX_inv @ design.T @ Y.T            # (n_design_cols, n_feat)

    # 加权"总均值"(各批次的截距按样本数加权)
    gamma_intercepts = B[:n_batches]        # n_batches × n_feat
    alpha_hat = (n_per_batch @ gamma_intercepts) / n_samp  # n_feat

    # 生物学拟合值 (不含批次)
    if mod_mat is None:
        mod_fitted = np.zeros_like(Y)
    else:
        mod_coefs = B[n_batches:]            # n_mod × n_feat
        mod_fitted = (mod_mat @ mod_coefs).T  # n_feat × n_samp

    # 池化残差方差
    fitted_full = (design @ B).T
    resid = Y - fitted_full
    sigma_sq = (resid ** 2).sum(axis=1) / n_samp
    sigma = np.sqrt(sigma_sq) + 1e-12

    # 标准化: 去掉总均值与生物学贡献后, 按特征 SD 缩放
    Z = (Y - mod_fitted - alpha_hat[:, None]) / sigma[:, None]

    # 各批次的批次效应 (γ_hat, δ²_hat)
    gamma_hat = np.zeros((n_batches, n_feat))
    delta_sq_hat = np.zeros((n_batches, n_feat))
    for g, b in enumerate(batches):
        Zg = Z[:, (batch == b).values]
        gamma_hat[g] = Zg.mean(axis=1)
        delta_sq_hat[g] = Zg.var(axis=1, ddof=1) + 1e-12

    # EB 先验 (跨特征)
    gamma_bar = gamma_hat.mean(axis=1)              # n_batches
    tau_sq = gamma_hat.var(axis=1, ddof=1)          # n_batches
    m = delta_sq_hat.mean(axis=1)
    v = delta_sq_hat.var(axis=1, ddof=1)
    lambda_hat = (2 * v + m ** 2) / np.where(v > 0, v, 1e-12)
    theta_hat = (m * v + m ** 3) / np.where(v > 0, v, 1e-12)

    # 后验
    if eb:
        # EB 收缩 (Johnson 2007 固定点迭代)
        gamma_star = np.zeros_like(gamma_hat)
        delta_sq_star = np.zeros_like(delta_sq_hat)
        for g in range(n_batches):
            gamma_star[g], delta_sq_star[g] = _it_sol(
                gamma_hat[g], delta_sq_hat[g], n_per_batch[g],
                gamma_bar[g], tau_sq[g], lambda_hat[g], theta_hat[g],
            )
    else:
        # 直接用经验估计 (适合小特征数, 避免 EB 过度收缩)
        gamma_star = gamma_hat.copy()
        delta_sq_star = delta_sq_hat.copy()

    # 校正并还原
    Y_adj = Y.copy()
    for g, b in enumerate(batches):
        mask = (batch == b).values
        Z_g = Z[:, mask]
        if mean_only:
            Z_g_adj = Z_g - gamma_star[g][:, None]
        else:
            Z_g_adj = (Z_g - gamma_star[g][:, None]) / np.sqrt(delta_sq_star[g][:, None])
        Y_adj[:, mask] = sigma[:, None] * Z_g_adj + alpha_hat[:, None] + mod_fitted[:, mask]

    return pd.DataFrame(Y_adj, index=data.index, columns=data.columns)
