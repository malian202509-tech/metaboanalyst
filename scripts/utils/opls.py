"""OPLS-DA (Trygg & Wold 2002) 完整实现.

适用: 二分类响应 (y in {-1, +1} 或 {0, 1}).
- 单 predictive 成分 + N_ortho orthogonal 成分.
- 模型存全部参数 (含均值/方差/旋转矩阵), 可用 opls_predict 在新数据上预测.
- 支持 k-fold CV 计算 Q², SIMCA-P 风格置换检验.

代谢组学惯例报告:
- R²X_pred, R²X_ortho, R²Y, Q²
- 置换检验 (200+ 次): R²Y/Q² 实测值在置换分布中的 P 值
- 置换图截距判定: R²Y intercept < 0.4 且 Q² intercept < 0.05 -> 模型可信
"""
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import KFold


def _norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def opls_da(X, y, n_ortho=1, scale_X=True, eps=1e-12):
    """OPLS-DA fit. 返回包含完整模型参数的 dict (可用 opls_predict 应用到新数据).

    Args:
        X: ndarray (n_samples, n_features)
        y: ndarray (n_samples,) 二分类 0/1 或 -1/+1
        n_ortho: orthogonal 成分数

    Returns dict keys:
        t_pred, t_ortho, w_pred, p_pred, c_pred, w_ortho_list, p_ortho_list,
        r2x_pred, r2x_ortho, r2y, y_hat,
        X_mean, X_sd, y_mean (用于 predict)
    """
    X = np.asarray(X, dtype=float)
    y_arr = np.asarray(y, dtype=float).flatten()
    X_mean = X.mean(axis=0)
    Xc = X - X_mean
    if scale_X:
        X_sd = Xc.std(axis=0, ddof=1)
        X_sd = np.where(X_sd > eps, X_sd, 1.0)
        Xc = Xc / X_sd
    else:
        X_sd = np.ones(X.shape[1])
    y_mean = y_arr.mean()
    yc = (y_arr - y_mean).reshape(-1, 1)

    total_var_X = (Xc ** 2).sum()
    ss_y = float((yc ** 2).sum())
    X_work = Xc.copy()

    w_pred = _norm((X_work.T @ yc).flatten())

    w_ortho_list, p_ortho_list, t_ortho_list = [], [], []
    for _ in range(n_ortho):
        t = X_work @ w_pred
        p = (X_work.T @ t) / (t @ t + eps)
        w_ortho = p - ((w_pred @ p) / (w_pred @ w_pred + eps)) * w_pred
        w_ortho = _norm(w_ortho)
        t_ortho = X_work @ w_ortho
        p_ortho = (X_work.T @ t_ortho) / (t_ortho @ t_ortho + eps)
        X_work = X_work - np.outer(t_ortho, p_ortho)
        w_ortho_list.append(w_ortho)
        p_ortho_list.append(p_ortho)
        t_ortho_list.append(t_ortho)
        w_pred = _norm((X_work.T @ yc).flatten())

    t_pred = X_work @ w_pred
    p_pred = (X_work.T @ t_pred) / (t_pred @ t_pred + eps)
    c_pred_arr = (yc.T @ t_pred) / (t_pred @ t_pred + eps)
    c_pred = float(np.asarray(c_pred_arr).flatten()[0])
    y_hat = t_pred * c_pred + y_mean

    r2x_pred = float((t_pred @ t_pred) * (p_pred @ p_pred) / total_var_X) if total_var_X > 0 else 0
    r2x_ortho = float(sum((to @ to) * (po @ po) for to, po in zip(t_ortho_list, p_ortho_list)) / total_var_X) if total_var_X > 0 else 0
    r2y = float(1 - ((y_arr - y_hat) ** 2).sum() / ss_y) if ss_y > 0 else 0

    return {
        't_pred': t_pred.flatten(),
        't_ortho': np.column_stack(t_ortho_list) if t_ortho_list else np.zeros((len(y_arr), 0)),
        'w_pred': w_pred,
        'p_pred': p_pred,
        'c_pred': c_pred,
        'w_ortho_list': w_ortho_list,
        'p_ortho_list': p_ortho_list,
        'r2x_pred': r2x_pred,
        'r2x_ortho': r2x_ortho,
        'r2y': r2y,
        'y_hat': y_hat,
        'X_mean': X_mean,
        'X_sd': X_sd,
        'y_mean': y_mean,
    }


def opls_predict(model, X_new):
    """用训练好的 OPLS-DA 模型预测新样本的 y."""
    X_c = (np.asarray(X_new, dtype=float) - model['X_mean']) / model['X_sd']
    for w_o, p_o in zip(model['w_ortho_list'], model['p_ortho_list']):
        t_o = X_c @ w_o
        X_c = X_c - np.outer(t_o, p_o)
    t_pred = X_c @ model['w_pred']
    return t_pred * model['c_pred'] + model['y_mean']


def opls_da_q2(X, y, n_ortho=1, n_folds=7, random_state=42):
    """OPLS-DA 的 Q² 经 k-fold 交叉验证.
    Q² = 1 - PRESS / SS_total, PRESS 用 out-of-fold 预测."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).flatten()
    n = len(y)
    n_folds = min(n_folds, n)
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=random_state)
    y_pred = np.zeros_like(y)
    for tr, te in kf.split(X):
        if len(np.unique(y[tr])) < 2:
            y_pred[te] = y[tr].mean()
            continue
        model = opls_da(X[tr], y[tr], n_ortho=n_ortho)
        y_pred[te] = opls_predict(model, X[te])
    ss_t = ((y - y.mean()) ** 2).sum()
    press = ((y - y_pred) ** 2).sum()
    return float(1 - press / ss_t) if ss_t > 0 else 0


def permutation_test_opls(X, y, n_ortho=1, n_perm=200, n_folds=7, random_state=42):
    """OPLS-DA SIMCA-P 风格置换检验.

    返回 dict:
      r2y_obs, q2_obs           - 实测
      r2y_perm, q2_perm         - 各 n_perm 次置换后的 R²Y / Q²
      correlations              - 置换 y 与原 y 的 |相关系数|
      p_r2y, p_q2               - 置换 P 值 (置换分布中 >= 实测的比例)
      r2y_intercept, q2_intercept - 在 |corr|=0 处的截距 (SIMCA 判定阈值)
    """
    rng = np.random.default_rng(random_state)
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).flatten()

    r2y_obs = opls_da(X, y, n_ortho=n_ortho)['r2y']
    q2_obs = opls_da_q2(X, y, n_ortho=n_ortho, n_folds=n_folds, random_state=random_state)

    r2y_perm, q2_perm, correlations = [], [], []
    for _ in range(n_perm):
        y_shuf = rng.permutation(y)
        try:
            r2y_perm.append(opls_da(X, y_shuf, n_ortho=n_ortho)['r2y'])
            q2_perm.append(opls_da_q2(X, y_shuf, n_ortho=n_ortho, n_folds=n_folds, random_state=42))
            correlations.append(abs(np.corrcoef(y_shuf, y)[0, 1]))
        except Exception:
            pass

    r2y_perm = np.asarray(r2y_perm)
    q2_perm = np.asarray(q2_perm)
    correlations = np.asarray(correlations)

    p_r2y = float((np.sum(r2y_perm >= r2y_obs) + 1) / (len(r2y_perm) + 1))
    p_q2 = float((np.sum(q2_perm >= q2_obs) + 1) / (len(q2_perm) + 1))

    # SIMCA-P 截距: 把 (corr=1, 实测) 和 (corr=|corr_i|, 置换值) 一起线性回归
    x_all = np.concatenate([correlations, [1.0]])
    y_all_r2y = np.concatenate([r2y_perm, [r2y_obs]])
    y_all_q2 = np.concatenate([q2_perm, [q2_obs]])
    slope_r2y, intercept_r2y = np.polyfit(x_all, y_all_r2y, 1)
    slope_q2, intercept_q2 = np.polyfit(x_all, y_all_q2, 1)

    return {
        'r2y_obs': float(r2y_obs),
        'q2_obs': float(q2_obs),
        'r2y_perm': r2y_perm,
        'q2_perm': q2_perm,
        'correlations': correlations,
        'p_r2y': p_r2y,
        'p_q2': p_q2,
        'r2y_intercept': float(intercept_r2y),
        'q2_intercept': float(intercept_q2),
    }


def pls_da(X, y, n_components=2):
    """普通 PLS-DA (sklearn 包装)."""
    y_bin = np.asarray(y).astype(float).reshape(-1, 1)
    pls = PLSRegression(n_components=n_components, scale=True)
    pls.fit(X, y_bin)
    t = pls.x_scores_   # n_samp × n_components
    y_pred = pls.predict(X).flatten()
    ss_t = ((y_bin - y_bin.mean()) ** 2).sum()
    ss_r = ((y_bin.flatten() - y_pred) ** 2).sum()
    r2y = float(1 - ss_r / ss_t)
    # R²X by component
    Xc = X - X.mean(axis=0)
    total_var = (Xc ** 2).sum()
    r2x = []
    for k in range(n_components):
        p = pls.x_loadings_[:, k]
        tk = pls.x_scores_[:, k]
        r2x.append(float((tk @ tk) * (p @ p) / total_var))
    return {
        't': t,
        'r2x_per_comp': r2x,
        'r2y': r2y,
        'y_pred': y_pred,
        'pls_obj': pls,
    }


def permutation_test_pls(X, y, n_components=2, n_perm=200, random_state=42):
    """置换检验 PLS-DA 的 R²Y. 返回 (R²Y_obs, p_value, perm_distribution)."""
    rng = np.random.default_rng(random_state)
    obs = pls_da(X, y, n_components)['r2y']
    perms = []
    for _ in range(n_perm):
        y_shuf = rng.permutation(y)
        try:
            perms.append(pls_da(X, y_shuf, n_components)['r2y'])
        except Exception:
            pass
    perms = np.asarray(perms)
    p = float((np.sum(perms >= obs) + 1) / (len(perms) + 1))
    return obs, p, perms
