"""
Between-study variance (tau^2) estimators for the second stage of IPD-QMA.

The shipped :class:`ipd_qma.IPDQMA` engine pools study-level quantile treatment
effects with the DerSimonian-Laird (DL) estimator. The IPD-QMA manuscript notes
(Limitations / Future directions) that a sensitivity analysis comparing DL,
restricted maximum likelihood (REML), and Paule-Mandel (PM) estimators "was not
performed and is left for future work"; the pre-registered protocol in fact names
REML as the intended primary estimator.

This module implements that documented-but-unbuilt feature as an ADDITIVE,
standalone layer:

  * ``dl_tau2``   - reproduces the engine's exact DL moment estimator (cross-check).
  * ``pm_tau2``   - Paule-Mandel / empirical-Bayes iterative moment estimator.
  * ``reml_tau2`` - restricted maximum likelihood (Viechtbauer 2005 iteration).
  * ``pool_random_effects`` - inverse-variance RE pool given any tau^2.
  * ``tau2_sensitivity`` - all three, side by side, on the same (est, se) inputs.

It does NOT modify the shipped engine and does NOT change the published DL-based
results. Validated against R ``metafor::rma`` to 1e-6 (see test_tau2_estimators.py).

Formulas / gotchas honored (per project stats rules):
  * DL for k < 10 is biased low for tau^2; REML/PM reduce that bias.
  * All estimators clamp tau^2 at 0 (no negative variance).
  * Homogeneous data (Q <= k-1) => tau^2 = 0 for all three.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy import stats

_EPS = 1e-12


def _clean(est, se):
    est = np.asarray(est, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    if est.shape != se.shape:
        raise ValueError("est and se must have the same shape")
    if est.ndim != 1:
        raise ValueError("est and se must be 1-D arrays")
    if len(est) < 1:
        raise ValueError("need at least 1 study")
    if np.any(~np.isfinite(est)) or np.any(~np.isfinite(se)):
        raise ValueError("est and se must be finite")
    if np.any(se <= 0):
        raise ValueError("all standard errors must be strictly positive")
    return est, se, len(est)


def _weighted_mean(est, w):
    return float(np.sum(est * w) / np.sum(w))


def _pm_residual(t2, est, se, k):
    """Paule-Mandel weighted residual sum of squares minus its expectation (k-1).

    Monotonically decreasing in ``t2``; the PM estimate is its positive root.
    """
    w = 1.0 / (se ** 2 + t2)
    mu = _weighted_mean(est, w)
    return float(np.sum(w * (est - mu) ** 2) - (k - 1))


def dl_tau2(est, se) -> float:
    """DerSimonian-Laird moment estimator of tau^2 (clamped at 0).

    Identical formula to ``IPDQMA._pool_dl`` so it can serve as a cross-check.
    """
    est, se, k = _clean(est, se)
    if k == 1:
        return 0.0
    w = 1.0 / se ** 2
    mu = _weighted_mean(est, w)
    Q = float(np.sum(w * (est - mu) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    if c <= 0:
        return 0.0
    return max(0.0, (Q - (k - 1)) / c)


def pm_tau2(est, se, tol: float = 1e-10, max_iter: int = 200) -> float:
    """Paule-Mandel (empirical-Bayes) estimator via monotone bisection.

    Solves ``sum w_i (y_i - mu_w)^2 = k - 1`` with ``w_i = 1/(se_i^2 + tau^2)``.
    """
    est, se, k = _clean(est, se)
    if k == 1:
        return 0.0
    # If even tau^2 = 0 already under-shoots the target, heterogeneity is 0.
    if _pm_residual(0.0, est, se, k) <= 0:
        return 0.0
    lo, hi = 0.0, 1.0
    # Expand the bracket until the residual turns negative (it must: as tau^2 -> inf
    # the weighted RSS -> 0 < k-1). Cap the expansion defensively.
    guard = 0
    while _pm_residual(hi, est, se, k) > 0 and guard < 200:
        hi *= 2.0
        guard += 1
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        r = _pm_residual(mid, est, se, k)
        if abs(r) < tol or (hi - lo) < tol:
            return mid
        if r > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def reml_tau2(est, se, tol: float = 1e-10, max_iter: int = 1000) -> float:
    """Restricted maximum likelihood estimator (Viechtbauer 2005 fixed-point).

    Iterates ``tau2 <- sum w^2[(y-mu)^2 + 1/sum(w) - se^2] / sum(w^2)`` with
    ``w = 1/(se^2 + tau2)`` and ``mu`` the current RE-weighted mean; clamps at 0.
    Seeded from the DL estimate for fast, stable convergence.
    """
    est, se, k = _clean(est, se)
    if k == 1:
        return 0.0
    t2 = dl_tau2(est, se)
    for _ in range(max_iter):
        w = 1.0 / (se ** 2 + t2)
        sw = float(np.sum(w))
        mu = _weighted_mean(est, w)
        num = float(np.sum(w ** 2 * ((est - mu) ** 2 + 1.0 / sw - se ** 2)))
        den = float(np.sum(w ** 2))
        new = max(0.0, num / den)
        if abs(new - t2) < tol:
            t2 = new
            break
        t2 = new
    return max(0.0, t2)


_ESTIMATORS = {"DL": dl_tau2, "REML": reml_tau2, "PM": pm_tau2}


@dataclass
class REPool:
    method: str
    tau2: float
    estimate: float
    se: float
    ci_lower: float
    ci_upper: float
    I2: float
    Q: float
    k: int


def pool_random_effects(est, se, tau2: float, conf_level: float = 0.95,
                        use_hksj: bool = True, method: str = "") -> REPool:
    """Inverse-variance random-effects pool given a fixed ``tau2``.

    Mirrors the inference in ``IPDQMA._pool_dl`` (HKSJ ``max(1, Q*/(k-1))`` floor
    with ``t_{k-1}`` when ``use_hksj`` and ``k >= 3``; normal otherwise) so a
    sensitivity table stays comparable to the shipped engine's DL cell.
    """
    est, se, k = _clean(est, se)
    alpha = 1.0 - conf_level
    if k == 1:
        z = stats.norm.ppf(1 - alpha / 2)
        s = max(float(se[0]), _EPS)
        return REPool(method or "-", 0.0, float(est[0]), s,
                      float(est[0]) - z * s, float(est[0]) + z * s, 0.0, 0.0, 1)

    w_fe = 1.0 / se ** 2
    mu_fe = _weighted_mean(est, w_fe)
    Q = float(np.sum(w_fe * (est - mu_fe) ** 2))
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0

    w_re = 1.0 / (se ** 2 + tau2)
    p_est = _weighted_mean(est, w_re)
    p_se_iv = float(np.sqrt(1.0 / np.sum(w_re)))

    if use_hksj and k >= 3:
        q_hksj = float(np.sum(w_re * (est - p_est) ** 2) / (k - 1))
        p_se = p_se_iv * np.sqrt(max(1.0, q_hksj))
        crit = stats.t.ppf(1 - alpha / 2, df=k - 1)
    else:
        p_se = p_se_iv
        crit = stats.norm.ppf(1 - alpha / 2)

    return REPool(method or "-", float(tau2), float(p_est), float(p_se),
                  float(p_est - crit * p_se), float(p_est + crit * p_se),
                  float(I2), Q, k)


def tau2_sensitivity(est, se, conf_level: float = 0.95,
                     use_hksj: bool = True) -> Dict[str, REPool]:
    """DL / REML / PM sensitivity table on the same study-level inputs.

    Returns a dict keyed by estimator name; each value is the RE pool under that
    estimator's tau^2. DL is the shipped default and is reported first.
    """
    out: Dict[str, REPool] = {}
    for name, fn in _ESTIMATORS.items():
        t2 = fn(est, se)
        out[name] = pool_random_effects(est, se, t2, conf_level=conf_level,
                                        use_hksj=use_hksj, method=name)
    return out


if __name__ == "__main__":  # pragma: no cover - manual smoke
    y = [0.20, 0.80, -0.10, 0.55, 0.35, 0.90, 0.05, 0.65]
    s = [0.10, 0.12, 0.11, 0.13, 0.10, 0.14, 0.09, 0.12]
    print(f"{'method':>6} {'tau2':>10} {'estimate':>10} {'ci_lower':>10} {'ci_upper':>10} {'I2':>6}")
    for name, r in tau2_sensitivity(y, s).items():
        print(f"{name:>6} {r.tau2:>10.5f} {r.estimate:>10.5f} "
              f"{r.ci_lower:>10.5f} {r.ci_upper:>10.5f} {r.I2:>5.1f}%")
