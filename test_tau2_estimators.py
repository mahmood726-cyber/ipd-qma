"""
Tests for tau2_estimators.py (DL / REML / Paule-Mandel sensitivity layer).

Gold-standard reference values were produced with R metafor 4.x:

    library(metafor)
    y <- c(0.20,0.80,-0.10,0.55,0.35,0.90,0.05,0.65)
    s <- c(0.10,0.12,0.11,0.13,0.10,0.14,0.09,0.12)
    for (m in c("DL","REML","PM")) print(rma(yi=y, sei=s, method=m)$tau2)

    DL   = 0.107946469966814
    REML = 0.114977616645157
    PM   = 0.115554997408346   # control=list(tol=1e-12); metafor's default PM
                               # tolerance leaves ~1e-6 slack, so the high-precision
                               # root is used here (this module's bisection matches
                               # it to ~3e-11).
    mu(DL) = 0.417151002843966

Tolerance 1e-6 per the project's R-validation rule.
"""

import numpy as np
import pytest

import tau2_estimators as te

# Heterogeneous fixture matched to the metafor reference above.
Y = np.array([0.20, 0.80, -0.10, 0.55, 0.35, 0.90, 0.05, 0.65])
S = np.array([0.10, 0.12, 0.11, 0.13, 0.10, 0.14, 0.09, 0.12])

REF = {"DL": 0.107946469966814, "REML": 0.114977616645157, "PM": 0.115554997408346}
REF_MU_DL = 0.417151002843966

TOL = 1e-6


@pytest.mark.parametrize("method,fn", [
    ("DL", te.dl_tau2), ("REML", te.reml_tau2), ("PM", te.pm_tau2)
])
def test_matches_metafor(method, fn):
    assert abs(fn(Y, S) - REF[method]) < TOL


def test_dl_reproduces_engine_moment_formula():
    # Independently recompute the DL moment estimator and compare.
    k = len(Y)
    w = 1.0 / S ** 2
    mu = np.sum(w * Y) / np.sum(w)
    Q = np.sum(w * (Y - mu) ** 2)
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    expected = max(0.0, (Q - (k - 1)) / c)
    assert abs(te.dl_tau2(Y, S) - expected) < 1e-12


def test_dl_pool_matches_metafor_mu():
    pool = te.pool_random_effects(Y, S, te.dl_tau2(Y, S), use_hksj=False, method="DL")
    assert abs(pool.estimate - REF_MU_DL) < TOL


def test_pm_solves_moment_equation():
    # At the PM root the weighted residual sum of squares equals k-1 exactly.
    t2 = te.pm_tau2(Y, S)
    k = len(Y)
    assert abs(te._pm_residual(t2, Y, S, k)) < 1e-6


def test_reml_and_pm_reduce_dl_bias_upward():
    # For this heterogeneous set REML and PM sit above DL (DL under-estimates).
    dl, reml, pm = te.dl_tau2(Y, S), te.reml_tau2(Y, S), te.pm_tau2(Y, S)
    assert dl < reml < pm  # strict for this fixture; all within 0.008 of each other


def test_homogeneous_data_gives_zero_tau2_all_methods():
    # Perfectly consistent estimates => Q <= k-1 => tau^2 = 0 for every estimator.
    y = np.array([0.30, 0.30, 0.30, 0.30, 0.30])
    s = np.array([0.10, 0.12, 0.09, 0.11, 0.10])
    for fn in (te.dl_tau2, te.reml_tau2, te.pm_tau2):
        assert fn(y, s) == 0.0


def test_low_heterogeneity_clamps_to_zero_not_negative():
    # Small scatter, large SEs => moment numerator negative => clamp, never < 0.
    rng = np.random.default_rng(7)
    y = 0.5 + rng.normal(0, 0.01, 6)
    s = np.full(6, 0.25)
    for fn in (te.dl_tau2, te.reml_tau2, te.pm_tau2):
        assert fn(y, s) >= 0.0


def test_k1_returns_zero_and_passthrough_pool():
    for fn in (te.dl_tau2, te.reml_tau2, te.pm_tau2):
        assert fn([0.4], [0.1]) == 0.0
    pool = te.pool_random_effects([0.4], [0.1], 0.0)
    assert pool.k == 1 and abs(pool.estimate - 0.4) < 1e-12


def test_all_tau2_nonnegative_random_stress():
    rng = np.random.default_rng(2026)
    for _ in range(200):
        k = int(rng.integers(3, 20))
        tau_true = float(rng.uniform(0, 0.3))
        se = rng.uniform(0.05, 0.25, k)
        y = rng.normal(0.2, np.sqrt(se ** 2 + tau_true), k)
        for fn in (te.dl_tau2, te.reml_tau2, te.pm_tau2):
            v = fn(y, se)
            assert np.isfinite(v) and v >= 0.0


def test_sensitivity_table_shape_and_default_first():
    out = te.tau2_sensitivity(Y, S)
    assert list(out.keys()) == ["DL", "REML", "PM"]  # DL (shipped default) first
    for name, pool in out.items():
        assert pool.method == name
        assert pool.ci_lower < pool.estimate < pool.ci_upper
        assert abs(pool.tau2 - REF[name]) < TOL


def test_input_validation():
    with pytest.raises(ValueError):
        te.dl_tau2([0.1, 0.2], [0.1])          # mismatched length
    with pytest.raises(ValueError):
        te.dl_tau2([0.1, 0.2], [0.1, -0.1])    # non-positive SE
    with pytest.raises(ValueError):
        te.dl_tau2([0.1, np.nan], [0.1, 0.2])  # non-finite
