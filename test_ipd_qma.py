"""
Comprehensive pytest test suite for IPD-QMA.

Tests cover:
  - DerSimonian-Laird pooling (hand-calculated, homogeneous, single-study)
  - HKSJ correction (wider CI than DL)
  - Stage 1 analyze_study (normal location shift, scale shift)
  - Ground truth formulas (normal, t5)
  - Full pipeline (fit + slope test + lnVR test)
  - Prediction intervals (K>=3 vs K=2)
  - Coverage probability (null simulation sanity check)
"""

import numpy as np
import pytest
from scipy import stats

from ipd_qma import IPDQMA, PooledResult


# ============================================================
# Helper: standalone DL pooling (no HKSJ) for coverage test
# ============================================================

def _pool_dl_simple(est, se, conf_level=0.95):
    """
    Minimal DL pooling matching run_simulation_fast.pool_dl (HKSJ variant).
    Returns (estimate, se, p_value, ci_lower, ci_upper).
    """
    est = np.asarray(est, dtype=float)
    se = np.asarray(se, dtype=float)
    k = len(est)
    se = np.maximum(se, 1e-12)
    alpha = 1 - conf_level

    if k == 1:
        z_crit = stats.norm.ppf(1 - alpha / 2)
        z = est[0] / se[0]
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        return est[0], se[0], p, est[0] - z_crit * se[0], est[0] + z_crit * se[0]

    w_fe = 1.0 / se**2
    mu_fe = np.sum(est * w_fe) / np.sum(w_fe)
    Q = np.sum(w_fe * (est - mu_fe) ** 2)
    c = np.sum(w_fe) - np.sum(w_fe**2) / np.sum(w_fe)
    tau2 = max(0.0, (Q - (k - 1)) / c)

    w_re = 1.0 / (se**2 + tau2)
    p_est = np.sum(est * w_re) / np.sum(w_re)
    p_se = np.sqrt(1.0 / np.sum(w_re))

    if k >= 3:
        q_hksj = np.sum(w_re * (est - p_est) ** 2) / (k - 1)
        p_se_adj = p_se * np.sqrt(max(1.0, q_hksj))
        t_crit = stats.t.ppf(1 - alpha / 2, df=k - 1)
        z = p_est / p_se_adj
        p_value = 2 * (1 - stats.t.cdf(abs(z), df=k - 1))
        ci_lo = p_est - t_crit * p_se_adj
        ci_hi = p_est + t_crit * p_se_adj
        return p_est, p_se_adj, p_value, ci_lo, ci_hi
    else:
        z_crit = stats.norm.ppf(1 - alpha / 2)
        z = p_est / p_se
        p_value = 2 * (1 - stats.norm.cdf(abs(z)))
        ci_lo = p_est - z_crit * p_se
        ci_hi = p_est + z_crit * p_se
        return p_est, p_se, p_value, ci_lo, ci_hi


# ============================================================
# 1. DL pooling — hand-calculated known values
# ============================================================

class TestPoolDL:

    def test_pool_dl_known(self):
        """Pool 3 studies with hand-computed DL result."""
        est = np.array([0.5, 0.3, 0.7])
        se = np.array([0.1, 0.15, 0.12])

        # Hand calculation:
        w_fe = 1.0 / se**2  # [100, 44.444, 69.444]
        sum_w = np.sum(w_fe)  # 213.889
        mu_fe = np.sum(est * w_fe) / sum_w
        Q = np.sum(w_fe * (est - mu_fe) ** 2)
        c = sum_w - np.sum(w_fe**2) / sum_w
        tau2 = max(0.0, (Q - 2) / c)

        w_re = 1.0 / (se**2 + tau2)
        expected_est = np.sum(est * w_re) / np.sum(w_re)
        expected_se_dl = np.sqrt(1.0 / np.sum(w_re))

        # Use IPDQMA._pool_dl with HKSJ disabled for clean comparison
        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], seed=42)
        result = qma._pool_dl(est, se, use_hksj=False)

        assert isinstance(result, PooledResult)
        assert result.k == 3
        assert np.isclose(result.estimate, expected_est, atol=1e-10)
        assert np.isclose(result.se, expected_se_dl, atol=1e-10)
        assert np.isclose(result.tau2, tau2, atol=1e-10)
        assert result.p_value >= 0
        assert result.p_value <= 1
        assert result.ci_lower < result.estimate < result.ci_upper

    def test_pool_dl_homogeneous(self):
        """Three identical estimates -> tau2=0, estimate=1.0."""
        est = np.array([1.0, 1.0, 1.0])
        se = np.array([0.1, 0.1, 0.1])

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], seed=42)
        result = qma._pool_dl(est, se, use_hksj=False)

        assert np.isclose(result.tau2, 0.0, atol=1e-12)
        assert np.isclose(result.estimate, 1.0, atol=1e-10)
        assert np.isclose(result.I2, 0.0, atol=1e-10)
        assert np.isclose(result.Q, 0.0, atol=1e-10)

    def test_pool_dl_single_study(self):
        """K=1 -> returns study estimate and SE directly."""
        est = np.array([2.5])
        se = np.array([0.3])

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], seed=42)
        result = qma._pool_dl(est, se)

        assert result.k == 1
        assert np.isclose(result.estimate, 2.5, atol=1e-10)
        assert np.isclose(result.se, 0.3, atol=1e-10)
        assert np.isclose(result.tau2, 0.0, atol=1e-10)
        assert np.isclose(result.I2, 0.0, atol=1e-10)


# ============================================================
# 2. HKSJ correction
# ============================================================

class TestHKSJ:

    def test_hksj_wider_ci(self):
        """HKSJ CI should be >= DL CI width on heterogeneous data."""
        est = np.array([0.1, 0.5, 0.3, 0.9, 0.2])
        se = np.array([0.1, 0.15, 0.12, 0.08, 0.2])

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], seed=42)

        result_dl = qma._pool_dl(est, se, use_hksj=False)
        result_hksj = qma._pool_dl(est, se, use_hksj=True)

        width_dl = result_dl.ci_upper - result_dl.ci_lower
        width_hksj = result_hksj.ci_upper - result_hksj.ci_lower

        # HKSJ should produce wider (or equal) CI
        assert width_hksj >= width_dl - 1e-12, (
            f"HKSJ CI width ({width_hksj:.6f}) should be >= DL CI width ({width_dl:.6f})"
        )

        # Both should give same point estimate (same weights)
        assert np.isclose(result_dl.estimate, result_hksj.estimate, atol=1e-10)


# ============================================================
# 3. Stage 1: analyze_study
# ============================================================

class TestAnalyzeStudy:

    def test_analyze_study_normal(self):
        """N(0,1) vs N(0.5,1): median QTE ~ 0.5, slope ~ 0, lnVR ~ 0."""
        rng = np.random.default_rng(42)
        n = 200
        control = rng.standard_normal(n)
        treatment = rng.standard_normal(n) + 0.5

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], n_boot=200, seed=42)
        res = qma.analyze_study(control, treatment, label="test")

        # Median QTE should be near 0.5 (pure location shift)
        idx_50 = qma.quantiles.index(0.5)
        assert abs(res.quantile_effects[idx_50] - 0.5) < 0.2, (
            f"Median QTE = {res.quantile_effects[idx_50]:.4f}, expected ~0.5"
        )

        # Slope (Q90 - Q10) should be near 0 for pure location shift
        # Tolerance 0.45: quantile regression on n=200 has non-trivial
        # sampling variability at extreme quantiles (0.1, 0.9), so slopes
        # up to ~0.4 can occur by chance. A true scale shift (VR=2+)
        # produces slopes > 0.8.
        assert abs(res.slope) < 0.45, (
            f"Slope = {res.slope:.4f}, expected ~0 for location shift"
        )

        # lnVR should be near 0 (equal variances)
        assert abs(res.lnvr) < 0.15, (
            f"lnVR = {res.lnvr:.4f}, expected ~0 for equal variances"
        )

        # Mean difference should be near 0.5
        assert abs(res.mean_diff - 0.5) < 0.2, (
            f"Mean diff = {res.mean_diff:.4f}, expected ~0.5"
        )

    def test_analyze_study_scale_shift(self):
        """N(0,1) vs N(0,2): Q90 effect positive, Q10 negative, slope positive."""
        rng = np.random.default_rng(42)
        n = 200
        control = rng.standard_normal(n)
        treatment = rng.standard_normal(n) * 2.0  # scale factor = 2

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], n_boot=200, seed=42)
        res = qma.analyze_study(control, treatment, label="scale_test")

        idx_10 = qma.quantiles.index(0.1)
        idx_90 = qma.quantiles.index(0.9)

        # Q90 effect should be positive (treatment upper tail stretched)
        assert res.quantile_effects[idx_90] > 0.5, (
            f"Q90 effect = {res.quantile_effects[idx_90]:.4f}, expected > 0.5"
        )

        # Q10 effect should be negative (treatment lower tail stretched)
        assert res.quantile_effects[idx_10] < -0.5, (
            f"Q10 effect = {res.quantile_effects[idx_10]:.4f}, expected < -0.5"
        )

        # Slope should be strongly positive
        assert res.slope > 0.8, (
            f"Slope = {res.slope:.4f}, expected > 0.8"
        )


# ============================================================
# 4. Ground truth formulas (inline, matching compute_ground_truth)
# ============================================================

class TestGroundTruth:

    def test_ground_truth_normal(self):
        """
        Normal, VR=2.0, quantiles=[0.1, 0.5, 0.9], mean_shift=0.
        QTE at Q50 = 0 exactly.
        QTE at Q90 = z_0.9 * (sqrt(2) - 1) = 1.2816 * 0.4142 = 0.5310.
        Slope = Q90 - Q10 = 2 * 0.5310 = 1.0621.
        lnVR = 0.5 * log(2) = 0.3466.
        """
        vr = 2.0
        sf = np.sqrt(vr)
        quantiles = [0.1, 0.5, 0.9]
        mean_shift = 0.0

        # Compute QTE for each quantile
        qte = []
        for tau in quantiles:
            z_tau = stats.norm.ppf(tau)
            qte.append(z_tau * (sf - 1) + mean_shift)
        qte = np.array(qte)

        # Q50 = 0 exactly (z_0.5 = 0)
        assert np.isclose(qte[1], 0.0, atol=1e-10), f"Q50 QTE = {qte[1]}, expected 0"

        # Q90 = z_0.9 * (sqrt(2) - 1)
        expected_q90 = stats.norm.ppf(0.9) * (np.sqrt(2) - 1)
        assert np.isclose(qte[2], expected_q90, atol=1e-3), (
            f"Q90 QTE = {qte[2]:.6f}, expected {expected_q90:.6f}"
        )
        assert np.isclose(expected_q90, 0.5310, atol=1e-3), (
            f"Q90 expected = {expected_q90:.4f}, should be ~0.5310"
        )

        # Slope = Q90 - Q10
        slope = qte[2] - qte[0]
        expected_slope = 2 * expected_q90  # symmetric for normal
        assert np.isclose(slope, expected_slope, atol=1e-3), (
            f"Slope = {slope:.6f}, expected {expected_slope:.6f}"
        )
        assert np.isclose(expected_slope, 1.0621, atol=1e-3), (
            f"Slope expected = {expected_slope:.4f}, should be ~1.0621"
        )

        # lnVR = 0.5 * log(VR)
        lnvr = 0.5 * np.log(vr)
        assert np.isclose(lnvr, 0.3466, atol=1e-3), (
            f"lnVR = {lnvr:.4f}, expected ~0.3466"
        )

    def test_ground_truth_t5(self):
        """
        t5 distribution, VR=3.0, quantiles=[0.1, 0.5, 0.9].
        Q90 = t5_inv(0.9) * (sqrt(3) - 1) = 1.4759 * 0.7321 = 1.0805.
        """
        vr = 3.0
        sf = np.sqrt(vr)
        quantiles = [0.1, 0.5, 0.9]

        t5_q90 = stats.t.ppf(0.9, df=5)
        expected_q90 = t5_q90 * (sf - 1)

        assert np.isclose(t5_q90, 1.4759, atol=1e-3), (
            f"t5_inv(0.9) = {t5_q90:.4f}, expected ~1.4759"
        )
        assert np.isclose(expected_q90, 1.0805, atol=1e-3), (
            f"Q90 QTE = {expected_q90:.4f}, expected ~1.0805"
        )

        # Also verify Q50 = 0 (t distribution is symmetric)
        t5_q50 = stats.t.ppf(0.5, df=5)
        assert np.isclose(t5_q50, 0.0, atol=1e-10)


# ============================================================
# 5. Full pipeline (IPDQMA.fit)
# ============================================================

class TestFullPipeline:

    def test_full_pipeline_small(self):
        """
        5 studies with known scale shift (VR=3). Slope and lnVR tests
        should be significant.
        """
        rng = np.random.default_rng(42)
        vr = 3.0
        sf = np.sqrt(vr)

        studies_data = []
        for _ in range(5):
            n = rng.integers(80, 200)
            c = rng.standard_normal(n)
            t = rng.standard_normal(n) * sf
            studies_data.append((c, t))

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], n_boot=50, seed=42)
        results = qma.fit(studies_data)

        # Slope test should detect scale shift
        assert results['slope_test']['p_value'] < 0.05, (
            f"Slope p-value = {results['slope_test']['p_value']:.4f}, expected < 0.05"
        )

        # lnVR test should detect variance ratio
        assert results['lnvr_test']['p_value'] < 0.05, (
            f"lnVR p-value = {results['lnvr_test']['p_value']:.4f}, expected < 0.05"
        )

        # Profile should have correct structure
        profile = results['profile']
        assert len(profile) == 3  # 3 quantiles
        assert list(profile.columns) == [
            'Quantile', 'Effect', 'SE', 'P', 'CI_Lower', 'CI_Upper',
            'tau2', 'I2', 'Q', 'PI_Lower', 'PI_Upper'
        ]

        # n_studies should match
        assert results['n_studies'] == 5


# ============================================================
# 6. Prediction intervals
# ============================================================

class TestPredictionInterval:

    def test_prediction_interval_exists(self):
        """After .fit() with K>=3, PI bounds should not be NaN."""
        rng = np.random.default_rng(99)
        studies_data = []
        for _ in range(5):
            n = rng.integers(50, 100)
            c = rng.standard_normal(n)
            t = rng.standard_normal(n) * 1.5
            studies_data.append((c, t))

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], n_boot=50, seed=99)
        results = qma.fit(studies_data)

        profile = results['profile']
        # At least one quantile should have non-NaN PI
        assert not profile['PI_Lower'].isna().all(), "All PI_Lower are NaN for K=5"
        assert not profile['PI_Upper'].isna().all(), "All PI_Upper are NaN for K=5"

        # PI should be wider than CI
        for _, row in profile.iterrows():
            if np.isfinite(row['PI_Lower']) and np.isfinite(row['PI_Upper']):
                pi_width = row['PI_Upper'] - row['PI_Lower']
                ci_width = row['CI_Upper'] - row['CI_Lower']
                assert pi_width >= ci_width - 1e-10, (
                    f"PI width ({pi_width:.4f}) should be >= CI width ({ci_width:.4f})"
                )

    def test_k2_no_prediction_interval(self):
        """With K=2, PI bounds should be NaN."""
        rng = np.random.default_rng(77)
        studies_data = []
        for _ in range(2):
            n = rng.integers(50, 100)
            c = rng.standard_normal(n)
            t = rng.standard_normal(n) * 2.0
            studies_data.append((c, t))

        qma = IPDQMA(quantiles=[0.1, 0.5, 0.9], n_boot=50, seed=77)
        results = qma.fit(studies_data)

        profile = results['profile']
        assert profile['PI_Lower'].isna().all(), (
            f"PI_Lower should be NaN for K=2, got {profile['PI_Lower'].values}"
        )
        assert profile['PI_Upper'].isna().all(), (
            f"PI_Upper should be NaN for K=2, got {profile['PI_Upper'].values}"
        )


# ============================================================
# 7. Coverage probability (null simulation sanity check)
# ============================================================

class TestCoverage:

    def test_coverage_probability_null(self):
        """
        Quick check: 100 null simulations (VR=1, Normal, K=5, n=50-100).
        Pool MD with DL+HKSJ, compute coverage of true MD=0.
        Coverage should be >= 0.85 (loose bound for 100 sims).
        """
        n_sims = 100
        quantiles = [0.1, 0.5, 0.9]
        n_boot = 50
        K = 5
        true_md = 0.0
        covered = 0

        for sim in range(n_sims):
            rng_data = np.random.default_rng(10000 + sim)
            rng_boot = np.random.default_rng(10000 + sim + 1_000_000)

            # Generate K null studies (VR=1, no mean shift)
            md_ests = []
            md_ses = []
            for _ in range(K):
                n = rng_data.integers(50, 101)
                c = rng_data.standard_normal(n)
                t = rng_data.standard_normal(n)  # VR=1, shift=0

                md = np.mean(t) - np.mean(c)
                se_md = np.sqrt(np.var(c, ddof=1) / n + np.var(t, ddof=1) / n)
                md_ests.append(md)
                md_ses.append(se_md)

            # Pool with DL + HKSJ
            _, _, _, ci_lo, ci_hi = _pool_dl_simple(md_ests, md_ses, conf_level=0.95)

            if ci_lo <= true_md <= ci_hi:
                covered += 1

        coverage = covered / n_sims
        assert coverage >= 0.85, (
            f"Coverage = {coverage:.2f} ({covered}/{n_sims}), expected >= 0.85"
        )
