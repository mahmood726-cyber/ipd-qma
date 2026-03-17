"""
IPD-QMA v2: Individual Participant Data Quantile Meta-Analysis
==============================================================
Corrected implementation addressing all P0 issues from multi-persona review:
  - Actual quantile regression (statsmodels.QuantReg) instead of np.percentile
  - DerSimonian-Laird random-effects pooling with tau^2 estimation
  - Correct lnVR SE formula (Nakagawa et al. 2015)
  - Configurable confidence level (no hardcoded z=1.96)
  - HKSJ correction option for small K
  - Seeded PRNG for reproducibility
  - Heterogeneity diagnostics (Q, I^2, tau^2, prediction interval)
  - Multivariate-aware slope test
  - Input validation and edge-case guards
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from dataclasses import dataclass
from typing import Optional
import warnings


@dataclass
class PooledResult:
    """Result of a random-effects meta-analysis pooling."""
    estimate: float
    se: float
    p_value: float
    ci_lower: float
    ci_upper: float
    tau2: float
    I2: float
    Q: float
    k: int
    pi_lower: float = np.nan  # prediction interval
    pi_upper: float = np.nan


@dataclass
class StudyResult:
    """Stage 1 results for a single study."""
    quantile_effects: np.ndarray
    se_quantiles: np.ndarray
    slope: float
    se_slope: float
    lnvr: float
    se_lnvr: float
    mean_diff: float
    se_mean_diff: float
    n_control: int
    n_treatment: int
    boot_effects: np.ndarray  # (n_quantiles, n_boot) for covariance
    study_label: str = ""


class IPDQMA:
    """
    IPD-QMA: Individual Participant Data Quantile Meta-Analysis.

    Two-Stage framework:
      Stage 1: Quantile regression at each tau to estimate unconditional QTEs
               (conditional and unconditional QTEs coincide in the two-sample
               no-covariate setting).
      Stage 2: DerSimonian-Laird (or REML) random-effects pooling across studies.

    Parameters
    ----------
    quantiles : list of float
        Quantile levels to estimate (default: [0.1, 0.25, 0.5, 0.75, 0.9]).
    n_boot : int
        Number of bootstrap resamples for SE estimation (default: 500).
    conf_level : float
        Confidence level for CIs (default: 0.95).
    use_hksj : bool
        Use Hartung-Knapp-Sidik-Jonkman correction (default: True for K>=3).
    seed : int or None
        Random seed for reproducibility.
    """

    def __init__(
        self,
        quantiles=None,
        n_boot: int = 500,
        conf_level: float = 0.95,
        use_hksj: bool = True,
        seed: Optional[int] = None,
    ):
        self.quantiles = quantiles if quantiles is not None else [0.1, 0.25, 0.5, 0.75, 0.9]
        self.n_boot = n_boot
        self.conf_level = conf_level
        self.use_hksj = use_hksj
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.results = None
        self._study_results = None

        # Validate quantiles contain 0.1 and 0.9 for slope test
        if 0.1 not in self.quantiles or 0.9 not in self.quantiles:
            warnings.warn(
                "Slope test requires quantiles 0.1 and 0.9. "
                "They will be added to the quantile grid."
            )
            if 0.1 not in self.quantiles:
                self.quantiles = [0.1] + self.quantiles
            if 0.9 not in self.quantiles:
                self.quantiles = self.quantiles + [0.9]
        self.quantiles = sorted(self.quantiles)

    def analyze_study(
        self,
        control: np.ndarray,
        treatment: np.ndarray,
        label: str = "",
    ) -> StudyResult:
        """
        Stage 1: Quantile regression + bootstrap SEs for a single study.

        Uses statsmodels.QuantReg to estimate unconditional quantile treatment
        effects (QTEs) at each specified quantile level.

        Parameters
        ----------
        control : array-like
            Outcome values for control arm.
        treatment : array-like
            Outcome values for treatment arm.
        label : str
            Optional study label.

        Returns
        -------
        StudyResult
        """
        control = np.asarray(control, dtype=np.float64)
        treatment = np.asarray(treatment, dtype=np.float64)

        # Remove NaN
        control = control[~np.isnan(control)]
        treatment = treatment[~np.isnan(treatment)]

        n_c, n_t = len(control), len(treatment)
        if n_c < 10 or n_t < 10:
            raise ValueError(
                f"Need >= 10 observations per arm for reliable quantile regression "
                f"(got n_c={n_c}, n_t={n_t})"
            )

        # Construct design matrix: Y = alpha + beta * Treatment
        y = np.concatenate([control, treatment])
        t_indicator = np.concatenate([np.zeros(n_c), np.ones(n_t)])
        X = sm.add_constant(t_indicator)

        # --- Point estimates via quantile regression ---
        n_q = len(self.quantiles)
        obs_effects = np.zeros(n_q)
        se_effects_asymp = np.zeros(n_q)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i, q in enumerate(self.quantiles):
                try:
                    model = sm.QuantReg(y, X)
                    res = model.fit(q=q, max_iter=1000, p_tol=1e-6)
                    obs_effects[i] = res.params[1]  # treatment coefficient
                    se_effects_asymp[i] = res.bse[1]
                except Exception as e:
                    warnings.warn(f"QuantReg failed at tau={q}: {e}. Using percentile fallback.")
                    obs_effects[i] = np.percentile(treatment, q * 100) - np.percentile(control, q * 100)
                    se_effects_asymp[i] = np.nan

        # --- Bootstrap SEs (preserves within-study cross-quantile covariance) ---
        boot_effects = np.zeros((n_q, self.n_boot))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for b in range(self.n_boot):
                # Resample within each arm (parallel-group design)
                bc = control[self.rng.integers(0, n_c, n_c)]
                bt = treatment[self.rng.integers(0, n_t, n_t)]
                yb = np.concatenate([bc, bt])
                Xb = sm.add_constant(np.concatenate([np.zeros(n_c), np.ones(n_t)]))

                for i, q in enumerate(self.quantiles):
                    try:
                        res_b = sm.QuantReg(yb, Xb).fit(q=q, max_iter=500, p_tol=1e-5)
                        boot_effects[i, b] = res_b.params[1]
                    except Exception:
                        # Fallback to percentile difference
                        boot_effects[i, b] = (
                            np.percentile(bt, q * 100) - np.percentile(bc, q * 100)
                        )

        se_effects = np.std(boot_effects, axis=1, ddof=1)

        # Use bootstrap SE (more robust at extreme quantiles)
        # Fall back to asymptotic SE if bootstrap is degenerate
        for i in range(n_q):
            if se_effects[i] < 1e-12 and np.isfinite(se_effects_asymp[i]):
                se_effects[i] = se_effects_asymp[i]

        # --- Slope (Q90 - Q10) computed within each bootstrap iteration ---
        idx_10 = self.quantiles.index(0.1)
        idx_90 = self.quantiles.index(0.9)
        obs_slope = obs_effects[idx_90] - obs_effects[idx_10]
        boot_slopes = boot_effects[idx_90, :] - boot_effects[idx_10, :]
        se_slope = np.std(boot_slopes, ddof=1)

        # --- Mean Difference (for comparison) ---
        mean_diff = np.mean(treatment) - np.mean(control)
        se_mean_diff = np.sqrt(np.var(control, ddof=1) / n_c + np.var(treatment, ddof=1) / n_t)

        # --- lnVR as log SD ratio (Nakagawa et al. 2015) ---
        # lnVR = log(SD_t / SD_c) = 0.5 * log(Var_t / Var_c)
        # Bias correction: + 1/(2*(n_t-1)) - 1/(2*(n_c-1))
        # SE = sqrt(1/(2*(n_t-1)) + 1/(2*(n_c-1)))
        var_t = np.var(treatment, ddof=1)
        var_c = np.var(control, ddof=1)
        if var_c > 0 and var_t > 0:
            lnvr = (0.5 * np.log(var_t / var_c)
                    + 1.0 / (2 * (n_t - 1)) - 1.0 / (2 * (n_c - 1)))
            se_lnvr = np.sqrt(1.0 / (2 * (n_t - 1)) + 1.0 / (2 * (n_c - 1)))
        else:
            lnvr = np.nan
            se_lnvr = np.nan

        return StudyResult(
            quantile_effects=obs_effects,
            se_quantiles=se_effects,
            slope=obs_slope,
            se_slope=se_slope,
            lnvr=lnvr,
            se_lnvr=se_lnvr,
            mean_diff=mean_diff,
            se_mean_diff=se_mean_diff,
            n_control=n_c,
            n_treatment=n_t,
            boot_effects=boot_effects,
            study_label=label,
        )

    def _pool_dl(self, est, se, use_hksj=None) -> PooledResult:
        """
        DerSimonian-Laird random-effects pooling.

        Parameters
        ----------
        est : array-like
            Study-level point estimates.
        se : array-like
            Study-level standard errors.
        use_hksj : bool or None
            Override HKSJ setting (default: use self.use_hksj).

        Returns
        -------
        PooledResult
        """
        est = np.asarray(est, dtype=np.float64)
        se = np.asarray(se, dtype=np.float64)
        k = len(est)

        if use_hksj is None:
            use_hksj = self.use_hksj

        if k < 1:
            raise ValueError("Need at least 1 study to pool")

        if k == 1:
            alpha = 1 - self.conf_level
            z_crit = stats.norm.ppf(1 - alpha / 2)
            se_safe = max(se[0], 1e-12)
            return PooledResult(
                estimate=est[0], se=se_safe,
                p_value=2 * (1 - stats.norm.cdf(abs(est[0] / se_safe))),
                ci_lower=est[0] - z_crit * se_safe,
                ci_upper=est[0] + z_crit * se_safe,
                tau2=0.0, I2=0.0, Q=0.0, k=1,
            )

        # Guard against zero SE
        se = np.maximum(se, 1e-12)

        # --- Fixed-effect quantities for DL estimation ---
        w_fe = 1.0 / se**2
        mu_fe = np.sum(est * w_fe) / np.sum(w_fe)

        # Cochran's Q
        Q = np.sum(w_fe * (est - mu_fe) ** 2)

        # DL tau^2 estimator
        c = np.sum(w_fe) - np.sum(w_fe**2) / np.sum(w_fe)
        tau2 = max(0.0, (Q - (k - 1)) / c)

        # I^2
        I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0

        # --- Random-effects weights ---
        w_re = 1.0 / (se**2 + tau2)
        p_est = np.sum(est * w_re) / np.sum(w_re)
        p_se_dl = np.sqrt(1.0 / np.sum(w_re))  # DL SE (used for PI)

        # --- Inference ---
        alpha = 1 - self.conf_level

        if use_hksj and k >= 3:
            # HKSJ correction: use t-distribution with k-1 df
            # and adjust SE by sqrt(q_hksj)
            q_hksj = np.sum(w_re * (est - p_est) ** 2) / (k - 1)
            p_se = p_se_dl * np.sqrt(max(1.0, q_hksj))
            t_crit = stats.t.ppf(1 - alpha / 2, df=k - 1)
            z_stat = p_est / p_se
            p_value = 2 * (1 - stats.t.cdf(abs(z_stat), df=k - 1))
            ci_lower = p_est - t_crit * p_se
            ci_upper = p_est + t_crit * p_se
        else:
            p_se = p_se_dl
            z_crit = stats.norm.ppf(1 - alpha / 2)
            z_stat = p_est / p_se
            p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))
            ci_lower = p_est - z_crit * p_se
            ci_upper = p_est + z_crit * p_se

        # Prediction interval (Riley et al. 2011) — uses DL SE, not HKSJ
        if k >= 3:
            t_pi = stats.t.ppf(1 - alpha / 2, df=k - 2)
            pi_se = np.sqrt(tau2 + p_se_dl**2)
            pi_lower = p_est - t_pi * pi_se
            pi_upper = p_est + t_pi * pi_se
        else:
            pi_lower = np.nan
            pi_upper = np.nan

        return PooledResult(
            estimate=p_est, se=p_se, p_value=p_value,
            ci_lower=ci_lower, ci_upper=ci_upper,
            tau2=tau2, I2=I2, Q=Q, k=k,
            pi_lower=pi_lower, pi_upper=pi_upper,
        )

    def fit(self, studies_data, labels=None):
        """
        Stage 2: Pool study-level QTEs using DerSimonian-Laird.

        Parameters
        ----------
        studies_data : list of (control_array, treatment_array) tuples
        labels : list of str, optional
            Study labels.

        Returns
        -------
        dict with keys 'profile', 'slope_test', 'lnvr_test', 'md_test',
             'study_results', 'heterogeneity'
        """
        if labels is None:
            labels = [f"Study {i+1}" for i in range(len(studies_data))]

        # Stage 1: Analyze each study
        print(f"Stage 1: Analyzing {len(studies_data)} studies...")
        self._study_results = []
        for idx, (c, t) in enumerate(studies_data):
            lbl = labels[idx] if idx < len(labels) else f"Study {idx+1}"
            print(f"  [{idx+1}/{len(studies_data)}] {lbl} (n_c={len(c)}, n_t={len(t)})")
            res = self.analyze_study(c, t, label=lbl)
            self._study_results.append(res)

        # Stage 2: Pool each quantile
        print(f"\nStage 2: Pooling with DerSimonian-Laird random effects...")
        n_q = len(self.quantiles)
        q_summary = []

        for i, q in enumerate(self.quantiles):
            est = [s.quantile_effects[i] for s in self._study_results]
            se = [s.se_quantiles[i] for s in self._study_results]
            pooled = self._pool_dl(est, se)
            q_summary.append({
                'Quantile': q,
                'Effect': pooled.estimate,
                'SE': pooled.se,
                'P': pooled.p_value,
                'CI_Lower': pooled.ci_lower,
                'CI_Upper': pooled.ci_upper,
                'tau2': pooled.tau2,
                'I2': pooled.I2,
                'Q': pooled.Q,
                'PI_Lower': pooled.pi_lower,
                'PI_Upper': pooled.pi_upper,
            })

        # Pool slope
        s_est = [s.slope for s in self._study_results]
        s_se = [s.se_slope for s in self._study_results]
        slope_pooled = self._pool_dl(s_est, s_se)

        # Pool lnVR
        l_est = [s.lnvr for s in self._study_results]
        l_se = [s.se_lnvr for s in self._study_results]
        # Remove NaN entries
        valid_lnvr = [(e, s) for e, s in zip(l_est, l_se) if np.isfinite(e) and np.isfinite(s)]
        if len(valid_lnvr) >= 1:
            lnvr_pooled = self._pool_dl(
                [v[0] for v in valid_lnvr],
                [v[1] for v in valid_lnvr],
            )
        else:
            lnvr_pooled = None

        # Pool mean difference (for comparison)
        md_est = [s.mean_diff for s in self._study_results]
        md_se = [s.se_mean_diff for s in self._study_results]
        md_pooled = self._pool_dl(md_est, md_se)

        profile_df = pd.DataFrame(q_summary)

        self.results = {
            'profile': profile_df,
            'slope_test': {
                'estimate': slope_pooled.estimate,
                'se': slope_pooled.se,
                'p_value': slope_pooled.p_value,
                'ci_lower': slope_pooled.ci_lower,
                'ci_upper': slope_pooled.ci_upper,
                'tau2': slope_pooled.tau2,
                'I2': slope_pooled.I2,
            },
            'lnvr_test': {
                'estimate': lnvr_pooled.estimate if lnvr_pooled else np.nan,
                'se': lnvr_pooled.se if lnvr_pooled else np.nan,
                'p_value': lnvr_pooled.p_value if lnvr_pooled else np.nan,
                'ci_lower': lnvr_pooled.ci_lower if lnvr_pooled else np.nan,
                'ci_upper': lnvr_pooled.ci_upper if lnvr_pooled else np.nan,
                'tau2': lnvr_pooled.tau2 if lnvr_pooled else np.nan,
                'I2': lnvr_pooled.I2 if lnvr_pooled else np.nan,
            },
            'md_test': {
                'estimate': md_pooled.estimate,
                'se': md_pooled.se,
                'p_value': md_pooled.p_value,
                'ci_lower': md_pooled.ci_lower,
                'ci_upper': md_pooled.ci_upper,
                'tau2': md_pooled.tau2,
                'I2': md_pooled.I2,
            },
            'study_results': self._study_results,
            'n_studies': len(studies_data),
            'conf_level': self.conf_level,
            'n_boot': self.n_boot,
            'seed': self.seed,
        }

        # Multivariate Wald test for flat QTE profile
        wald_result = self.wald_test()
        self.results['wald_test'] = wald_result

        self._print_summary()
        return self.results

    def _print_summary(self):
        """Print a formatted summary of results."""
        if self.results is None:
            return

        r = self.results
        cl_pct = int(self.conf_level * 100)
        print(f"\n{'='*70}")
        print(f"IPD-QMA RESULTS ({r['n_studies']} studies, {cl_pct}% CI, B={r['n_boot']})")
        print(f"{'='*70}")

        print(f"\n--- Quantile Treatment Effect Profile ---")
        profile = r['profile']
        print(f"{'Quantile':>10} {'Effect':>10} {'SE':>8} {'P-value':>10} "
              f"{'CI Lower':>10} {'CI Upper':>10} {'I2':>6} {'tau2':>8}")
        for _, row in profile.iterrows():
            sig = "*" if row['P'] < 0.05 else ""
            print(f"{row['Quantile']:>10.2f} {row['Effect']:>10.4f} {row['SE']:>8.4f} "
                  f"{row['P']:>10.4f}{sig:1s} {row['CI_Lower']:>10.4f} "
                  f"{row['CI_Upper']:>10.4f} {row['I2']:>5.1f}% {row['tau2']:>8.4f}")

        print(f"\n--- Slope Test (Q90 - Q10) ---")
        s = r['slope_test']
        sig = " ***" if s['p_value'] < 0.001 else (" **" if s['p_value'] < 0.01 else (" *" if s['p_value'] < 0.05 else ""))
        print(f"  Estimate: {s['estimate']:.4f} (SE: {s['se']:.4f})")
        print(f"  P-value:  {s['p_value']:.6f}{sig}")
        print(f"  CI:       [{s['ci_lower']:.4f}, {s['ci_upper']:.4f}]")
        print(f"  I2:       {s['I2']:.1f}%  |  tau2: {s['tau2']:.4f}")

        print(f"\n--- Wald Test (H0: flat QTE profile) ---")
        if 'wald_test' in r:
            wt = r['wald_test']
            if np.isfinite(wt['statistic']):
                sig = " ***" if wt['p_value'] < 0.001 else (" **" if wt['p_value'] < 0.01 else (" *" if wt['p_value'] < 0.05 else ""))
                print(f"  Chi-squared: {wt['statistic']:.4f}  (df = {wt['df']})")
                print(f"  P-value:     {wt['p_value']:.6f}{sig}")
                if wt['p_value'] < 0.05:
                    print(f"  => REJECT H0: treatment effects differ across quantiles.")
                else:
                    print(f"  => Cannot reject H0: no evidence of non-flat QTE profile.")
            else:
                print(f"  Not computable (singular covariance matrix)")

        print(f"\n--- Log SD Ratio (lnVR = log(SD_t/SD_c)) ---")
        lv = r['lnvr_test']
        if np.isfinite(lv['p_value']):
            sig = " ***" if lv['p_value'] < 0.001 else (" **" if lv['p_value'] < 0.01 else (" *" if lv['p_value'] < 0.05 else ""))
            print(f"  Estimate: {lv['estimate']:.4f} (SE: {lv['se']:.4f})")
            print(f"  P-value:  {lv['p_value']:.6f}{sig}")
        else:
            print(f"  Not computable (degenerate variance)")

        print(f"\n--- Standard Meta-Analysis (Mean Difference) ---")
        md = r['md_test']
        sig = " ***" if md['p_value'] < 0.001 else (" **" if md['p_value'] < 0.01 else (" *" if md['p_value'] < 0.05 else ""))
        print(f"  Estimate: {md['estimate']:.4f} (SE: {md['se']:.4f})")
        print(f"  P-value:  {md['p_value']:.6f}{sig}")
        print(f"  CI:       [{md['ci_lower']:.4f}, {md['ci_upper']:.4f}]")
        print(f"  I2:       {md['I2']:.1f}%  |  tau2: {md['tau2']:.4f}")

        # Interpretation
        print(f"\n--- Interpretation ---")
        if s['p_value'] < 0.05:
            print(f"  The Slope Test is SIGNIFICANT (p={s['p_value']:.4f}).")
            print(f"  Treatment effects vary across the outcome distribution.")
            print(f"  Standard meta-analysis of mean differences may be misleading.")
        else:
            print(f"  The Slope Test is NOT significant (p={s['p_value']:.4f}).")
            print(f"  No evidence of heterogeneous treatment effects across quantiles.")
        print(f"{'='*70}\n")

    def wald_test(self):
        """
        Multivariate Wald test for H0: all quantile treatment effects are equal.

        Tests whether the QTE profile is flat (i.e., treatment effect does not
        vary across quantiles) using the full bootstrap covariance matrix.

        The test constructs adjacent-difference contrasts C (a (q-1) x q matrix)
        and computes W = (C @ theta)' @ inv(C @ V @ C') @ (C @ theta), which
        follows chi2(q-1) under H0.

        Returns
        -------
        dict with keys:
            'statistic' : float — Wald chi-squared statistic
            'df' : int — degrees of freedom (q - 1)
            'p_value' : float — p-value from chi2 distribution
            'contrast_matrix' : np.ndarray — the (q-1) x q contrast matrix C

        Raises
        ------
        ValueError
            If .fit() has not been called.
        """
        if self._study_results is None:
            raise ValueError("Must call .fit() before wald_test()")

        n_q = len(self.quantiles)
        k = len(self._study_results)

        if n_q < 2:
            raise ValueError("Need at least 2 quantiles for Wald test")

        # --- Construct pooled QTE vector theta_hat ---
        # Get tau2 for each quantile from the profile results
        profile = self.results['profile']
        theta_hat = np.array(profile['Effect'].values, dtype=np.float64)
        tau2_vec = np.array(profile['tau2'].values, dtype=np.float64)

        # --- Compute RE weights for each study-quantile pair ---
        # w_ki = 1 / (SE_ki^2 + tau2_i)
        W = np.zeros((k, n_q))
        for s_idx, study in enumerate(self._study_results):
            for q_idx in range(n_q):
                se_ki = max(study.se_quantiles[q_idx], 1e-12)
                W[s_idx, q_idx] = 1.0 / (se_ki ** 2 + tau2_vec[q_idx])

        # Sum of weights per quantile
        W_sum = W.sum(axis=0)  # shape (n_q,)

        # --- Compute pooled covariance matrix V_pooled ---
        # V_pooled[i, j] = sum_k (w_ki * w_kj * Sigma_k[i, j]) / (W_sum[i] * W_sum[j])
        V_pooled = np.zeros((n_q, n_q))

        for s_idx, study in enumerate(self._study_results):
            # Study-level covariance from bootstrap: shape (n_q, n_q)
            if study.boot_effects.shape[1] < 2:
                warnings.warn(
                    f"Study '{study.study_label}' has < 2 bootstrap samples; "
                    "skipping in Wald covariance."
                )
                continue
            Sigma_k = np.cov(study.boot_effects)  # (n_q, n_q)
            # Handle case where np.cov returns scalar for n_q == 1
            Sigma_k = np.atleast_2d(Sigma_k)

            for i in range(n_q):
                for j in range(n_q):
                    V_pooled[i, j] += W[s_idx, i] * W[s_idx, j] * Sigma_k[i, j]

        # Normalize by product of weight sums
        for i in range(n_q):
            for j in range(n_q):
                denom = W_sum[i] * W_sum[j]
                if denom > 0:
                    V_pooled[i, j] /= denom
                else:
                    V_pooled[i, j] = 0.0

        # --- Contrast matrix C: adjacent differences, (q-1) x q ---
        C = np.zeros((n_q - 1, n_q))
        for j in range(n_q - 1):
            C[j, j] = 1.0
            C[j, j + 1] = -1.0

        # --- Wald statistic: W = (C @ theta)' @ inv(C @ V @ C') @ (C @ theta) ---
        C_theta = C @ theta_hat  # shape (q-1,)
        C_V_Ct = C @ V_pooled @ C.T  # shape (q-1, q-1)

        try:
            # Use solve instead of inv for numerical stability:
            # W = C_theta' @ inv(C_V_Ct) @ C_theta
            # = C_theta' @ x, where C_V_Ct @ x = C_theta
            x = np.linalg.solve(C_V_Ct, C_theta)
            W_stat = float(C_theta @ x)
        except np.linalg.LinAlgError:
            warnings.warn(
                "Wald test: contrast covariance matrix is singular. "
                "Cannot compute Wald statistic."
            )
            return {
                'statistic': np.nan,
                'df': n_q - 1,
                'p_value': np.nan,
                'contrast_matrix': C,
            }

        df = n_q - 1
        p_value = 1.0 - stats.chi2.cdf(W_stat, df)

        return {
            'statistic': W_stat,
            'df': df,
            'p_value': p_value,
            'contrast_matrix': C,
        }

    def plot_profile(self, ax=None, show=True, save_path=None):
        """
        Plot the Quantile Treatment Effect Profile with CI band.

        Parameters
        ----------
        ax : matplotlib Axes, optional
        show : bool
        save_path : str, optional
            Path to save the figure.

        Returns
        -------
        matplotlib Axes
        """
        if self.results is None:
            raise ValueError("Must call .fit() before plotting")

        df = self.results['profile']

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        else:
            fig = ax.figure

        # Zero line
        ax.axhline(0, color='gray', linestyle='--', alpha=0.7, linewidth=1)

        # Prediction interval (if available)
        if 'PI_Lower' in df.columns and not df['PI_Lower'].isna().all():
            ax.fill_between(
                df['Quantile'], df['PI_Lower'], df['PI_Upper'],
                color='#AED6F1', alpha=0.3, label=f'Prediction Interval'
            )

        # Confidence interval
        ax.fill_between(
            df['Quantile'], df['CI_Lower'], df['CI_Upper'],
            color='#2E86C1', alpha=0.3,
            label=f'{int(self.conf_level*100)}% CI'
        )

        # Point estimates
        ax.plot(
            df['Quantile'], df['Effect'], 'o-',
            color='#2E86C1', linewidth=2.5, markersize=8,
            label='Pooled QTE', zorder=5
        )

        # Mark significant quantiles
        for _, row in df.iterrows():
            if row['P'] < 0.05:
                ax.plot(row['Quantile'], row['Effect'], 'o',
                        color='#E74C3C', markersize=10, zorder=6, markeredgewidth=2,
                        markerfacecolor='none')

        ax.set_xlabel("Outcome Quantile (tau)", fontsize=12)
        ax.set_ylabel("Treatment Effect (QTE)", fontsize=12)
        ax.set_title("IPD-QMA: Quantile Treatment Effect Profile", fontsize=14)
        ax.legend(loc='upper left', fontsize=10)
        ax.grid(True, alpha=0.3)

        # Add slope test annotation
        s = self.results['slope_test']
        p_str = "<0.001" if s['p_value'] < 0.001 else f"{s['p_value']:.3f}"
        sig_text = f"Slope Test: {s['estimate']:.3f} (p={p_str})"
        ax.annotate(
            sig_text,
            xy=(0.98, 0.02), xycoords='axes fraction',
            ha='right', va='bottom', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='gray')
        )

        if ax is None or ax.figure is not None:
            try:
                ax.figure.tight_layout()
            except Exception:
                pass
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Figure saved to {save_path}")
        if show:
            plt.show()
        return ax

    def plot_forest(self, quantile_idx=None, ax=None, show=True, save_path=None):
        """
        Forest plot for a specific quantile across studies.

        Parameters
        ----------
        quantile_idx : int, optional
            Index into self.quantiles (default: median quantile).
        """
        if self.results is None or self._study_results is None:
            raise ValueError("Must call .fit() before plotting")

        if quantile_idx is None:
            # Default to median
            quantile_idx = len(self.quantiles) // 2

        q = self.quantiles[quantile_idx]
        studies = self._study_results
        profile_row = self.results['profile'].iloc[quantile_idx]

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, max(4, len(studies) * 0.5 + 2)))

        y_positions = list(range(len(studies), 0, -1))

        for i, (s, yp) in enumerate(zip(studies, y_positions)):
            est = s.quantile_effects[quantile_idx]
            se = s.se_quantiles[quantile_idx]
            alpha = 1 - self.conf_level
            z_crit = stats.norm.ppf(1 - alpha / 2)
            ci_lo = est - z_crit * se
            ci_hi = est + z_crit * se

            ax.plot([ci_lo, ci_hi], [yp, yp], 'b-', linewidth=1.5)
            ax.plot(est, yp, 'bs', markersize=6)
            ax.text(-0.02, yp, s.study_label, ha='right', va='center',
                    transform=ax.get_yaxis_transform(), fontsize=9)

        # Pooled diamond
        p_est = profile_row['Effect']
        p_lo = profile_row['CI_Lower']
        p_hi = profile_row['CI_Upper']
        yp = 0
        diamond_x = [p_lo, p_est, p_hi, p_est]
        diamond_y = [yp, yp + 0.3, yp, yp - 0.3]
        ax.fill(diamond_x, diamond_y, color='red', alpha=0.6)
        ax.text(-0.02, yp, "Pooled", ha='right', va='center',
                transform=ax.get_yaxis_transform(), fontsize=9, fontweight='bold')

        ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel(f"Treatment Effect at Q{int(q*100)}", fontsize=12)
        ax.set_title(f"Forest Plot: Quantile {q:.2f}", fontsize=13)
        ax.set_yticks([])
        ax.grid(True, axis='x', alpha=0.3)

        plt.tight_layout()
        if save_path:
            ax.figure.savefig(save_path, dpi=150, bbox_inches='tight')
        if show:
            plt.show()
        return ax

    def plot_comparison(self, show=True, save_path=None):
        """
        Side-by-side comparison: MD vs QTE profile vs lnVR.
        """
        if self.results is None:
            raise ValueError("Must call .fit() before plotting")

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        # Panel 1: Standard MA (Mean Difference)
        ax1 = axes[0]
        md = self.results['md_test']
        ax1.bar(['Mean\nDifference'], [md['estimate']],
                yerr=[[md['estimate'] - md['ci_lower']], [md['ci_upper'] - md['estimate']]],
                color='#3498DB', capsize=10, width=0.4)
        ax1.axhline(0, color='gray', linestyle='--')
        p_str = "<0.001" if md['p_value'] < 0.001 else f"{md['p_value']:.3f}"
        ax1.set_title(f"Standard MA\np={p_str}", fontsize=12)
        ax1.set_ylabel("Effect Size")

        # Panel 2: QTE Profile
        ax2 = axes[1]
        self.plot_profile(ax=ax2, show=False)
        ax2.set_title("IPD-QMA Profile", fontsize=12)

        # Panel 3: lnVR
        ax3 = axes[2]
        lv = self.results['lnvr_test']
        if np.isfinite(lv['estimate']):
            lv_lo = lv['ci_lower']
            lv_hi = lv['ci_upper']
            ax3.bar(['lnVR'], [lv['estimate']],
                    yerr=[[lv['estimate'] - lv_lo], [lv_hi - lv['estimate']]],
                    color='#E67E22', capsize=10, width=0.4)
            p_str = "<0.001" if lv['p_value'] < 0.001 else f"{lv['p_value']:.3f}"
            ax3.set_title(f"Log SD Ratio\np={p_str}", fontsize=12)
        else:
            ax3.set_title("lnVR\n(not computable)")
        ax3.axhline(0, color='gray', linestyle='--')
        ax3.set_ylabel("lnVR (log SD ratio)")

        fig.suptitle("IPD-QMA vs Standard Meta-Analysis", fontsize=14, y=1.02)
        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
        if show:
            plt.show()
        return fig


def simulate_location_scale(
    K: int = 20,
    n_range: tuple = (100, 300),
    variance_ratio: float = 3.0,
    mean_shift: float = 0.0,
    distribution: str = "normal",
    seed: int = 42,
):
    """
    Generate simulated IPD for a location-scale scenario.

    Parameters
    ----------
    K : int
        Number of studies.
    n_range : tuple
        (min, max) sample size per arm per study.
    variance_ratio : float
        Treatment variance / Control variance.
    mean_shift : float
        True mean difference (0 for pure scale shift).
    distribution : str
        "normal", "exponential", "lognormal", or "t5" (t with 5 df).
    seed : int
        Random seed.

    Returns
    -------
    studies_data : list of (control, treatment) arrays
    labels : list of study labels
    true_params : dict of ground-truth values
    """
    rng = np.random.default_rng(seed)
    studies_data = []
    labels = []

    scale_factor = np.sqrt(variance_ratio)

    for k in range(K):
        n = rng.integers(n_range[0], n_range[1] + 1)
        n_c = n
        n_t = n

        if distribution == "normal":
            control = rng.standard_normal(n_c)
            treatment = rng.standard_normal(n_t) * scale_factor + mean_shift
        elif distribution == "exponential":
            # Center by POPULATION mean (Exp(s) has mean=s) to preserve
            # sampling variability in the mean difference
            control = rng.exponential(1.0, n_c) - 1.0
            treatment = rng.exponential(scale_factor, n_t) - scale_factor + mean_shift
        elif distribution == "lognormal":
            # Center by POPULATION mean to preserve sampling variability
            sigma_ln = 0.5
            pop_mean_ln = np.exp(sigma_ln**2 / 2)  # E[LogN(0, sigma)]
            control = rng.lognormal(0, sigma_ln, n_c) - pop_mean_ln
            treatment = rng.lognormal(0, sigma_ln, n_t) * scale_factor - scale_factor * pop_mean_ln + mean_shift
        elif distribution == "t5":
            control = rng.standard_t(5, n_c)
            treatment = rng.standard_t(5, n_t) * scale_factor + mean_shift
        else:
            raise ValueError(f"Unknown distribution: {distribution}")

        studies_data.append((control, treatment))
        labels.append(f"Study {k+1} (n={n_c}+{n_t})")

    # Compute true quantile differences (analytical ground truth)
    true_q_diffs = {}
    for q in [0.1, 0.25, 0.5, 0.75, 0.9]:
        if distribution == "normal":
            z_q = stats.norm.ppf(q)
            true_diff = z_q * (scale_factor - 1) + mean_shift
        elif distribution == "exponential":
            # Exp(s) quantile = -s*ln(1-q), mean = s.
            # After population centering: Q(p) - s = s*(-ln(1-q) - 1)
            # True QTE = (sf-1)*(-ln(1-q) - 1)
            q_exp = -np.log(1 - q)
            true_diff = (scale_factor - 1) * (q_exp - 1.0) + mean_shift
        elif distribution == "lognormal":
            # Control: LogN(0, 0.5) centered; Treatment: LogN(0, 0.5) * sf centered
            # Q_LN(q) = exp(mu + sigma * z_q); after centering: Q - E[X]
            # True QTE = sf * (Q_LN(q) - E[LN]) - (Q_LN(q) - E[LN])
            #          = (sf - 1) * (Q_LN(q) - E[LN])
            sigma_ln = 0.5
            q_ln = np.exp(sigma_ln * stats.norm.ppf(q))  # mu=0
            e_ln = np.exp(sigma_ln**2 / 2)  # E[LogN(0, sigma)]
            true_diff = (scale_factor - 1) * (q_ln - e_ln) + mean_shift
        elif distribution == "t5":
            t_q = stats.t.ppf(q, df=5)
            true_diff = t_q * (scale_factor - 1) + mean_shift
        else:
            true_diff = None
        true_q_diffs[q] = true_diff

    true_params = {
        'mean_diff': mean_shift,
        'variance_ratio': variance_ratio,
        'quantile_diffs': true_q_diffs,
        'distribution': distribution,
        'K': K,
        'n_range': n_range,
    }

    return studies_data, labels, true_params
