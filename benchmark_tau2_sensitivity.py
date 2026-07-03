#!/usr/bin/env python
"""
Worked example + reproducible benchmark: tau^2 estimator sensitivity in IPD-QMA.
================================================================================

The IPD-QMA engine pools each quantile's study-level treatment effects with the
DerSimonian-Laird (DL) moment estimator. The pre-registered protocol names REML
as the intended primary estimator and calls for a DL / REML / Paule-Mandel (PM)
sensitivity analysis (manuscript Limitations / Future directions). This script is
the runnable realisation of that analysis and doubles as the package's worked
example.

It is FULLY SELF-CONTAINED and DETERMINISTIC: it generates a seeded location-scale
IPD scenario (no NHANES download required), fits IPD-QMA, and prints the
per-quantile DL/REML/PM sensitivity table produced by
``IPDQMA.tau2_sensitivity_profile()``. The DL rows reproduce the engine's
published DL profile exactly, so the table is a diagnostic overlay, not a
re-estimation of the headline result.

Usage
-----
    python benchmark_tau2_sensitivity.py                 # default scenario, pretty table
    python benchmark_tau2_sensitivity.py --k 12 --vr 2.0 # 12 studies, variance ratio 2
    python benchmark_tau2_sensitivity.py --no-hksj       # classical (non-HKSJ) inference
    python benchmark_tau2_sensitivity.py --csv out.csv   # also write the table to CSV
    python benchmark_tau2_sensitivity.py --self-check    # assert DL cell == .fit() profile

The ``--self-check`` mode is a deterministic reproducibility gate: it fails
(exit 1) if the DL sensitivity cell ever diverges from the published DL profile.
"""

from __future__ import annotations

import argparse
import io
import contextlib
import sys

import numpy as np
import pandas as pd

from ipd_qma import IPDQMA, simulate_location_scale


def run_scenario(
    k: int = 15,
    vr: float = 1.5,
    mean_shift: float = 0.0,
    distribution: str = "normal",
    n_boot: int = 300,
    seed: int = 2026,
    use_hksj: bool = True,
    quantiles=None,
    quiet: bool = True,
):
    """Fit IPD-QMA on a seeded location-scale scenario and return the fitted model.

    Deterministic for a fixed ``seed`` (both the data generator and the IPD-QMA
    bootstrap are seeded).
    """
    if quantiles is None:
        quantiles = [0.1, 0.25, 0.5, 0.75, 0.9]

    studies_data, labels, true_params = simulate_location_scale(
        K=k,
        variance_ratio=vr,
        mean_shift=mean_shift,
        distribution=distribution,
        seed=seed,
    )

    qma = IPDQMA(
        quantiles=quantiles,
        n_boot=n_boot,
        conf_level=0.95,
        use_hksj=use_hksj,
        seed=seed,
    )

    if quiet:
        with contextlib.redirect_stdout(io.StringIO()):
            qma.fit(studies_data, labels)
    else:
        qma.fit(studies_data, labels)

    return qma, true_params


def format_table(sens: pd.DataFrame) -> str:
    """Render the sensitivity DataFrame as a fixed-width text table."""
    lines = []
    header = (
        f"{'Quantile':>8} {'Method':>6} {'tau2':>10} {'Effect':>10} "
        f"{'SE':>8} {'CI_Lower':>10} {'CI_Upper':>10} {'I2 (%)':>7}"
    )
    lines.append(header)
    lines.append("-" * len(header))
    for q in sorted(sens['Quantile'].unique()):
        sub = sens[sens['Quantile'] == q]
        for _, r in sub.iterrows():
            lines.append(
                f"{r['Quantile']:>8.2f} {r['Method']:>6} {r['tau2']:>10.5f} "
                f"{r['Effect']:>10.5f} {r['SE']:>8.5f} {r['CI_Lower']:>10.5f} "
                f"{r['CI_Upper']:>10.5f} {r['I2']:>7.1f}"
            )
        lines.append("")
    return "\n".join(lines)


def self_check(qma: IPDQMA, sens: pd.DataFrame) -> bool:
    """Assert the DL sensitivity cell reproduces the published DL profile exactly.

    Returns True on success. This is the reproducibility guarantee: the additive
    sensitivity layer must never perturb the headline DL result.
    """
    profile = qma.results['profile'].reset_index(drop=True)
    dl = sens[sens['Method'] == 'DL'].reset_index(drop=True)
    ok = True
    for i, q in enumerate(qma.quantiles):
        for col in ('Effect', 'tau2', 'CI_Lower', 'CI_Upper'):
            a = float(profile.iloc[i][col])
            b = float(dl.iloc[i][col])
            if abs(a - b) > 1e-10:
                print(
                    f"  MISMATCH at q={q:.2f} {col}: profile={a:.12g} sens_DL={b:.12g}",
                    file=sys.stderr,
                )
                ok = False
    return ok


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Reproducible DL/REML/PM tau^2 sensitivity benchmark for IPD-QMA."
    )
    parser.add_argument("--k", type=int, default=15, help="number of studies (default 15)")
    parser.add_argument("--vr", type=float, default=1.5,
                        help="treatment/control variance ratio (default 1.5)")
    parser.add_argument("--mean-shift", type=float, default=0.0,
                        help="true mean difference (default 0.0, pure scale shift)")
    parser.add_argument("--distribution", default="normal",
                        choices=["normal", "exponential", "lognormal", "t5"],
                        help="outcome distribution (default normal)")
    parser.add_argument("--n-boot", type=int, default=300,
                        help="bootstrap resamples per study (default 300)")
    parser.add_argument("--seed", type=int, default=2026, help="random seed (default 2026)")
    parser.add_argument("--no-hksj", action="store_true",
                        help="use classical (non-HKSJ) inference in the sensitivity pools")
    parser.add_argument("--csv", default=None, help="optional path to write the table as CSV")
    parser.add_argument("--self-check", action="store_true",
                        help="assert DL cell == published DL profile, then exit")
    args = parser.parse_args(argv)

    if args.k < 1:
        parser.error("--k must be >= 1")
    if args.vr <= 0:
        parser.error("--vr must be > 0")
    if args.n_boot < 2:
        parser.error("--n-boot must be >= 2")

    use_hksj = not args.no_hksj

    qma, true_params = run_scenario(
        k=args.k,
        vr=args.vr,
        mean_shift=args.mean_shift,
        distribution=args.distribution,
        n_boot=args.n_boot,
        seed=args.seed,
        use_hksj=use_hksj,
    )
    sens = qma.tau2_sensitivity_profile(use_hksj=use_hksj)

    if args.self_check:
        ok = self_check(qma, sens)
        print("SELF-CHECK PASSED: DL sensitivity cell == published DL profile."
              if ok else "SELF-CHECK FAILED.")
        return 0 if ok else 1

    print("=" * 72)
    print("IPD-QMA tau^2 sensitivity benchmark (DL / REML / Paule-Mandel)")
    print("=" * 72)
    print(f"Scenario: K={args.k}  variance_ratio={args.vr}  mean_shift={args.mean_shift}  "
          f"dist={args.distribution}")
    print(f"Inference: {'HKSJ (t, k-1 df)' if use_hksj else 'classical (normal)'}  "
          f"seed={args.seed}  n_boot={args.n_boot}")
    print("Note: the DL rows below reproduce the engine's published quantile profile;")
    print("      REML/PM are the protocol-mandated sensitivity overlay.\n")
    print(format_table(sens))

    # Compact heterogeneity comparison across estimators (median quantile).
    med = 0.5 if 0.5 in list(qma.quantiles) else qma.quantiles[len(qma.quantiles) // 2]
    sub = sens[sens['Quantile'] == med].set_index('Method')
    print(f"tau^2 at median quantile (q={med:.2f}): "
          f"DL={sub.loc['DL', 'tau2']:.5f}  "
          f"REML={sub.loc['REML', 'tau2']:.5f}  "
          f"PM={sub.loc['PM', 'tau2']:.5f}")

    if args.csv:
        sens.to_csv(args.csv, index=False)
        print(f"\nSensitivity table written to {args.csv}")

    # Always confirm the additive-diagnostic guarantee at the end of a run.
    ok = self_check(qma, sens)
    print("\n[reproducibility] DL cell matches published DL profile: "
          + ("YES" if ok else "NO -- INVESTIGATE"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
