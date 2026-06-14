"""
harness.py -- Truth-recovery yardstick for ipd-qma.

IPD-QMA estimates quantile treatment effects (QTEs) per study (quantile regression
+ bootstrap SE) and pools them across studies (DL + HKSJ). The honest test: inject
a KNOWN QTE truth and measure whether the pooled per-quantile CI covers it -- with
special attention to the TAIL quantiles (q=0.1, 0.9), where quantile estimation is
hardest, vs the median.

Two truths:
  - location shift: treatment = control + delta  ->  QTE(q) = delta for ALL q.
  - scale change:   treatment ~ N(0, s^2)         ->  QTE(q) = (s-1)*Phi^{-1}(q),
    i.e. the QTE varies by quantile (zero at the median, large in the tails).

Uses the app's own IPDQMA, run unchanged (n_boot reduced for simulation speed;
coverage is insensitive to n_boot at this resolution). Truth-first: every number
is produced from seeded simulation here.
Run:  python truth-recovery/harness.py --reps 60
"""
import sys, os, argparse, io, contextlib
import numpy as np
from scipy.stats import norm
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ipd_qma import IPDQMA

BASE_SEED = 20260613
QUANTILES = [0.1, 0.5, 0.9]
N_BOOT = 40


def true_qte(q, scenario, delta, s):
    if scenario == "location":
        return delta
    return (s - 1.0) * norm.ppf(q)   # scale change


def gen(rng, k, scenario, delta, s, n_arm=90, tau=0.0):
    studies = []
    for _ in range(k):
        shift = rng.normal(0, tau)              # between-study heterogeneity in the effect
        control = rng.standard_normal(n_arm)
        if scenario == "location":
            treatment = rng.standard_normal(n_arm) + delta + shift
        else:
            treatment = rng.standard_normal(n_arm) * s + shift
        studies.append((control, treatment))
    return studies


def run_cell(k, scenario, delta, s, reps, seed0, tau=0.0):
    cov = {q: 0 for q in QUANTILES}
    n = 0
    for r in range(reps):
        rng = np.random.default_rng(seed0 + r)
        studies = gen(rng, k, scenario, delta, s, tau=tau)
        eng = IPDQMA(quantiles=QUANTILES, n_boot=N_BOOT, use_hksj=True, seed=seed0 + r)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                res = eng.fit(studies)
        except Exception:
            continue
        prof = {float(row["Quantile"]): row for _, row in res["profile"].iterrows()}
        n += 1
        for q in QUANTILES:
            row = prof.get(q)
            if row is None:
                continue
            t = true_qte(q, scenario, delta, s)
            if row["CI_Lower"] <= t <= row["CI_Upper"]:
                cov[q] += 1
    return {q: round(cov[q] / max(1, n), 3) for q in QUANTILES}, n


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--reps", type=int, default=60)
    reps = ap.parse_args().reps
    import time; t0 = time.time()
    print(f"\n# Truth-recovery yardstick -- ipd-qma (quantile treatment effects)")
    print(f"reps={reps}/cell  k=5  n_boot={N_BOOT}  quantiles={QUANTILES}  seed={BASE_SEED}\n")
    print("coverage of the TRUE QTE at each quantile (should be ~0.95)\n")
    cells = [
        ("location shift d=0.5, homog", dict(k=5, scenario="location", delta=0.5, s=1.0, tau=0.0)),
        ("location shift d=0.5, heterog", dict(k=5, scenario="location", delta=0.5, s=1.0, tau=0.15)),
        ("scale change s=1.5 (tail QTEs)", dict(k=5, scenario="scale", delta=0.0, s=1.5, tau=0.0)),
    ]
    print(f"{'scenario':34s} | " + " ".join(f"q={q}" for q in QUANTILES))
    for label, cfg in cells:
        cov, n = run_cell(reps=reps, seed0=BASE_SEED, **cfg)
        print(f"{label:34s} | " + "  ".join(f"{cov[q]:.3f}" for q in QUANTILES) + f"   (n={n})")
    print(f"\n(tail quantiles q=0.1/0.9 are the hard cases. {time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
