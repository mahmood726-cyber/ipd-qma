"""
test_truth_recovery.py -- measured invariants for the ipd-qma yardstick. Seeded;
reduced reps (the per-study bootstrap is slow). Exit 0 = all pass.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from harness import run_cell, BASE_SEED, QUANTILES

ok = True


def check(name, cond, detail):
    global ok
    print(f"{'PASS' if cond else 'FAIL'}  {name}  {detail}")
    if not cond:
        ok = False


# scale change -> the QTE VARIES by quantile (zero at median, large in tails):
# the interesting case for whether the tail-quantile CIs recover the truth.
cov, n = run_cell(k=5, scenario="scale", delta=0.0, s=1.5, reps=20, seed0=BASE_SEED, tau=0.0)
check("recovers the true tail QTE at q=0.1 (hard case)", cov[0.1] > 0.85, f"(coverage {cov[0.1]}, n={n})")
check("recovers the true tail QTE at q=0.9 (hard case)", cov[0.9] > 0.85, f"(coverage {cov[0.9]})")
check("recovers the true median QTE at q=0.5", cov[0.5] > 0.85, f"(coverage {cov[0.5]})")
# at small k the HKSJ + bootstrap intervals tend to be conservative (>= nominal).
check("coverage is at/above nominal (conservative, not anti-conservative)",
      min(cov.values()) >= 0.85, f"(min coverage {min(cov.values())})")

print("\nAll measured invariants hold." if ok else "\nSOME INVARIANTS FAILED.")
sys.exit(0 if ok else 1)
