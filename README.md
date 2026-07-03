# IPD-QMA: Individual Participant Data Quantile Meta-Analysis

A two-stage framework combining within-study quantile regression with DerSimonian-Laird random-effects pooling and Hartung-Knapp-Sidik-Jonkman (HKSJ) correction. IPD-QMA detects treatment effect heterogeneity across the outcome distribution that standard mean-difference meta-analysis misses.

## Installation

Requires Python 3.11+.

```bash
pip install -r requirements.txt
```

## Reproduction Steps

1. **Download NHANES data** (~5 min):
   ```bash
   python fetch_nhanes.py
   ```

2. **Run Monte Carlo simulation** (~30 min, outputs to `output/`):
   ```bash
   python run_simulation_fast.py
   ```

3. **Run NHANES analysis and generate figures**:
   ```bash
   python run_analysis.py
   ```

4. **Generate manuscript .docx**:
   ```bash
   python build_docx.py
   ```

Or use the Makefile to run all steps:

```bash
make all
```

## Key Files

| File | Description |
|------|-------------|
| `ipd_qma.py` | Core library: IPDQMA class, DL pooling, HKSJ correction, quantile regression |
| `tau2_estimators.py` | Between-study variance estimators (DL / REML / Paule-Mandel), R-validated |
| `benchmark_tau2_sensitivity.py` | Worked example + reproducible DL/REML/PM sensitivity benchmark |
| `run_simulation_fast.py` | Monte Carlo simulation (Type I error, power, coverage) |
| `run_analysis.py` | NHANES applied analysis and figure generation |
| `fetch_nhanes.py` | Downloads and preprocesses NHANES data |
| `build_docx.py` | Generates manuscript .docx from analysis outputs |
| `test_ipd_qma.py` / `test_tau2_estimators.py` | Pytest test suite |

## tau^2 Estimator Sensitivity (DL / REML / Paule-Mandel)

The engine pools each quantile's study-level effects with the DerSimonian-Laird
(DL) moment estimator. DL is known to underestimate between-study variance for
small K; the pre-registered protocol names REML as the intended primary
estimator and calls for a DL/REML/PM sensitivity analysis. This is available as
an **additive diagnostic** — it does not change the published DL-based estimates.

After fitting, obtain the per-quantile sensitivity table directly from the engine:

```python
from ipd_qma import IPDQMA, simulate_location_scale

studies, labels, _ = simulate_location_scale(K=15, variance_ratio=1.5, seed=2026)
qma = IPDQMA(seed=2026)
qma.fit(studies, labels)

sens = qma.tau2_sensitivity_profile()   # DataFrame: one row per (quantile, DL/REML/PM)
print(sens)
```

The `DL` rows reproduce the `.fit()` quantile profile to machine precision, so
the table is a comparison overlay rather than a re-estimation of the headline
result.

### Reproducible benchmark

`benchmark_tau2_sensitivity.py` is a self-contained, deterministic worked example
(no NHANES download needed). It generates a seeded location-scale scenario, fits
IPD-QMA, and prints the DL/REML/PM sensitivity table:

```bash
python benchmark_tau2_sensitivity.py                  # default scenario, pretty table
python benchmark_tau2_sensitivity.py --k 12 --vr 2.0  # 12 studies, variance ratio 2
python benchmark_tau2_sensitivity.py --no-hksj        # classical (non-HKSJ) inference
python benchmark_tau2_sensitivity.py --csv table.csv  # also write the table to CSV
python benchmark_tau2_sensitivity.py --self-check     # reproducibility gate (exit 1 on mismatch)
```

`--self-check` asserts the DL sensitivity cell still equals the published DL
profile and exits non-zero if it ever diverges — a machine-checkable guarantee
that the sensitivity layer stays additive.

## Running Tests

Run the full suite (core engine + tau^2 estimators):

```bash
python -m pytest -v
```

## Citation

If you use this software, please cite:

> [Author Names]. IPD-QMA: Individual Participant Data Quantile Meta-Analysis. [Journal]. [Year]. [DOI].

## License

MIT License. See [LICENSE](LICENSE) for details.
