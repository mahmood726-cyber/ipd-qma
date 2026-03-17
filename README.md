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
| `run_simulation_fast.py` | Monte Carlo simulation (Type I error, power, coverage) |
| `run_analysis.py` | NHANES applied analysis and figure generation |
| `fetch_nhanes.py` | Downloads and preprocesses NHANES data |
| `build_docx.py` | Generates manuscript .docx from analysis outputs |
| `test_ipd_qma.py` | Pytest test suite |

## Running Tests

```bash
python -m pytest test_ipd_qma.py -v
```

## Citation

If you use this software, please cite:

> [Author Names]. IPD-QMA: Individual Participant Data Quantile Meta-Analysis. [Journal]. [Year]. [DOI].

## License

MIT License. See [LICENSE](LICENSE) for details.
