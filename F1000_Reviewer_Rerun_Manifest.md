# IPD-QMA: reviewer rerun manifest

This manifest is the shortest reviewer-facing rerun path for the local software package. It lists the files that should be sufficient to recreate one worked example, inspect saved outputs, and verify that the manuscript claims remain bounded to what the repository actually demonstrates.

## Reviewer Entry Points
- Project directory: `C:\Models\IPD_QMA`.
- Preferred documentation start points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected public repository root: `https://github.com/mahmood726-cyber/ipd-qma`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/ipd-qma/tree/0ef940a893d03381b074949bf9a2d05d8b19a458`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture files: `requirements.txt`.
- Validation/test artifacts: `f1000_artifacts/validation_summary.md`, `test_ipd_qma.py`.

## Worked Example Inputs
- Manuscript-named example paths: `run_simulation_fast.py` for the main Monte Carlo stress tests; `run_analysis.py` plus `data/nhanes_combined.csv` for the applied example; `IPD_QMA_corrected_draft.md` and `paper_text.txt` for project-specific narrative sources; f1000_artifacts/example_dataset.csv.
- Auto-detected sample/example files: `f1000_artifacts/example_dataset.csv`.

## Expected Outputs To Inspect
- Quantile treatment-effect profiles across the outcome distribution.
- Slope-test summaries and simulation tables comparing IPD-QMA to mean-based alternatives.
- An applied analysis workflow using public NHANES data.

## Minimal Reviewer Rerun Sequence
- Start with the README/tutorial files listed below and keep the manuscript paths synchronized with the public archive.
- Create the local runtime from the detected environment capture files if available: `requirements.txt`.
- Run at least one named example path from the manuscript and confirm that the generated outputs match the saved validation materials.
- Quote one concrete numeric result from the local validation snippets below when preparing the final software paper.

## Local Numeric Evidence Available
- `IPD_QMA_corrected_draft.md` reports Type I error was controlled (0.8-4.9%) and 95% CI coverage was verified against analytical ground truths (all >= 95.1%).
- `paper_text.txt` reports Table 1: Summary of Simulation Results (Power Analysis) Method Normal Data Power Skewed Data Power Interpretation Standard MA (Mean) 5.3% (Random Noise) 14.0% (Erratic) Failed.
