# IPD-QMA: concrete submission fixes

This file converts the multi-persona review into repository-side actions that should be checked before external submission of the F1000 software paper for `IPD_QMA`.

## Detectable Local State
- Documentation files detected: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Environment lock or container files detected: `requirements.txt`.
- Package manifests detected: none detected.
- Example data files detected: `f1000_artifacts/example_dataset.csv`.
- Validation artifacts detected: `f1000_artifacts/validation_summary.md`, `test_ipd_qma.py`.
- Detected public repository root: `https://github.com/mahmood726-cyber/ipd-qma`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/ipd-qma/tree/0ef940a893d03381b074949bf9a2d05d8b19a458`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.

## High-Priority Fixes
- Check that the manuscript's named example paths exist in the public archive and can be run without repository archaeology.
- Confirm that the cited repository root (`https://github.com/mahmood726-cyber/ipd-qma`) resolves to the same fixed public source snapshot used for submission.
- Archive the tagged release and insert the Zenodo DOI or record URL once it has been minted; no project-specific archive DOI was detected locally.
- Reconfirm the quoted benchmark or validation sentence after the final rerun so the narrative text matches the shipped artifacts.

## Numeric Evidence Available To Quote
- `IPD_QMA_corrected_draft.md` reports Type I error was controlled (0.8-4.9%) and 95% CI coverage was verified against analytical ground truths (all >= 95.1%).
- `paper_text.txt` reports Table 1: Summary of Simulation Results (Power Analysis) Method Normal Data Power Skewed Data Power Interpretation Standard MA (Mean) 5.3% (Random Noise) 14.0% (Erratic) Failed.

## Manuscript Files To Keep In Sync
- `F1000_Software_Tool_Article.md`
- `F1000_Reviewer_Rerun_Manifest.md`
- `F1000_MultiPersona_Review.md`
- `F1000_Submission_Checklist_RealReview.md` where present
- README/tutorial files and the public repository release metadata
