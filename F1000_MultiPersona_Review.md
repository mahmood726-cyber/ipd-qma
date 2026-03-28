# IPD-QMA: multi-persona peer review

This memo applies the recurring concerns in the supplied peer-review document to the current F1000 draft for this project (`IPD_QMA`). It distinguishes changes already made in the draft from repository-side items that still need to hold in the released repository and manuscript bundle.

## Detected Local Evidence
- Detected documentation files: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected environment capture or packaging files: `requirements.txt`.
- Detected validation/test artifacts: `f1000_artifacts/validation_summary.md`, `test_ipd_qma.py`.
- Detected browser deliverables: no HTML file detected.
- Detected public repository root: `https://github.com/mahmood726-cyber/ipd-qma`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/ipd-qma/tree/0ef940a893d03381b074949bf9a2d05d8b19a458`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.

## Reviewer Rerun Companion
- `F1000_Reviewer_Rerun_Manifest.md` consolidates the shortest reviewer-facing rerun path, named example files, environment capture, and validation checkpoints.

## Detected Quantitative Evidence
- `IPD_QMA_corrected_draft.md` reports Type I error was controlled (0.8-4.9%) and 95% CI coverage was verified against analytical ground truths (all >= 95.1%).
- `paper_text.txt` reports Table 1: Summary of Simulation Results (Power Analysis) Method Normal Data Power Skewed Data Power Interpretation Standard MA (Mean) 5.3% (Random Noise) 14.0% (Erratic) Failed.

## Current Draft Strengths
- States the project rationale and niche explicitly: Average treatment effects can miss clinically important heterogeneity when interventions alter the spread of an outcome distribution rather than its mean. Reviewers therefore need a software paper that explains not only the estimator, but also when the method should be used and how it differs from simpler heterogeneity screens.
- Names concrete worked-example paths: `run_simulation_fast.py` for the main Monte Carlo stress tests; `run_analysis.py` plus `data/nhanes_combined.csv` for the applied example; `IPD_QMA_corrected_draft.md` and `paper_text.txt` for project-specific narrative sources.
- Points reviewers to local validation materials: `test_ipd_qma.py` for core library tests; `build_docx.py` for manuscript build reproducibility; `f1000_artifacts/validation_summary.md` and `tutorial_walkthrough.md` for reviewer-facing documentation.
- Moderates conclusions and lists explicit limitations for IPD-QMA.

## Remaining High-Priority Fixes
- Keep one minimal worked example public and ensure the manuscript paths match the released files.
- Ensure README/tutorial text, software availability metadata, and public runtime instructions stay synchronized with the manuscript.
- Confirm that the cited repository root resolves to the same fixed public source snapshot used for the submission package.
- Mint and cite a Zenodo DOI or record URL for the tagged release; none was detected locally.
- Reconfirm the quoted benchmark or validation sentence after the final rerun so the narrative text stays synchronized with the shipped artifacts.

## Persona Reviews

### Reproducibility Auditor
- Review question: Looks for a frozen computational environment, a fixed example input, and an end-to-end rerun path with saved outputs.
- What the revised draft now provides: The revised draft names concrete rerun assets such as `run_simulation_fast.py` for the main Monte Carlo stress tests; `run_analysis.py` plus `data/nhanes_combined.csv` for the applied example and ties them to validation files such as `test_ipd_qma.py` for core library tests; `build_docx.py` for manuscript build reproducibility.
- What still needs confirmation before submission: Before submission, freeze the public runtime with `requirements.txt` and keep at least one minimal example input accessible in the external archive.

### Validation and Benchmarking Statistician
- Review question: Checks whether the paper shows evidence that outputs are accurate, reproducible, and compared against known references or stress tests.
- What the revised draft now provides: The manuscript now cites concrete validation evidence including `test_ipd_qma.py` for core library tests; `build_docx.py` for manuscript build reproducibility; `f1000_artifacts/validation_summary.md` and `tutorial_walkthrough.md` for reviewer-facing documentation and frames conclusions as being supported by those materials rather than by interface availability alone.
- What still needs confirmation before submission: Concrete numeric evidence detected locally is now available for quotation: `IPD_QMA_corrected_draft.md` reports Type I error was controlled (0.8-4.9%) and 95% CI coverage was verified against analytical ground truths (all >= 95.1%); `paper_text.txt` reports Table 1: Summary of Simulation Results (Power Analysis) Method Normal Data Power Skewed Data Power Interpretation Standard MA (Mean) 5.3% (Random Noise) 14.0% (Erratic) Failed.

### Methods-Rigor Reviewer
- Review question: Examines modeling assumptions, scope conditions, and whether method-specific caveats are stated instead of implied.
- What the revised draft now provides: The architecture and discussion sections now state the method scope explicitly and keep caveats visible through limitations such as The framework requires individual participant data and is therefore less portable than aggregate-data meta-analysis; Quantile profiles describe heterogeneity but do not automatically solve confounding in observational demonstrations.
- What still needs confirmation before submission: Retain method-specific caveats in the final Results and Discussion and avoid collapsing exploratory thresholds or heuristics into universal recommendations.

### Comparator and Positioning Reviewer
- Review question: Asks what gap the tool fills relative to existing software and whether the manuscript avoids unsupported superiority claims.
- What the revised draft now provides: The introduction now positions the software against an explicit comparator class: The paper explicitly compares the framework with standard mean-difference meta-analysis and with variance-ratio screens such as lnVR. This makes the method's niche clearer to reviewers and avoids overstating it as a replacement for all conventional meta-analysis.
- What still needs confirmation before submission: Keep the comparator discussion citation-backed in the final submission and avoid phrasing that implies blanket superiority over better-established tools.

### Documentation and Usability Reviewer
- Review question: Looks for a README, tutorial, worked example, input-schema clarity, and short interpretation guidance for outputs.
- What the revised draft now provides: The revised draft points readers to concrete walkthrough materials such as `run_simulation_fast.py` for the main Monte Carlo stress tests; `run_analysis.py` plus `data/nhanes_combined.csv` for the applied example; `IPD_QMA_corrected_draft.md` and `paper_text.txt` for project-specific narrative sources and spells out expected outputs in the Methods section.
- What still needs confirmation before submission: Make sure the public archive exposes a readable README/tutorial bundle: currently detected files include `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.

### Software Engineering Hygiene Reviewer
- Review question: Checks for evidence of testing, deployment hygiene, browser/runtime verification, secret handling, and removal of obvious development leftovers.
- What the revised draft now provides: The draft now foregrounds regression and validation evidence via `f1000_artifacts/validation_summary.md`, `test_ipd_qma.py`, and browser-facing projects are described as self-validating where applicable.
- What still needs confirmation before submission: Before submission, remove any dead links, exposed secrets, or development-stage text from the public repo and ensure the runtime path described in the manuscript matches the shipped code.

### Claims-and-Limitations Editor
- Review question: Verifies that conclusions are bounded to what the repository actually demonstrates and that limitations are explicit.
- What the revised draft now provides: The abstract and discussion now moderate claims and pair them with explicit limitations, including The framework requires individual participant data and is therefore less portable than aggregate-data meta-analysis; Quantile profiles describe heterogeneity but do not automatically solve confounding in observational demonstrations; Reviewers should interpret it as a complementary diagnostic layer rather than a universal first-line estimator.
- What still needs confirmation before submission: Keep the conclusion tied to documented functions and artifacts only; avoid adding impact claims that are not directly backed by validation, benchmarking, or user-study evidence.

### F1000 and Editorial Compliance Reviewer
- Review question: Checks for manuscript completeness, software/data availability clarity, references, and reviewer-facing support files.
- What the revised draft now provides: The revised draft is more complete structurally and now points reviewers to software availability, data availability, and reviewer-facing support files.
- What still needs confirmation before submission: Confirm repository/archive metadata, figure/export requirements, and supporting-file synchronization before release.
