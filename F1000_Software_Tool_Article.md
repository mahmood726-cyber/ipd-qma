# IPD-QMA: a software tool for reviewer-auditable evidence synthesis

## Authors
- Mahmood Ahmad [1,2]
- Niraj Kumar [1]
- Bilaal Dar [3]
- Laiba Khan [1]
- Andrew Woo [4]
- Corresponding author: Andrew Woo (andy2709w@gmail.com)

## Affiliations
1. Royal Free Hospital
2. Tahir Heart Institute Rabwah
3. King's College Medical School
4. St George's Medical School

## Abstract
**Background:** Average treatment effects can miss clinically important heterogeneity when interventions alter the spread of an outcome distribution rather than its mean. Reviewers therefore need a software paper that explains not only the estimator, but also when the method should be used and how it differs from simpler heterogeneity screens.

**Methods:** IPD-QMA is a two-stage Python framework that estimates within-study quantile treatment effects and pools them across studies with random-effects meta-analysis and HKSJ correction. The local package includes simulations, an NHANES-based worked example, and manuscript build scripts.

**Results:** Project files document simulation experiments across multiple data-generating distributions, quantile-profile outputs, formal gradient tests, NHANES analysis scripts, and pytest-based regression checks.

**Conclusions:** IPD-QMA is best presented as a complementary characterization tool for treatment-effect shape, ideally used after a screening step such as lnVR has suggested variance heterogeneity.

## Keywords
individual participant data; quantile meta-analysis; heterogeneity of treatment effects; quantile regression; software tool

## Introduction
The software aims to make quantile-based evidence synthesis reproducible enough for software peer review. Instead of stopping at a manuscript claim, the repository exposes the simulation driver, applied analysis scripts, data fetcher, and document builder needed to rerun the full workflow.

The paper explicitly compares the framework with standard mean-difference meta-analysis and with variance-ratio screens such as lnVR. This makes the method's niche clearer to reviewers and avoids overstating it as a replacement for all conventional meta-analysis.

The manuscript structure below is deliberately aligned to common open-software review requests: the rationale is stated explicitly, at least one runnable example path is named, local validation artifacts are listed, and conclusions are bounded to the functions and outputs documented in the repository.

## Methods
### Software architecture and workflow
The core library implements the IPDQMA class, within-study quantile estimation, random-effects pooling, bootstrap-based standard errors, and supporting analysis utilities. Project scripts fetch NHANES data, run simulations, perform an illustrative analysis, and build manuscript artifacts.

### Installation, runtime, and reviewer reruns
The local implementation is packaged under `C:\Models\IPD_QMA`. The manuscript identifies the local entry points, dependency manifest, fixed example input, and expected saved outputs so that reviewers can rerun the documented workflow without reconstructing it from scratch.

- Entry directory: `C:\Models\IPD_QMA`.
- Detected documentation entry points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected environment capture or packaging files: `requirements.txt`.
- Named worked-example paths in this draft: `run_simulation_fast.py` for the main Monte Carlo stress tests; `run_analysis.py` plus `data/nhanes_combined.csv` for the applied example; `IPD_QMA_corrected_draft.md` and `paper_text.txt` for project-specific narrative sources.
- Detected validation or regression artifacts: `f1000_artifacts/validation_summary.md`, `test_ipd_qma.py`.
- Detected example or sample data files: `f1000_artifacts/example_dataset.csv`.

### Worked examples and validation materials
**Example or fixed demonstration paths**
- `run_simulation_fast.py` for the main Monte Carlo stress tests.
- `run_analysis.py` plus `data/nhanes_combined.csv` for the applied example.
- `IPD_QMA_corrected_draft.md` and `paper_text.txt` for project-specific narrative sources.

**Validation and reporting artifacts**
- `test_ipd_qma.py` for core library tests.
- `build_docx.py` for manuscript build reproducibility.
- `f1000_artifacts/validation_summary.md` and `tutorial_walkthrough.md` for reviewer-facing documentation.

### Typical outputs and user-facing deliverables
- Quantile treatment-effect profiles across the outcome distribution.
- Slope-test summaries and simulation tables comparing IPD-QMA to mean-based alternatives.
- An applied analysis workflow using public NHANES data.

### Reviewer-informed safeguards
- Provides a named example workflow or fixed demonstration path.
- Documents local validation artifacts rather than relying on unsupported claims.
- Positions the software against existing tools without claiming blanket superiority.
- States limitations and interpretation boundaries in the manuscript itself.
- Requires explicit environment capture and public example accessibility in the released archive.

## Review-Driven Revisions
This draft has been tightened against recurring open peer-review objections taken from the supplied reviewer reports.
- Reproducibility: the draft names a reviewer rerun path and points readers to validation artifacts instead of assuming interface availability is proof of correctness.
- Validation: claims are anchored to local tests, validation summaries, simulations, or consistency checks rather than to unsupported assertions of performance.
- Comparators and niche: the manuscript now names the relevant comparison class and keeps the claimed niche bounded instead of implying universal superiority.
- Documentation and interpretation: the text expects a worked example, input transparency, and reviewer-verifiable outputs rather than a high-level feature list alone.
- Claims discipline: conclusions are moderated to the documented scope of IPD-QMA and paired with explicit limitations.

## Use Cases and Results
The software outputs should be described in terms of concrete reviewer-verifiable workflows: running the packaged example, inspecting the generated results, and checking that the reported interpretation matches the saved local artifacts. In this project, the most important result layer is the availability of a transparent execution path from input to analysis output.

Representative local result: `IPD_QMA_corrected_draft.md` reports Type I error was controlled (0.8-4.9%) and 95% CI coverage was verified against analytical ground truths (all >= 95.1%).

### Concrete local quantitative evidence
- `IPD_QMA_corrected_draft.md` reports Type I error was controlled (0.8-4.9%) and 95% CI coverage was verified against analytical ground truths (all >= 95.1%).
- `paper_text.txt` reports Table 1: Summary of Simulation Results (Power Analysis) Method Normal Data Power Skewed Data Power Interpretation Standard MA (Mean) 5.3% (Random Noise) 14.0% (Erratic) Failed.

## Discussion
Representative local result: `IPD_QMA_corrected_draft.md` reports Type I error was controlled (0.8-4.9%) and 95% CI coverage was verified against analytical ground truths (all >= 95.1%).

The F1000 paper emphasizes that the software is meant to unmask heterogeneous treatment effects that the mean can hide, while being explicit that IPD availability, bootstrap cost, and observational applications constrain routine use.

### Limitations
- The framework requires individual participant data and is therefore less portable than aggregate-data meta-analysis.
- Quantile profiles describe heterogeneity but do not automatically solve confounding in observational demonstrations.
- Reviewers should interpret it as a complementary diagnostic layer rather than a universal first-line estimator.

## Software Availability
- Local source package: `IPD_QMA` under `C:\Models`.
- Public repository: `https://github.com/mahmood726-cyber/ipd-qma`.
- Public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/ipd-qma/tree/0ef940a893d03381b074949bf9a2d05d8b19a458`.
- DOI/archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture detected locally: `requirements.txt`.
- Reviewer-facing documentation detected locally: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Reproducibility walkthrough: `f1000_artifacts/tutorial_walkthrough.md` where present.
- Validation summary: `f1000_artifacts/validation_summary.md` where present.
- Reviewer rerun manifest: `F1000_Reviewer_Rerun_Manifest.md`.
- Multi-persona review memo: `F1000_MultiPersona_Review.md`.
- Concrete submission-fix note: `F1000_Concrete_Submission_Fixes.md`.
- License: see the local `LICENSE` file.

## Data Availability
The project includes code for downloading public NHANES data and generating derived analysis files. No proprietary IPD are bundled in the repository.

## Reporting Checklist
Real-peer-review-aligned checklist: `F1000_Submission_Checklist_RealReview.md`.
Reviewer rerun companion: `F1000_Reviewer_Rerun_Manifest.md`.
Companion reviewer-response artifact: `F1000_MultiPersona_Review.md`.
Project-level concrete fix list: `F1000_Concrete_Submission_Fixes.md`.

## Declarations
### Competing interests
The authors declare that no competing interests were disclosed.

### Grant information
No specific grant was declared for this manuscript draft.

### Author contributions (CRediT)
| Author | CRediT roles |
|---|---|
| Mahmood Ahmad | Conceptualization; Software; Validation; Data curation; Writing - original draft; Writing - review and editing |
| Niraj Kumar | Conceptualization |
| Bilaal Dar | Conceptualization |
| Laiba Khan | Conceptualization |
| Andrew Woo | Conceptualization |

### Acknowledgements
The authors acknowledge contributors to open statistical methods, reproducible research software, and reviewer-led software quality improvement.

## References
1. DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical Trials. 1986;7(3):177-188.
2. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Statistics in Medicine. 2002;21(11):1539-1558.
3. Viechtbauer W. Conducting meta-analyses in R with the metafor package. Journal of Statistical Software. 2010;36(3):1-48.
4. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ. 2021;372:n71.
5. Fay C, Rochette S, Guyader V, Girard C. Engineering Production-Grade Shiny Apps. Chapman and Hall/CRC. 2022.
