# E156 Protocol — `ipd-qma`

This repository is the source code and dashboard backing an E156 micro-paper on the [E156 Student Board](https://mahmood726-cyber.github.io/e156/students.html).

---

## `[73]` IPD-QMA: Detecting Heterogeneous Treatment Effects via Quantile Meta-Analysis of Individual Participant Data

**Type:** methods  |  ESTIMAND: Quantile treatment effect  
**Data:** Monte Carlo simulation (4 distributions, ADEMP) + NHANES application

### 156-word body

Can quantile meta-analysis of individual participant data detect treatment effect heterogeneity that standard mean-difference pooling misses? IPD-QMA combines within-study quantile regression with DerSimonian-Laird random-effects pooling and Hartung-Knapp-Sidik-Jonkman correction across four distributional scenarios following the ADEMP simulation reporting framework. The two-stage framework estimates unconditional quantile treatment effects at five grid points via random-effects models, with a novel Slope Test formally contrasting effects at the 90th versus 10th percentile. Under pure scale shift with variance ratio at least 2.0, IPD-QMA detected the mean difference with 100 percent power and valid 95% CI coverage at K equals 10 while standard meta-analysis correctly retained the null. Type I error remained controlled between 0.8 and 4.9 percent, and coverage was verified against analytical ground truths at 95.1 percent or above. IPD-QMA complements existing tools by characterising treatment effect shape across the full outcome distribution. However, a limitation is that rank invariance required for individual-level causal interpretation remains untestable in observational settings.

### Submission metadata

```
Corresponding author: Mahmood Ahmad <mahmood.ahmad2@nhs.net>
ORCID: 0000-0001-9107-3704
Affiliation: Tahir Heart Institute, Rabwah, Pakistan

Links:
  Code:      https://github.com/mahmood726-cyber/ipd-qma
  Protocol:  https://github.com/mahmood726-cyber/ipd-qma/blob/main/E156-PROTOCOL.md
  Dashboard: https://mahmood726-cyber.github.io/ipd-qma/

References (topic pack: individual participant data (IPD) meta-analysis):
  1. Riley RD, Lambert PC, Abo-Zaid G. 2010. Meta-analysis of individual participant data: rationale, conduct, and reporting. BMJ. 340:c221. doi:10.1136/bmj.c221
  2. Burke DL, Ensor J, Riley RD. 2017. Meta-analysis using individual participant data: one-stage and two-stage approaches, and why they may differ. Stat Med. 36(5):855-875. doi:10.1002/sim.7141

Data availability: No patient-level data used. Analysis derived exclusively
  from publicly available aggregate records. All source identifiers are in
  the protocol document linked above.

Ethics: Not required. Study uses only publicly available aggregate data; no
  human participants; no patient-identifiable information; no individual-
  participant data. No institutional review board approval sought or required
  under standard research-ethics guidelines for secondary methodological
  research on published literature.

Funding: None.

Competing interests: MA serves on the editorial board of Synthēsis (the
  target journal); MA had no role in editorial decisions on this
  manuscript, which was handled by an independent editor of the journal.

Author contributions (CRediT):
  [STUDENT REWRITER, first author] — Writing – original draft, Writing –
    review & editing, Validation.
  [SUPERVISING FACULTY, last/senior author] — Supervision, Validation,
    Writing – review & editing.
  Mahmood Ahmad (middle author, NOT first or last) — Conceptualization,
    Methodology, Software, Data curation, Formal analysis, Resources.

AI disclosure: Computational tooling (including AI-assisted coding via
  Claude Code [Anthropic]) was used to develop analysis scripts and assist
  with data extraction. The final manuscript was human-written, reviewed,
  and approved by the author; the submitted text is not AI-generated. All
  quantitative claims were verified against source data; cross-validation
  was performed where applicable. The author retains full responsibility for
  the final content.

Preprint: Not preprinted.

Reporting checklist: PRISMA 2020 (methods-paper variant — reports on review corpus).

Target journal: ◆ Synthēsis (https://www.synthesis-medicine.org/index.php/journal)
  Section: Methods Note — submit the 156-word E156 body verbatim as the main text.
  The journal caps main text at ≤400 words; E156's 156-word, 7-sentence
  contract sits well inside that ceiling. Do NOT pad to 400 — the
  micro-paper length is the point of the format.

Manuscript license: CC-BY-4.0.
Code license: MIT.

SUBMITTED: [ ]
```


---

_Auto-generated from the workbook by `C:/E156/scripts/create_missing_protocols.py`. If something is wrong, edit `rewrite-workbook.txt` and re-run the script — it will overwrite this file via the GitHub API._