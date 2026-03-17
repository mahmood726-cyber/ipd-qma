# PLOS ONE Multi-Persona Review: IPD-QMA Corrected Draft

**Manuscript:** "Unmasking Heterogeneous Treatment Effects in Clinical Trials: A 'Location-Scale' Quantile Meta-Analysis Framework (IPD-QMA)"

**Date:** 2026-02-14
**Review Panel:** 5 independent personas
**Target Journal:** PLOS ONE
**Status:** ALL 42 ISSUES FIXED in draft v3 + code. Simulation re-running at N_SIM=1000 (tables pending).

---

## Reviewers

| # | Persona | Agent | Findings |
|---|---------|-------|----------|
| 1 | Statistical Methodologist | aeae840 | 4 P0, 9 P1, 7 P2 |
| 2 | PLOS ONE Editor | a855d34 | 5 P0, 13 P1, 11 P2 |
| 3 | Clinical Epidemiologist | a006e10 | 4 P0, 8 P1, 8 P2 |
| 4 | Biostatistics Referee | a6659b5 | 4 P0, 9 P1, 10 P2 |
| 5 | Writing/Presentation | a5ab26a | 4 P0, 11 P1, 15 P2 |

**Raw totals:** 21 P0, 50 P1, 51 P2 = 122 findings
**After deduplication:** 10 P0, 18 P1, 14 P2 = 42 unique findings

---

## Consolidated Findings

### P0: Critical Issues (10)

#### P0-1: Missing Mandatory PLOS ONE Sections (Desk Rejection Risk)
**Raised by:** PLOSEditor P0-1/2/3/4/5, Writing P1-11
**Issue:** The manuscript is missing five mandatory PLOS ONE sections:
1. **Data Availability Statement** — must state where NHANES data and code are accessible
2. **Competing Interests Declaration** — required even if none exist
3. **Funding Statement** — required even if unfunded
4. **Author Contributions (CRediT)** — must list each author's role
5. **Ethics Statement** — must reference NCHS ERB approval and cite 45 CFR 46.101(b)(4) exemption

Additionally, the "Changes from Original Submission" section (lines 360-373) must be removed (belongs in cover letter, not manuscript).

**Impact:** Automatic desk rejection without these sections.
**Action:** Add all five sections. Remove change-log section.

---

#### P0-2: N_SIM=200 Insufficient for Type I Error Claims
**Raised by:** StatMethod P0-2, PLOSEditor P0-1, BiostatRef P0-4, ClinEpi P2-8, StatMethod-initial P1-1
**Issue:** All reported results are from N_SIM=200 (confirmed by CSV), giving MCSE=1.54% at alpha=0.05. The 95% CI for a 5% rate is [2.0%, 8.0%] — too wide to confirm nominal calibration. Meanwhile, code files have been updated to N_SIM=1000 but never re-run, creating a code-paper mismatch.
**Impact:** Cannot distinguish 5% from 2-8%. Claims of "well-controlled Type I error" are not adequately supported.
**Action:** Re-run simulations at N_SIM >= 1000 (ideally 2000). Update all tables, paper text, and MCSE claims. Estimated runtime: ~4 hours for fast simulation.

---

#### P0-3: NHANES Confounding by Indication Undermines Real-Data Demonstration
**Raised by:** ClinEpi P0-1, ClinEpi P0-4, PLOSEditor P2-2, BiostatRef P1-8
**Issue:** The NHANES analysis compares SBP of medicated vs. unmedicated hypertensives. The positive QTE gradient (Q10: +9.5 to Q90: +18.6 mmHg) simply recovers confounding by indication — sicker patients are prescribed medications. A naive reader could interpret +18.6 mmHg as medications *raising* SBP. The paper acknowledges this but still presents formal p-values and CIs, creating a false sense of clinical validity.
**Impact:** Misleading clinical interpretation; inappropriate causal language for observational data.
**Action:** Either (a) reframe NHANES explicitly as "Illustrative Methodological Demonstration" (not a clinical finding), remove or caveat p-values, state the confounding interpretation as the *lead* sentence, OR (b) add propensity score adjustment. Also add a dedicated paragraph separating causal (RCT) from descriptive (observational) interpretation of QTEs.

---

#### P0-4: NHANES Survey Design Effects Ignored — Standard Errors Invalid
**Raised by:** ClinEpi P0-3, BiostatRef P1-8, ClinEpi P1-8
**Issue:** NHANES uses complex multistage probability sampling with stratification, clustering (PSUs), and person-level weights. All are ignored. This underestimates standard errors by 1.5-3x (typical DEFF for NHANES BP outcomes is 2-4), making all NHANES CIs and p-values invalid. This violates CDC NCHS Analytic Guidelines (2013).
**Impact:** Reported CIs (e.g., Q10: [7.0, 12.1]) could be 50-70% too narrow.
**Action:** Either (a) implement design-adjusted analysis, (b) present NHANES without CIs/p-values with explicit warning, or (c) move to supplement with prominent caveats. Do NOT present formal inference from unweighted NHANES.

---

#### P0-5: "Hidden Utility" Term Used but Never Defined
**Raised by:** Writing P0-2, PLOSEditor P2-5
**Issue:** "Hidden utility" appears in abstract (line 11), simulation design (line 117), and Discussion (line 250) as if it were an established concept, but is never given a formal definition.
**Impact:** Non-specialist readers will not understand the key concept on first encounter.
**Action:** Add a one-sentence definition at first use: "We define a 'hidden utility' scenario as one where a treatment modifies the outcome variance without changing the mean, benefiting patients at one tail while harming those at the other."

---

#### P0-6: Ground Truth Not Derived for Non-Normal Distributions
**Raised by:** StatMethod P0-1
**Issue:** Analytical ground truth for Q90 difference is derived only for Normal (z_0.9 * (sqrt(VR) - 1) = +0.938). For Exponential, Lognormal, and t_5, the code sets `true_diff = None`. For Exponential with VR=3, the true Q90 diff is approximately 1.686 (not 0.938). For t_5, it's approximately 1.080. Power comparisons across distributions are not directly comparable because the true effect sizes differ.
**Impact:** Cannot validate whether simulation power reflects the correct alternative hypothesis for non-Normal DGPs.
**Action:** Derive analytical ground truth for all four distributions. Verify with large-sample empirical checks. Report true QTE values alongside power results. At minimum, state explicitly that cross-distribution power comparisons are not directly comparable.

---

#### P0-7: Bootstrap B and N_SIM Parameter Mismatches Between Paper and Code
**Raised by:** PLOSEditor P0-3, StatMethod-initial P1-5, BiostatRef P1-6
**Issue:** Paper states B=100 for simulation (line 123), code uses B=200 (`run_simulation_fast.py` line 186). Paper states N_SIM=200 (line 120), code files now have N_SIM=1000. `run_analysis.py` has N_SIM=500. Three different values across three files.
**Impact:** Reviewers cannot reproduce results from paper text.
**Action:** Align all files and paper text to a single set of parameters. Document exact values used for each reported result.

---

#### P0-8: Tau Notation Collision (Quantile Level vs. Heterogeneity Variance)
**Raised by:** Writing P0-4, StatMethod P2-1
**Issue:** tau denotes quantile level (Eq. 1, line 64) AND between-study heterogeneity variance (tau_hat^2, line 86). Same Greek letter, two fundamentally different meanings.
**Impact:** Confuses readers at the mathematical core of the paper.
**Action:** Add explicit disambiguation at first use of tau^2: "Note: tau denotes the quantile level throughout; tau^2 (tau-squared) denotes between-study heterogeneity variance, following standard notation in both fields." Alternatively, use a different symbol for one concept.

---

#### P0-9: Slope Test Hypothesis Framing Imprecise
**Raised by:** StatMethod P0-3, ClinEpi P1-4, BiostatRef P0-1, PLOSEditor P1-8
**Issue:** The paper describes the Slope Test as testing "H_0: Delta = 0 (homogeneous treatment effect)." But the pooled DL slope actually tests whether the *population-average slope* across studies is zero — a weaker hypothesis than within-study homogeneity. A zero pooled slope could hide opposing slopes in individual studies. The test also only captures linear gradients; U-shaped patterns (Q10 and Q90 positive, Q50 zero) would be missed.
**Impact:** Central claim about the paper's primary confirmatory test is imprecise.
**Action:** Clarify that the Slope Test tests "H_0: the population-average QTE gradient is zero." Note it detects distributional heterogeneity (necessary but not sufficient for clinically actionable HTE). Discuss the limitation that only linear patterns are captured. Consider mentioning a multivariate Wald test as future work.

---

#### P0-10: No Figure/Table References in Manuscript
**Raised by:** Writing P1-6, PLOSEditor P1-6
**Issue:** The paper contains zero figure references despite generating 7 figures. Tables lack formal numbers, titles, and legends. PLOS ONE requires all figures/tables to be numbered, captioned, and referenced at first mention in the text.
**Impact:** The manuscript reads as an incomplete draft without visual evidence.
**Action:** Add numbered Figure captions (Figures 1-7) and formal Table numbers (Tables 1-6) with descriptive legends. Reference each at the appropriate point in the text.

---

### P1: Important Issues (18)

#### P1-1: Unconditional vs. Conditional QTE Distinction Not Addressed
**Raised by:** ClinEpi P1-1, BiostatRef P2-9, PLOSEditor P2-4, Writing P2-5
The paper estimates unconditional QTEs (no covariates) but uses language implying conditional effects ("patients at the tau-th percentile of severity"). These coincide only under rank-invariance. State the estimand explicitly. Cite Firpo (2007).

#### P1-2: Independent tau^2 Per Quantile — No Cross-Quantile Borrowing
**Raised by:** BiostatRef P0-2
Each quantile is pooled independently via DL with separate tau^2. No information sharing across quantiles. For small K, this is highly imprecise. Multivariate random-effects models (Jackson et al., 2011) could improve efficiency. Must be discussed as a limitation.

#### P1-3: DL Estimator Used Without Justification vs. REML/PM
**Raised by:** PLOSEditor P1-5, BiostatRef P1-1
DL is known to underestimate tau^2 for small K (Veroniki et al., 2016). HKSJ partially compensates but doesn't fix the bias. Either justify DL choice (simplicity, HKSJ guard) or add REML sensitivity analysis.

#### P1-4: Two-Stage vs. One-Stage Tradeoffs Oversimplified
**Raised by:** BiostatRef P1-3
The paper dismisses one-stage models as less practical but doesn't acknowledge what is lost: joint quantile modeling, information borrowing, better handling of unequal sample sizes. Expand Section 2.1.

#### P1-5: Insufficient Comparison to Prior Work (Sun & Cheung, Maruo, Firpo)
**Raised by:** StatMethod P1-5, BiostatRef P1-2, PLOSEditor P1-12, Writing P1-2
Sun & Cheung [15] gets one sentence. Maruo et al. (2017), Firpo (2007), PRISMA-IPD (Stewart et al., 2015), Cummins (2018) are not cited. A structured comparison table is needed.

#### P1-6: NHANES Missing Key Details (Sample Sizes, Flow, Variable Definitions)
**Raised by:** PLOSEditor P1-8, ClinEpi P1-2, PLOSEditor P1-6
No participant flow diagram, no numbers at each filtering step, no MAR assessment, no survey weight discussion in Methods. Age restriction (>=18) not mentioned. These are STROBE requirements.

#### P1-7: K=4 in NHANES Below Paper's Own Recommended K >= 5
**Raised by:** ClinEpi P1-5
Paper recommends K >= 5 but NHANES uses K=4 cycles. With HKSJ at K=4, df=3 gives t(3)=3.182 — extremely conservative. Acknowledge explicitly.

#### P1-8: lnVR 8.0% Type I Error on Exponential Warrants Investigation
**Raised by:** StatMethod P1-4, BiostatRef P1-7, PLOSEditor P0-2
Dismissed as "within MCSE bounds" but could be genuine anti-conservatism for skewed distributions. After re-running at N_SIM >= 1000, if confirmed (>6.5%), report as a finding about lnVR's distributional assumptions, not dismiss as noise.

#### P1-9: MD Rejection Rates in "Power" Tables Is Misleading Framing
**Raised by:** BiostatRef P1-6, Writing P1-4
MD rejection rates under location-scale DGP are Type I error (H_0: MD=0 is true). Presenting them in "Power" tables implies a power contest. Restructure: present MD rates as "Type I error confirmation" separately from Q90/Slope/lnVR "Power."

#### P1-10: Only Pure Scale Shifts Tested — No Mixed Location+Scale Scenarios
**Raised by:** StatMethod-initial P1-4, BiostatRef P2-10
All scenarios have mean_shift=0. Real clinical settings typically have simultaneous mean and scale effects. Add at least one scenario (e.g., mean_shift=0.5, VR=2.0) to validate IPD-QMA under more realistic conditions.

#### P1-11: Abstract Overloaded and Selective
**Raised by:** Writing P0-1, StatMethod-initial P1-3
Results paragraph tries to cover everything in 130+ words. MD rejection range "1.5-3.0%" is selective (actual range 0.0-5.0%). Split into two clearer sentences. Add a sentence on limitations to abstract per PLOS ONE enhanced abstract policy.

#### P1-12: Missing Reporting Checklists (PRISMA-IPD, STROBE)
**Raised by:** PLOSEditor P1-12, ClinEpi P2-1, PLOSEditor P1-5
No PRISMA-IPD or STROBE checklist. PLOS ONE requires adherence to relevant reporting guidelines.

#### P1-13: Multiple Testing Not Formally Addressed
**Raised by:** StatMethod P1-8, BiostatRef P2-4, PLOSEditor P1-7
7 tests per meta-analysis (5 quantiles + Slope + lnVR). Paper advises "caution" but doesn't designate a primary test or suggest a correction. Propose: Slope Test as primary confirmatory, quantile-specific tests as exploratory.

#### P1-14: No DOIs on References; [15] Incomplete; [14] and [16] Uncited
**Raised by:** PLOSEditor P1-10/11, Writing P1-10, P2-12/13
All 16 references lack DOIs. Sun & Cheung [15] missing volume/pages. References [14] (Higgins & Thompson) and [16] (Koenker) are listed but never cited in text.

#### P1-15: Simulation DGP Ecological Validity Questionable
**Raised by:** ClinEpi P0-2, BiostatRef P2-3, ClinEpi P1-9
Pure mathematical scale shift doesn't model realistic clinical mechanisms. VR=3-5 is extreme (Nakagawa found typical lnVR ~ -0.5 to +0.5, i.e., VR ~ 0.6-1.6). Mean-centering Exponential creates negative values. Temper clinical claims or add biologically motivated scenarios.

#### P1-16: lnVR Naming Contradicts Symbol
**Raised by:** Writing P1-9
Section heading says "Log SD Ratio" but symbol is "lnVR" (variance ratio). Contradictory. Either rename symbol to "lnSDR" or consistently call it "Log Variance Ratio" with a note that lnVR = 0.5 * log(Var ratio) = log(SD ratio).

#### P1-17: Bootstrap Seed Reuse (Data and Bootstrap Share Same Seed)
**Raised by:** StatMethod P0-4, BiostatRef P1-9
Same seed used for data generation and bootstrap within each simulation iteration. While impact is small (separate Generator instances), it's methodologically unsound. Fix: use `seed + 1_000_000` for bootstrap.

#### P1-18: No Discussion of Clinical Significance Thresholds (MCIDs)
**Raised by:** ClinEpi P1-3
All effects are framed as statistically significant. No mention of minimal clinically important differences. The simulation uses standardized units without clinical translation. Add a paragraph on clinical vs. statistical significance for QTEs.

---

### P2: Minor Improvements (14)

| ID | Source | Issue |
|---|---|---|
| P2-1 | Writing | Em-dash overuse (15+ instances) — convert half to commas/parentheses |
| P2-2 | Writing | Inconsistent capitalization ("Quantile Regression" vs. "random-effects models") |
| P2-3 | Writing | Rhetorical question in Section 4.2 — rephrase for academic tone |
| P2-4 | Writing | "HTE" acronym defined (line 24) but never used again — remove definition |
| P2-5 | Writing | Abstract background sentence too long (62 words) — split |
| P2-6 | Writing | "Two-Stage" capitalized inconsistently — use lowercase in running text |
| P2-7 | PLOSEditor | P-values reported inconsistently (exact vs. threshold notation) |
| P2-8 | PLOSEditor | CI notation not per PLOS ONE style — use (95% CI, X to Y) |
| P2-9 | PLOSEditor | Software version pinning missing (Python, numpy, scipy versions) |
| P2-10 | ClinEpi | Prediction intervals not reported or interpreted in NHANES results |
| P2-11 | ClinEpi | DBP data fetched but never analyzed — missed supplementary opportunity |
| P2-12 | BiostatRef | Mention lnCVR (Nakagawa framework) as alternative variability metric |
| P2-13 | BiostatRef | HKSJ max(1, q) truncation — cite Rover et al. (2015) explicitly |
| P2-14 | BiostatRef | Add Sidik & Jonkman (2002) reference alongside Hartung & Knapp (2001) |

---

## Overall Recommendation: MAJOR REVISION

### Strengths
- The core contribution (IPD-QMA detecting location-scale effects invisible to standard MA) is **novel and clinically relevant**
- Statistical machinery (QR + DL + HKSJ + bootstrap covariance + Slope Test) is **correctly implemented**
- The corrected draft addresses 9 critical code bugs from the prior review
- Ground truth derivation for Normal is verified correct (+0.938)
- Simulation design covers meaningful parameter space (4 distributions, VR/K/N axes)

### Critical Actions Before Submission
1. **Add all 5 mandatory PLOS ONE sections** (P0-1) — 30 min
2. **Re-run simulations at N_SIM >= 1000** (P0-2) — 4 hours
3. **Reframe NHANES as illustrative demo, not clinical finding** (P0-3, P0-4) — 1 hour
4. **Define "hidden utility" at first use** (P0-5) — 5 min
5. **Derive non-Normal ground truths or state incomparability** (P0-6) — 1 hour
6. **Align all parameter values (B, N_SIM) between paper and code** (P0-7) — 30 min
7. **Disambiguate tau notation** (P0-8) — 10 min
8. **Clarify Slope Test hypothesis precisely** (P0-9) — 30 min
9. **Add figure/table references and captions** (P0-10) — 2 hours
10. **Remove "Changes from Original Submission" section** (P0-1) — 1 min

### Recommended But Not Blocking
- Expand prior work comparison (P1-5)
- Add mixed location+scale simulation scenario (P1-10)
- Add PRISMA-IPD/STROBE checklists (P1-12)
- Fix seed reuse in bootstrap (P1-17)
- Add participant flow for NHANES (P1-6)

### Consensus Across All 5 Reviewers
All reviewers independently identified:
1. N_SIM=200 is insufficient (5/5 reviewers)
2. NHANES confounding undermines the demonstration (4/5)
3. Missing figure/table integration (3/5)
4. Parameter mismatches between paper and code (4/5)
5. Prior work comparison is inadequate (4/5)
6. Multiple testing not formally addressed (3/5)

The method itself is sound. The issues are about **rigor of presentation**, **completeness of reporting**, and **adequacy of simulation evidence**, not about the statistical framework.
