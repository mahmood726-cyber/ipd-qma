Mahmood Ahmad
Tahir Heart Institute
author@example.com

IPD-QMA: Detecting Heterogeneous Treatment Effects via Quantile Meta-Analysis of Individual Participant Data

Can quantile meta-analysis of individual participant data detect treatment effect heterogeneity that standard mean-difference pooling misses? IPD-QMA combines within-study quantile regression with DerSimonian-Laird random-effects pooling and Hartung-Knapp-Sidik-Jonkman correction across four distributional scenarios following the ADEMP simulation reporting framework. The two-stage framework estimates unconditional quantile treatment effects at five grid points via random-effects models, with a novel Slope Test formally contrasting effects at the 90th versus 10th percentile. Under pure scale shift with variance ratio at least 2.0, IPD-QMA detected the mean difference with 100 percent power and valid 95% CI coverage at K equals 10 while standard meta-analysis correctly retained the null. Type I error remained controlled between 0.8 and 4.9 percent, and coverage was verified against analytical ground truths at 95.1 percent or above. IPD-QMA complements existing tools by characterising treatment effect shape across the full outcome distribution. However, a limitation is that rank invariance required for individual-level causal interpretation remains untestable in observational settings.

Outside Notes

Type: methods
Primary estimand: Quantile treatment effect
App: IPD-QMA Python library
Data: Monte Carlo simulation (4 distributions, ADEMP) + NHANES application
Code: https://github.com/mahmood726-cyber/ipd-qma
Version: 1.0
Validation: DRAFT

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI (Claude, Anthropic) was used as a constrained synthesis engine operating on structured inputs and predefined rules for infrastructure generation, not as an autonomous author. The 156-word body was written and verified by the author, who takes full responsibility for the content. This disclosure follows ICMJE recommendations (2023) that AI tools do not meet authorship criteria, COPE guidance on transparency in AI-assisted research, and WAME recommendations requiring disclosure of AI use. All analysis code, data, and versioned evidence capsules (TruthCert) are archived for independent verification.
