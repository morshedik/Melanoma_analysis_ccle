# Project Strengthening Plan

## Strong core idea

Connecting baseline molecular state to drug response across melanoma cell lines can
generate hypotheses about response biomarkers and resistant states.

## Why the present analysis is not yet sufficient

1. A correlation across a short list of drug-target pairs is not a pharmacogenomic
   response model.
2. CCLE and GDSC identifiers require an explicit, audited cell-line mapping.
3. A single A375 profile cannot support across-cell-line association statistics.
4. Pathway enrichment of the most highly expressed genes lacks a matched background
   and can reflect generic housekeeping programs.
5. Package installation, data download, analysis, and reporting are currently mixed.
6. Several scripts depend on objects created interactively by earlier scripts.

## Stronger analysis

1. Freeze the CCLE and GDSC releases and checksum the raw files.
2. Build an explicit mapping table using stable model identifiers; review ambiguous
   mappings manually.
3. For each drug, model response across melanoma cell lines using expression features,
   lineage covariates, and nested cross-validation.
4. Keep all preprocessing and feature selection inside training folds.
5. Compare elastic net against simple baselines.
6. Report effect sizes, confidence intervals, sample counts, missingness, and
   out-of-fold predictions.
7. Validate shortlisted biomarkers in a separate dataset or perturbation experiment.
8. Replace interactive object sharing with one parameterized entry script.

## Biological upgrade

Focus the project on one mechanistic question rather than all drugs and pathways. A
strong example structure is:

- registered pathway or resistance state;
- a small, biologically justified drug class;
- cross-cell-line prediction;
- held-out validation; and
- a proposed experiment that can falsify the mechanism.

No clinical or therapeutic claim should be made from the current exploratory scripts.
