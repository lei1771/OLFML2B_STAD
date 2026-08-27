# Reproducibility and statistical contract

## Inference hierarchy

The analysis attempts to preserve the patient as the principal inferential unit. Cells and spatial spots are nested measurements and are not counted as independent patient replicates.

## Bulk and clinical analysis

- OLFML2B is analysed continuously within cohort after standardisation for survival models.
- Median/high-low splits are used for visualisation only and are not promoted as clinical thresholds.
- Common-covariate Cox models use verified age, sex and overall stage where structurally available.
- Structurally absent covariates are not forced into imputation for displayed or pooled estimates.
- Proportional hazards are checked using Schoenfeld-based diagnostics; time-varying sensitivity analyses are retained where relevant.
- Random-effects synthesis uses cohort-level estimates and reports uncertainty/heterogeneity rather than only pooled P values.

## Ecological programmes

OLFML2B is removed from any programme containing it before programme scoring/association analysis, preventing mechanical self-correlation. Programme associations are interpreted as shared RNA/ecological state, not causal mediation or validated absolute cell fractions.

## Single-cell analysis

Datasets are not jointly re-clustered solely to create apparent replication. Cross-dataset conclusions rely on harmonised compartment concepts, dataset/patient summaries, paired pseudobulk where formal mappings exist, and leave-one-dataset-out checks.

## Spatial analysis

All sections may be displayed, but repeated sections are collapsed within patient for formal cross-patient summaries. Spatial results support local contextualisation/association, not cell identity, cell-cell interaction or causality.

## Proteomics

The PDC analysis preserves same-patient, same-plex pairing and audits plex sensitivity and leave-one-plex-out behaviour.

## Immunotherapy

PRJEB25780 is used to define the boundary of a treatment-prediction claim. The workflow does not convert a small single cohort into an unsupported multivariable clinical predictor.

## Multiplicity and negative results

Analysis-family FDR control, effect sizes, confidence intervals, prediction intervals, exact/small-sample sensitivity tests and negative results are retained. The project deliberately avoids post hoc rescue by alternative cut-points when the locked analysis is inconclusive.
