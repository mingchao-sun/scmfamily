# scmfamily 1.0.0 (2026-11-14)

First release. Eleven analyses in the menu group SCM:

- **Synthetic Control** (Abadie, Diamond and Hainmueller, 2010): predictor weights by nested optimisation, in-space placebo test, robustness checks (leave-one-out, treatment timing, pre-period length, in-time placebo). Same code as the standalone `synthcontrol` module.
- **Penalized Synthetic Control**, **Demeaned Synthetic Control**, **Elastic-Net Synthetic Control**, **Bayesian Synthetic Control**, **Interactive Fixed Effects (GSC)**, **Matrix Completion**: permutation (placebo) inference with a pre-fit filter, optional conformal inference where it applies.
- **Synthetic Difference-in-Differences**: port of the R package *synthdid*; SDID, SC and DID estimators, covariates, placebo, bootstrap and jackknife standard errors, several treated units. Reproduces the package exactly.
- **Augmented Synthetic Control**: port of the single-period part of *augsynth*; ridge augmentation with cross-validated penalty, covariates, unit fixed effects, conformal, jackknife+ and jackknife inference. Reproduces the package exactly.
- **Staggered Augmented Synthetic Control**: port of *multisynth*; partial pooling, unit fixed effects, time cohorts, covariates, wild bootstrap and jackknife inference. Reproduces the package exactly.
- **Compare Estimators**: runs up to nine of the estimators on the same panel (outcomes only) and reports ATT, pre- and post-treatment fit, optional uncertainty measures, an overlay of the counterfactual paths, the gaps and the ATT of each method. The numbers equal those of the stand-alone analyses run with their defaults.
- Covariates in Synthetic Control, Penalized, Demeaned, Augmented, Staggered Augmented SC (pre-treatment means) and in Interactive Fixed Effects and SDID (time-varying regressors).
- Every analysis offers *Method notes* and *Data and setup guide* (Output section, last box on each screen).
- Four example datasets: California Proposition 99 (with and without covariates), West Germany reunification, minimum wage and teen employment (staggered adoption).
- Long computations yield to jamovi so that they can be stopped.
- Input checks with specific messages: unbalanced panels, duplicate rows, non-numeric time, unknown excluded donors, first treated period that is not an observed period, too few pre-treatment periods.
- Verified against reference implementations where they exist: Demeaned (*augsynth*), Interactive Fixed Effects (*gsynth*), Matrix Completion (*MCPanel*), Penalized (*pensynth*). Matrix Completion now estimates the fixed effects jointly with the low-rank matrix, as in the paper. Elastic-Net now solves the paper's objective (glmnet's internal scaling had changed the ridge part) and chooses alpha and the penalty by the paper's cross-validation over donors. Penalized SC scales the outcomes by the SD across donors. The Bayesian sampler uses a numerically stable draw and standardised data.
