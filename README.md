# scmfamily — synthetic control in jamovi

A jamovi module with ten synthetic-control estimators and a comparison analysis for comparative case studies and panel data. Chinese summary: [README.zh.md](README.zh.md). Japanese summary: [README.ja.md](README.ja.md).

## Analyses (jamovi menu: SCM)

All analyses share the same inputs: a long-format panel (unit,
time, outcome), a treated unit and its first treated period.

| Analysis | Method | Reference | Inference |
|---|---|---|---|
| Synthetic Control | Classic SCM, nested-optimised predictor weights | Abadie, Diamond and Hainmueller (2010) | In-space placebos; leave-one-out, timing, pre-period and in-time robustness |
| Penalized Synthetic Control | Simplex weights with a distance penalty | Abadie and L'Hour (2021) | Placebos; conformal |
| Demeaned Synthetic Control | Weights fitted to demeaned outcomes | Ferman and Pinto (2021) | Placebos; conformal |
| Augmented Synthetic Control | Simplex SC with ridge correction (port of *augsynth*) | Ben-Michael, Feller and Rothstein (2021) | Conformal, time jackknife+, jackknife |
| Elastic-Net Synthetic Control | Unrestricted weights, elastic net | Doudchenko and Imbens (2016) | Placebos; conformal |
| Bayesian Synthetic Control | Horseshoe-shrinkage regression | Kim, Lee and Gupta (2020) | Credible interval; placebos |
| Interactive Fixed Effects (GSC) | Factor model with optional covariates | Xu (2017) | Placebos; conformal |
| Matrix Completion | Nuclear-norm regularised imputation | Athey et al. (2021) | Placebos |
| Synthetic Difference-in-Differences | Unit and time weights (port of *synthdid*) | Arkhangelsky et al. (2021) | Placebo, bootstrap, jackknife |
| Staggered Augmented Synthetic Control | Partial pooling over treated units (port of *multisynth*) | Ben-Michael, Feller and Rothstein (2022) | Wild bootstrap, jackknife |
| Compare Estimators | Runs up to nine of the estimators above on one panel and shows them side by side | as above | Each method's own measure, optional |

## Data

One row per unit and period, a **balanced** panel (every unit in every period, no duplicate
unit-period rows, no missing outcomes). Unit: Nominal or ID; time: numeric (year or period index).
The first treated period must be an observed period, with at least two pre-treatment periods and two
untreated donors. Every analysis has a checkbox *Data and setup guide* (Output section) that shows the
exact settings for the sample data, and a checkbox *Method notes* that explains the method and every
option on the screen.

Example data (`data/`, sources in `data/README.md`):

| Dataset | Use it for |
|---|---|
| California Proposition 99 (with and without covariates) | all single-treated-unit analyses |
| West Germany reunification | covariates in Penalized, Demeaned, Augmented SC |
| Minimum wage and teen employment (staggered) | Staggered Augmented SC (use time cohorts, switch off unit fixed effects) |

## Covariates

* Synthetic Control, Penalized, Demeaned, Augmented and Staggered Augmented SC: pre-treatment means,
  balanced together with the outcomes. A balance table is reported.
* Interactive Fixed Effects and SDID: time-varying regressors with common coefficients. They must be
  observed in every period for every unit and must not be affected by the treatment (otherwise they
  absorb part of the effect).
* Elastic-Net, Bayesian and Matrix Completion use outcomes only.

## Validation

### Agreement with the original packages

SDID, Augmented SC and Staggered Augmented SC are ports of *synthdid*, *augsynth* and *multisynth*
(`R/sdid_engine.R`, `R/augsynth_engine.R`, `R/multisynth_engine.R`). Estimates, weights, penalty choices
and standard errors were compared with the official packages on the California data (SDID, Augmented SC)
and on simulated staggered panels (Staggered Augmented SC); in every setting checked, including bootstrap
and jackknife results with the same random seed, the differences are 0.
`augsynth` options that are not ported: `progfunc = "gsyn"`, permutation inference, multiple outcomes,
matching on covariates (exact / knn) and the time-regression outcome model of `multisynth`.

The other analyses were compared with reference implementations on the California data
(`tests/test_methods.R` keeps the reference values):

| Analysis | Reference | Result |
|---|---|---|
| Demeaned SC | *augsynth* (unit fixed effects, no outcome model) | weights and ATT equal (differences below 1e-8) |
| Interactive Fixed Effects | *gsynth* (two-way fixed effects, fixed r, no covariates) | ATT and counterfactual path equal (below 1e-12) |
| Matrix Completion | *MCPanel* `mcnnm_fit` (same penalty) | counterfactual equal (below 1e-4); the earlier two-step version did not match and was replaced |
| Penalized SC | *pensynth* (its lambda is half of ours) | weights agree to solver tolerance (about 1e-5) |
| Elastic-Net | exact solution of the paper's objective (independent solver) and the paper's cross-validation | weights equal (below 1e-4) |
| Bayesian SC | an independent implementation of the same Gibbs sampler | posterior means agree within Monte Carlo error; *bayesreg*'s horseshoe regression gave different posterior means when donors outnumber periods, a difference not explained |

Not compared with a reference package: the covariate versions of Penalized, Demeaned and IFE (checked on
simulated data only), and the conformal inference, which follows Chernozhukov, Wuthrich and Zhu (2021) but
has no reference run.

## Install
- **jamovi library**: in jamovi open the modules button (top right), *jamovi library*, find *Synthetic Control Family* and install (once the module is listed).
- **Sideload**: download the `.jmo` file and use *Install from file* in the same menu. A `.jmo` works only on the operating system, processor type and jamovi version it was built for.
- **From source** (needs R and the jamovi compiler `jmvtools`; adjust the path to your jamovi folder):
```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")
```

## Use from R (Rj editor, scripts)
Every analysis is also an R function, for example
```r
scmfamily::scm(data = d, outcome = "cigsale", unit = "state", time = "year", treatedUnit = "California", treatmentTime = 1989)
```

## References
Papers and R packages behind each analysis appear in the References section of the jamovi output (defined in `jamovi/00refs.yaml`).

## Known limitations

* Several analyses can take a minute or more: SDID standard errors (especially with covariates),
  conformal inference, the jackknife in Staggered Augmented SC. They can be stopped from jamovi.
* The Bayesian prior follows Kim et al. (2020) in spirit (horseshoe on unconstrained donor coefficients); check it
  against the paper before relying on it.
* Matrix Completion uses a penalty set by hand; the paper chooses it by cross-validation.
* The Interactive Fixed Effects analysis is the control-only version; it does not choose the number of factors.

## Notes for developers
- Engine tests (run from the module folder): `Rscript tests/smoke.R`, `tests/test_methods.R`, `tests/test_enet_bscm.R`, `tests/test_reference.R` (reference packages in the local R); required R packages: jmvcore, R6, ggplot2, quadprog, osqp, Matrix, glmnet.
- Layout: the ports of *synthdid*, *augsynth* and *multisynth* are in `R/sdid_engine.R`, `R/augsynth_engine.R`, `R/multisynth_engine.R`; `R/*.b.R` and `jamovi/*.yaml` are the jamovi front end.
- Ribbon icon: `icon/`; jmvtools 28.17 cannot attach an icon to a module, so `dev/icon_preview.ps1` (if present) adds one CSS rule for the menu group to jamovi's client stylesheet (administrator PowerShell, restart jamovi afterwards).

## How to cite
Sun, M. *scmfamily: Synthetic Control Family for jamovi* (jamovi module). https://github.com/mingchao-sun/scmfamily. Please also cite the papers behind the estimator you use (listed in the output).

## License
GPL (>= 2): version 2 of the GNU General Public License, or any later version (`LICENSE` holds version 2). This is compatible with *glmnet* (GPL-2 only), which the Elastic-Net analysis calls.
Third-party code and data are listed in `NOTICE.md`. Bug reports and suggestions: https://github.com/mingchao-sun/scmfamily/issues
