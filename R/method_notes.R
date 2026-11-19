# Method notes shown in the results when 'Method notes' is ticked (all analyses except the
# classic Synthetic Control, whose notes are in scm_help_text.R). Each note explains the
# method, every option on the analysis screen, and how to read the results.

fam_notes_vars <- function(covariates = NULL, treated_extra = FALSE) {
    paste0(
        "<h2>Data</h2>",
        "<p>The data must be a <strong>balanced panel in long format</strong>: one row per unit and period, ",
        "with an outcome value for every unit in every period.</p><ul>",
        "<li><strong>Outcome</strong>: numeric variable to analyse.</li>",
        "<li><strong>Unit</strong>: identifier of the cross-sectional unit (for example, state name).</li>",
        "<li><strong>Time</strong>: period variable, numeric years or consecutive period indices.</li>",
        "<li><strong>Treated unit</strong>: the unit that receives the treatment. All other units form the donor pool.</li>",
        if (treated_extra) "<li><strong>Additional treated units</strong>: further units treated in the same first period, comma-separated. Their outcomes are averaged with the main treated unit.</li>",
        "<li><strong>First treated period</strong>: the first period in which the treatment is in effect. Earlier periods are pre-treatment; this period and later ones are post-treatment.</li>",
        "<li><strong>Exclude donors</strong>: comma-separated units removed from the donor pool in all analyses (for example units with confounding policies). The treated unit cannot be excluded.</li>",
        covariates, "</ul>")
}

fam_notes_placebo <- function() paste0(
    "<h2>Permutation (placebo) inference</h2>",
    "<p>Each donor is treated as if it were the treated unit and the same estimator is re-run on the remaining donors. ",
    "The treated unit's result is compared with this placebo distribution.</p><ul>",
    "<li><strong>RMSPE-ratio p-value (one-sided)</strong>: share of units, counting the treated unit, whose post/pre-treatment RMSPE ratio is at least as large as the treated unit's.</li>",
    "<li><strong>|ATT| p-value (two-sided)</strong>: the same with the absolute average effect.</li>",
    "<li><strong>Pre-fit filter (max RMSPE ratio)</strong>: placebo units whose pre-treatment RMSPE exceeds this multiple of the treated unit's are left out of the filtered p-value, because a poor pre-treatment fit makes the ratio uninformative. The default is 2.</li></ul>",
    "<p>The smallest attainable p-value is 1 divided by the number of units, so power is low with few donors. The placebo table lists every unit's ratio, pre-treatment RMSPE and ATT.</p>")

fam_notes_conformal <- function(avail = "") paste0(
    "<h2>Conformal inference</h2>",
    "<p>Chernozhukov, Wuthrich and Zhu (2021). The effect is fixed at a null value, the counterfactual is re-estimated on all periods, ",
    "and the size of the post-treatment residuals is compared with cyclic rotations of the residual series. The reported p-value tests an average effect of zero, ",
    "and the 95% interval is obtained by inverting the test over a grid of effect values. ", avail, "</p>")

fam_notes_output <- function() paste0(
    "<h2>Results and output options</h2><ul>",
    "<li><strong>Estimate and fit</strong>: the average post-treatment effect (ATT), pre- and post-treatment RMSPE and their ratio.</li>",
    "<li><strong>Donor weights</strong>: weight on each donor, largest first.</li>",
    "<li><strong>Treated vs counterfactual</strong> and <strong>Gap</strong> plots: the outcome path against the estimated counterfactual, and their difference over time; placebo gaps appear in grey when permutation inference is on.</li>",
    "<li><strong>Period-by-period table</strong>: treated outcome, counterfactual and gap for every period.</li>",
    "<li><strong>Method notes</strong> (this text) and <strong>Data and setup guide</strong>: shown at the top of the results when ticked.</li></ul>")

fam_notes_html <- function(key) {
    H <- function(title, intro) paste0("<h1>", title, "</h1><p>", intro, "</p>")
    cov_pre <- paste0("<li><strong>Covariates</strong> (optional): each covariate is averaged over the pre-treatment periods for every unit (missing years are skipped), ",
                      "standardised across units and balanced together with the pre-treatment outcomes. <strong>Covariate weight</strong> sets the weight of the whole covariate block relative to the outcome block: ",
                      "0 ignores the covariates, 1 gives both blocks equal total weight.</li>")
    switch(key,

    penalized = paste0(
        H("Penalized Synthetic Control", "Abadie and L'Hour (2021). Donor weights are non-negative and sum to one, and they minimise the squared pre-treatment gap to the treated unit plus a penalty on the distance between the treated unit and each donor with positive weight. The penalty discourages weights on donors that resemble the treated unit only in combination, and makes the solution unique."),
        fam_notes_vars(cov_pre),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Penalty lambda (standardised units)</strong>: outcomes are standardised by their SD across the donors, and lambda is in those units. It multiplies the sum, over donors, of weight times squared distance to the treated unit, in the objective of the paper (squared gap plus lambda times that sum). The <em>pensynth</em> package halves the squared gap, so its lambda is half of this one. Lambda = 0 gives the ordinary synthetic control; larger values move towards nearest-neighbour matching. Report results for several values.</li>",
        "<li><strong>Covariate weight (relative to outcomes)</strong>: see Covariates above.</li></ul>",
        fam_notes_placebo(), fam_notes_conformal("Available for this estimator."), fam_notes_output(),
        "<h2>Notes and limitations</h2><ul><li>With many pre-treatment periods, outcome lags usually carry most of the information, and covariates change the weights little.</li><li>Permutation inference reruns the estimator for every donor and takes longer with many donors.</li></ul>"),

    demeaned = paste0(
        H("Demeaned Synthetic Control", "Ferman and Pinto (2021). Each unit's pre-treatment mean is removed before the simplex weights are fitted, so donors can match the treated unit's trend even when their levels differ. The counterfactual adds back the treated unit's pre-treatment mean."),
        fam_notes_vars(cov_pre),
        "<h2>Estimation</h2><ul>",
        "<li>There is no tuning parameter for the weights; they solve a constrained least-squares problem on the demeaned pre-treatment outcomes. Without covariates the result equals <em>augsynth</em> with unit fixed effects and no outcome model.</li>",
        "<li><strong>Covariate weight (relative to outcomes)</strong>: see Covariates above.</li></ul>",
        fam_notes_placebo(), fam_notes_conformal("Available for this estimator."), fam_notes_output(),
        "<h2>Notes and limitations</h2><ul><li>Use this version when the treated unit's level lies outside the range of the donors but its trend is similar.</li><li>The demeaned estimator is not invariant to which pre-treatment periods are used, so check that the result is stable when the window changes.</li></ul>"),

    ife = paste0(
        H("Interactive Fixed Effects (Generalized Synthetic Control)", "Xu (2017). Two-way fixed effects and r latent factors are estimated from the donor pool only. The treated unit's factor loadings come from a regression of its pre-treatment outcomes on the factors, and the counterfactual is built from them. This is the control-only version; it does not cross-validate the number of factors."),
        fam_notes_vars("<li><strong>Covariates</strong> (optional): time-varying regressors with one coefficient common to all units, estimated jointly with the factors by alternating least squares and SVD. They must be observed in <em>every period for every unit</em> and must not be affected by the treatment, otherwise they absorb part of the effect. Their estimates appear in the covariate coefficient table.</li>"),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Number of factors (r)</strong>: at least 1 and at most the number of pre-treatment periods minus 2. Compare several values; a value that is too small leaves common shocks unexplained, one that is too large over-fits the pre-treatment period. For a fixed r and no covariates the estimates equal those of the <em>gsynth</em> package (two-way fixed effects, no cross-validation).</li></ul>",
        fam_notes_placebo(), fam_notes_conformal("Available for this estimator."), fam_notes_output(),
        "<h2>Notes and limitations</h2><ul><li>No donor weights are reported because the estimator is a factor model, not a weighted average of donors.</li><li>Factors are only identified up to rotation, so compare the fit and the effect, not the factors themselves.</li></ul>"),

    enet = paste0(
        H("Elastic-Net Synthetic Control", "Doudchenko and Imbens (2016). Weights are unrestricted: they may be negative, need not sum to one, and an intercept is included. The weights come from an elastic-net regression of the treated unit's pre-treatment outcomes on the donors' outcomes."),
        fam_notes_vars(NULL),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Choose alpha by cross-validation (0.1 to 0.9)</strong>: on (default), alpha and the penalty lambda are chosen jointly as in the paper: each donor in turn is treated as the pseudo-treated unit, its pre-treatment outcomes are fitted on the other donors, and the squared prediction error in the post-treatment periods is averaged over donors; the pair with the smallest error is used. Off, alpha is fixed and only lambda is chosen this way.</li>",
"<li><strong>Mixing alpha (0 ridge, 1 lasso)</strong>: used when cross-validation of alpha is off. 0 is pure ridge (many small weights), 1 is pure lasso (few non-zero weights).</li></ul>",
"<p>The objective is the paper's: the sum of squared pre-treatment gaps plus lambda times (half of (1 minus alpha) times the squared L2 norm of the weights plus alpha times the L1 norm), with an unpenalised intercept and no scaling of the donors. The summary reports the selected alpha, the selected lambda on the paper's scale, and the intercept.</p>",
        fam_notes_placebo(), fam_notes_conformal("Available for this estimator; the penalty and alpha selected for the main fit are kept fixed inside the test."), fam_notes_output(),
        "<h2>Notes and limitations</h2><ul><li>This analysis uses the outcome only; it has no covariate option.</li><li>Negative weights extrapolate beyond the donors. With few pre-treatment periods the estimate can be unstable; compare with the simplex-weight methods.</li></ul>"),

    bscm = paste0(
        H("Bayesian Synthetic Control", "Kim, Lee and Gupta (2020). The treated unit's pre-treatment outcomes are regressed on the donors with horseshoe shrinkage priors on the coefficients, no simplex constraint and a free intercept. The posterior is sampled by Gibbs sampling."),
        fam_notes_vars(NULL),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Posterior draws (after burn-in)</strong>: number of retained draws; the same number is discarded as burn-in. More draws give smoother intervals and take longer.</li></ul>",
        "<p>The horseshoe prior (half-Cauchy global and local scales, Jeffreys prior on the error variance) is placed on the coefficients of the donors and the outcome standardised to unit SD, and the coefficients are transformed back to the donors' own scale. The ATT is the posterior mean of the average post-treatment gap, reported with a 95% credible interval; the interval reflects both the uncertainty of the coefficients and the sampling error of the post-treatment outcomes (error variance divided by the number of post-treatment periods). The weights table shows posterior-mean coefficients, which need not be non-negative or sum to one. The sampler was checked against an independent implementation of the same Gibbs scheme.</p>",
        fam_notes_placebo(), fam_notes_conformal("The option appears on every screen but is computed only for the penalized, demeaned, elastic-net and interactive fixed-effects estimators; for this analysis, ticking it only adds a note."), fam_notes_output(),
        "<h2>Notes and limitations</h2><ul><li>This analysis uses the outcome only; it has no covariate option.</li><li>The prior specification follows the paper in spirit; check it against the paper before publication.</li></ul>"),

    mc = paste0(
        H("Matrix Completion", "Athey, Bayati, Doudchenko, Imbens and Khosravi (2021). The treated unit's post-treatment outcomes are treated as missing and imputed from a low-rank matrix plus unit and time fixed effects, estimated by nuclear-norm regularisation (soft-impute)."),
        fam_notes_vars(NULL),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Nuclear-norm penalty (relative)</strong>: the penalty as a fraction of the smallest penalty that sets the low-rank part to zero (1 = no low-rank part, small values = close to a perfect fit of the observed cells). Larger values give a lower-rank fit. The paper chooses the penalty by cross-validation; here it is set by hand, so try several values. For a given value the estimate equals that of the authors' <em>MCPanel</em> package at convergence.</li></ul>",
        fam_notes_placebo(), fam_notes_conformal("The option appears on every screen but is computed only for the penalized, demeaned, elastic-net and interactive fixed-effects estimators; for this analysis, ticking it only adds a note."), fam_notes_output(),
        "<h2>Notes and limitations</h2><ul><li>No donor weights are reported because the estimator is not a weighted average of donors.</li><li>This analysis uses the outcome only; it has no covariate option.</li></ul>"),

    sdid = paste0(
        H("Synthetic Difference-in-Differences", "Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021). This analysis is a port of the R package synthdid and reproduces its estimates and standard errors. Unit weights (omega) make the controls resemble the treated units before treatment, time weights (lambda) make the pre-treatment periods resemble the post-treatment ones, and the effect is the weighted difference-in-differences."),
        fam_notes_vars("<li><strong>Covariates</strong> (optional): time-varying regressors with common coefficients, estimated jointly with the weights; the effect is computed on the outcome net of their linear effect. They must be observed in every period for every unit and must not respond to the treatment.</li>", treated_extra = TRUE),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Estimator</strong>: synthetic difference-in-differences; synthetic control (unit weights only, with a small ridge penalty and no intercept); or difference-in-differences (uniform weights). All three run in the same framework so that they can be compared.</li></ul>",
        "<h2>Inference</h2><ul>",
        "<li><strong>Standard error</strong>: <em>Placebo</em> (Algorithm 4) reassigns treatment to random controls and is the only choice with a single treated unit; <em>Bootstrap</em> (Algorithm 2) and <em>Jackknife</em> (Algorithm 3) need several treated units; <em>None</em> skips it.</li>",
        "<li><strong>Replications</strong>: number of placebo or bootstrap replications.</li>",
        "<li><strong>Random seed</strong>: makes the replications reproducible.</li></ul>",
        "<p>The confidence interval uses a normal approximation, estimate plus or minus 1.96 standard errors.</p>",
        "<h2>Results and output options</h2><ul>",
        "<li><strong>Unit weights (omega)</strong> and <strong>time weights (lambda)</strong>: only non-zero weights are listed.</li>",
        "<li><strong>Effect curve</strong>: the effect for each post-treatment period; the ATT is its average.</li>",
        "<li><strong>Covariate coefficients</strong>: shown when covariates are selected.</li>",
        "<li><strong>Period-by-period table</strong>, <strong>Method notes</strong> and <strong>Data and setup guide</strong>: optional output.</li></ul>",
        "<h2>Notes and limitations</h2><ul><li>Treated units must adopt in the same period. For staggered adoption use the Staggered Augmented SC analysis.</li><li>Placebo standard errors need more control units than treated units.</li><li>Standard errors re-estimate the weights in every replication and can be slow, especially with covariates (minutes for 200 replications). Try 20 replications first; the run can be stopped from jamovi.</li></ul>"),

    ascm = paste0(
        H("Augmented Synthetic Control", "Ben-Michael, Feller and Rothstein (2021). This analysis is a port of the single-period part of the R package augsynth and reproduces its results. Simplex synthetic control weights are corrected for imperfect pre-treatment fit by a ridge regression of outcomes on pre-treatment outcomes."),
        fam_notes_vars("<li><strong>Covariates</strong> (optional): auxiliary covariates averaged over the pre-treatment periods (missing values omitted). By default they are standardised to the scale of the outcomes and balanced together with them; with <em>Residualize covariates</em> they are instead regressed out of the pre-treatment outcomes.</li>", treated_extra = TRUE),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Outcome model</strong>: <em>Ridge augmentation</em> adds the ridge correction; <em>None</em> gives the plain synthetic control.</li>",
        "<li><strong>Use synthetic control weights</strong>: if unticked, uniform weights are used and only the ridge model adjusts the estimate.</li>",
        "<li><strong>Unit fixed effects</strong>: demeans each unit by its pre-treatment mean.</li>",
        "<li><strong>Residualize covariates</strong>: see Covariates above.</li>",
        "<li><strong>Ridge penalty</strong>: blank chooses it by leave-out cross-validation over the pre-treatment periods.</li>",
        "<li><strong>Largest lambda within one SE of the CV minimum</strong>: the more conservative choice; untick to use the minimiser.</li></ul>",
        "<h2>Inference</h2><ul>",
        "<li><strong>Inference</strong>: <em>Conformal</em> (default) gives period-by-period intervals and p-values; <em>Time jackknife+</em> gives prediction intervals by leaving out pre-treatment periods; <em>Jackknife</em> leaves out one donor at a time to get a standard error; <em>None</em> skips inference.</li>",
        "<li><strong>Significance level</strong>: 0.05 gives 95% intervals.</li>",
        "<li><strong>Conformal permutations</strong>: <em>iid</em> permutes residuals at random, <em>Moving block</em> rotates them cyclically.</li>",
        "<li><strong>Random permutations</strong> and <strong>Grid size</strong>: number of random permutations (iid) and the number of effect values tried when inverting the test.</li>",
        "<li><strong>Random seed</strong>: makes the random permutations reproducible.</li></ul>",
        "<h2>Results and output options</h2><ul>",
        "<li><strong>Average effect and fit</strong>, <strong>Donor weights</strong> (weights can be negative with ridge augmentation), <strong>Covariate balance</strong>, and <strong>Effect by post-treatment period</strong> with standard errors, intervals and p-values where available.</li>",
        "<li><strong>Period-by-period path table</strong>, <strong>Method notes</strong> and <strong>Data and setup guide</strong>: optional output.</li></ul>",
        "<h2>Notes and limitations</h2><ul><li>Treated units must adopt in the same period. For staggered adoption use the Staggered Augmented SC analysis.</li><li>The conformal interval for the average effect is not defined in the package, so the table reports its p-value only.</li></ul>"),

    stagascm = paste0(
        H("Staggered-Adoption Augmented Synthetic Control", "Ben-Michael, Feller and Rothstein (2022). This analysis is a port of the multisynth function of the R package augsynth and reproduces its results. Balancing weights for all treated units are estimated jointly: they trade off the imbalance of each treated unit against the imbalance of the pooled average, after removing fixed effects."),
        paste0("<h2>Data</h2><p>A balanced long-format panel with at least two never-treated donor units.</p><ul>",
               "<li><strong>Outcome</strong>, <strong>Unit</strong>, <strong>Time</strong>: as in the other analyses.</li>",
               "<li><strong>Treatment indicator (1 = treated)</strong>: a 0/1 variable that switches from 0 to 1 at a unit's adoption and stays 1. Use this <em>or</em> the next variable.</li>",
               "<li><strong>First treated period (optional)</strong>: each unit's adoption period, blank for never-treated units.</li>",
               "<li><strong>Covariates</strong> (optional): averaged over the periods before the earliest adoption (missing values omitted), standardised to the scale of the outcomes and balanced together with them.</li></ul>"),
        "<h2>Estimation</h2><ul>",
        "<li><strong>Unit fixed effects</strong>: removes each unit's pre-treatment mean; time effects from the never-treated units are always removed. Units with a single pre-treatment period are not allowed when this is on.</li>",
        "<li><strong>Aggregate units by adoption time (time cohorts)</strong>: one set of weights per adoption period instead of one per unit. Use it when many units adopt in the same period.</li>",
        "<li><strong>Pooling parameter nu</strong>: weight on the pooled balance relative to the separate balance. Blank chooses it automatically from the separate fits; 0 fits each unit on its own; 1 balances only the average.</li>",
        "<li><strong>Ridge penalty on weights</strong>: shrinks weights towards equal weights (default 0).</li>",
        "<li><strong>Post-treatment periods</strong> and <strong>Pre-treatment periods</strong>: how many periods after adoption to estimate and how many before to balance (blank = all).</li></ul>",
        "<h2>Inference</h2><ul>",
        "<li><strong>Inference</strong>: <em>Wild bootstrap</em> reweights units with random weights and reports percentile intervals; <em>Jackknife</em> leaves out one unit at a time (it refits once per unit, so it is slow with many units) and uses a normal approximation; <em>None</em> skips it.</li>",
        "<li><strong>Bootstrap replications</strong>, <strong>Significance level</strong>, <strong>Random seed</strong>: replications, interval level (0.05 gives 95%) and reproducibility.</li></ul>",
        "<h2>Results and output options</h2><ul>",
        "<li><strong>Overall estimate and balance</strong>: the average post-treatment effect with standard error and interval, the pooling parameter, and the pooled and individual pre-treatment imbalance, also scaled by the imbalance under uniform weights (smaller than 1 means the weights help).</li>",
        "<li><strong>Effect by treated unit (or cohort)</strong> and <strong>Average effect by event time</strong> (periods since adoption), with the matching event-study plot.</li>",
        "<li><strong>Covariate imbalance</strong>: shown when covariates are selected.</li>",
        "<li><strong>Donor weights table</strong>, <strong>Method notes</strong> and <strong>Data and setup guide</strong>: optional output.</li></ul>",
        "<h2>Notes and limitations</h2><ul><li>Periods after the last adoption are not used for balancing, and units adopting in the last period count as treated there.</li><li>Matching on covariates by exact or nearest-neighbour rules, multiple outcomes and the time-regression outcome model of multisynth are not included.</li></ul>"),

    compare = paste0(
        H("Compare Estimators", "Runs several synthetic-control estimators on the same panel and reports them side by side, so that you can see whether the conclusion depends on the estimator. Report all the results, not only the most favourable one."),
        "<h2>Data</h2><p>The data must be a <strong>balanced panel in long format</strong>: one row per unit and period, with an outcome for every unit in every period.</p><ul>",
        "<li><strong>Outcome</strong>, <strong>Unit</strong>, <strong>Time</strong>: the outcome variable, the unit identifier and the numeric period.</li>",
        "<li><strong>Treated unit</strong> and <strong>First treated period</strong>: the unit that receives the treatment and the first period in which it is treated. All other units form the donor pool.</li>",
        "<li><strong>Exclude donors</strong>: comma-separated units removed from the donor pool for all methods.</li></ul>",
        "<h2>Methods to compare</h2><ul>",
        "<li><strong>Synthetic Control (classic)</strong>: Abadie, Diamond and Hainmueller (2010), outcomes of all pre-treatment periods as predictors.</li>",
        "<li><strong>Penalized</strong>: Abadie and L'Hour (2021). <strong>Demeaned</strong>: Ferman and Pinto (2021). <strong>Augmented (ridge)</strong>: Ben-Michael, Feller and Rothstein (2021), cross-validated ridge penalty, no fixed effects.</li>",
        "<li><strong>Elastic-Net</strong>: Doudchenko and Imbens (2016), alpha and penalty by cross-validation. <strong>Bayesian</strong>: Kim, Lee and Gupta (2020). <strong>Interactive Fixed Effects</strong>: Xu (2017). <strong>Matrix Completion</strong>: Athey et al. (2021). <strong>SDID</strong>: Arkhangelsky et al. (2021), SDID estimator.</li></ul>",
        "<h2>Method settings</h2><ul>",
        "<li><strong>Penalized, penalty lambda</strong>, <strong>Interactive fixed effects, number of factors</strong>, <strong>Matrix completion, nuclear-norm penalty</strong>, <strong>Bayesian, posterior draws</strong>, <strong>Synthetic control, random starts</strong> and <strong>Random seed</strong> have the same meaning as in the analysis of each method. All other settings are the defaults of that analysis.</li></ul>",
        "<h2>Inference</h2><ul><li><strong>Add each method's own uncertainty measure</strong>: a permutation p-value (post/pre RMSPE ratio over all donors) for Synthetic Control, Penalized, Demeaned, Elastic-Net, Interactive Fixed Effects and Matrix Completion; a conformal p-value for Augmented (the conformal method gives intervals per period, not for the average); a placebo standard error with a normal-approximation interval for SDID. The Bayesian estimator always reports its 95% credible interval. The measures are of different kinds and are not directly comparable. Switched off by default because it can take several minutes.</li></ul>",
        "<h2>Results and output options</h2><ul>",
        "<li><strong>Estimates by method</strong>: ATT, pre- and post-treatment RMSPE and their ratio, the number of donors with non-zero weight (blank for Matrix Completion and Interactive Fixed Effects without weights), and the uncertainty measure when requested. For SDID the RMSPE uses the counterfactual shifted by the weighted pre-treatment gap, so it is not comparable with the other rows.</li>",
        "<li><strong>Treated vs counterfactuals</strong>, <strong>Gaps by method</strong> and <strong>ATT by method</strong>: the outcome path with every counterfactual, the gaps, and the ATT of each method with its interval.</li>",
        "<li><strong>Period-by-period table</strong>, <strong>Method notes</strong> and <strong>Data and setup guide</strong>: optional output.</li></ul>",
        "<h2>Notes and limitations</h2><ul><li>Methods with a poor pre-treatment fit (large pre RMSPE) deserve less weight, but choosing the best-fitting method after seeing the effects is specification search; decide on the set of methods in advance.</li>",
        "<li>Covariates, the staggered estimator, the placebo plots and the other method-specific options are available only in the analysis of each method.</li></ul>"),

    stop("Unknown notes: ", key))
}
