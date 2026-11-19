# "Data and setup guide": short usage notes shown at the top of each analysis when the
# option "Data and setup guide" is ticked. Each guide names the sample dataset to open,
# the exact variable settings for it, and the data layout the analysis expects.

fam_guide_table <- function(rows) {
    paste0("<table style='border-collapse:collapse;margin:4px 0 8px 0'>",
           "<tr><th style='text-align:left;padding:2px 14px 2px 0;border-bottom:1px solid #999'>Field</th>",
           "<th style='text-align:left;padding:2px 0;border-bottom:1px solid #999'>Set to</th></tr>",
           paste0("<tr><td style='padding:2px 14px 2px 0'>", names(rows), "</td><td style='padding:2px 0'>",
                  unlist(rows), "</td></tr>", collapse = ""),
           "</table>")
}

fam_guide_layout <- function(extra = NULL) {
    paste0(
        "<h4>Data layout</h4>",
        "<p>One row per unit and period (long format). Illustrative example:</p>",
        "<table style='border-collapse:collapse;margin:4px 0 8px 0'>",
        "<tr><th style='padding:2px 14px 2px 0;text-align:left;border-bottom:1px solid #999'>unit</th>",
        "<th style='padding:2px 14px 2px 0;text-align:left;border-bottom:1px solid #999'>year</th>",
        "<th style='padding:2px 14px 2px 0;text-align:left;border-bottom:1px solid #999'>outcome</th>",
        "<th style='padding:2px 0;text-align:left;border-bottom:1px solid #999'>covariate</th></tr>",
        "<tr><td>A</td><td>1988</td><td>90.1</td><td>5.2</td></tr>",
        "<tr><td>A</td><td>1989</td><td>88.7</td><td>5.3</td></tr>",
        "<tr><td>B</td><td>1988</td><td>101.4</td><td>6.0</td></tr>",
        "<tr><td>B</td><td>1989</td><td>102.2</td><td>6.1</td></tr></table>",
        "<ul>",
        "<li><b>Balanced panel</b>: every unit is observed in every period, with no duplicate unit-period rows and no missing outcomes.</li>",
        "<li><b>Variable types in jamovi</b>: set the unit variable to Nominal (or ID) and the time variable to Continuous (a numeric year or period index); the outcome and covariates are Continuous.</li>",
        "<li>The first treated period must be one of the observed time values, and there must be at least two pre-treatment periods and two never-treated donor units.</li>",
        extra,
        "</ul>")
}

fam_guide_sample <- list(
    prop99 = list(
        "Dataset" = "<b>California Proposition 99</b>",
        "Outcome" = "cigsale", "Unit" = "state", "Time" = "year",
        "Treated unit" = "California", "First treated period" = "1989"),
    prop99cov = list(
        "Dataset" = "<b>California Proposition 99 (with covariates)</b>",
        "Outcome" = "cigsale", "Unit" = "state", "Time" = "year",
        "Treated unit" = "California", "First treated period" = "1989",
        "Covariates" = "lnincome, beer, age15to24, retprice"),
    germany = list(
        "Dataset" = "<b>West Germany reunification</b>",
        "Outcome" = "gdp", "Unit" = "country", "Time" = "year",
        "Treated unit" = "West Germany", "First treated period" = "1990",
        "Covariates" = "infrate, trade, schooling, industry (some years are missing, which is fine here)")
)

fam_guide_raw <- function(key) {
    S <- fam_guide_sample
    head <- function(title, intro) paste0("<h3>", title, ": data and setup</h3><p>", intro, "</p>")
    sec <- function(title, body) paste0("<h4>", title, "</h4>", body)
    out <- switch(key,

    scm = paste0(
        head("Synthetic Control", "Classic synthetic control (Abadie, Diamond and Hainmueller, 2010) for one treated unit and a pool of untreated donors."),
        sec("Sample data", paste0(
            fam_guide_table(S$prop99),
            "<p>Leave <i>Predictors</i> empty to match on pre-treatment outcomes, or open the dataset <b>with covariates</b> and add lnincome, beer, age15to24 and retprice: their pre-treatment means (or median / last value, see <i>Predictor summary</i>) are then matched.</p>",
            "<p>A second example is <b>West Germany reunification</b>: outcome gdp, unit country, time year, treated unit West Germany, first treated period 1990; covariates such as trade and industry can be added (their missing years are skipped).</p>")),
        fam_guide_layout(),
        sec("Options worth knowing", paste0(
            "<ul><li><b>Exclude donors</b>: comma-separated units removed from the donor pool (for example states with similar policies).</li>",
            "<li><b>In-space placebo test</b> reruns the method on every donor to get permutation p-values; the robustness checks (leave-one-out, treatment timing, pre-period length, in-time placebo) need no extra data.</li></ul>")))
    ,

    penalized = paste0(
        head("Penalized Synthetic Control", "Abadie and L'Hour (2021): simplex weights with a penalty that favours donors close to the treated unit."),
        sec("Sample data", paste0(
            fam_guide_table(S$prop99),
            "<p>With covariates, add lnincome, beer, age15to24 and retprice. Another example is <b>West Germany reunification</b>: ",
            "outcome gdp, unit country, time year, treated unit West Germany, first treated period 1990, covariates trade and industry.</p>")),
        fam_guide_layout("<li><b>Covariates</b> enter through their pre-treatment means (missing years are skipped), standardised and balanced with the outcomes; <i>Covariate weight</i> sets the weight of the whole block (0 = ignore, 1 = equal to the outcomes).</li>"),
        sec("Options worth knowing", "<ul><li><b>Penalty lambda</b>: 0 gives the plain synthetic control; larger values favour similar donors. Try a few values.</li><li><b>Permutation inference</b> reruns the method on each donor; the pre-fit filter drops donors that fit badly before the treatment.</li></ul>"))
    ,

    demeaned = paste0(
        head("Demeaned Synthetic Control", "Ferman and Pinto (2021): simplex weights fitted to outcomes after removing each unit's pre-treatment mean, so donors can match trends when levels differ."),
        sec("Sample data", paste0(
            fam_guide_table(S$prop99),
            "<p>Optionally open the <b>with covariates</b> version and add lnincome, beer, age15to24 and retprice (pre-treatment means, see <i>Covariate weight</i>).</p>",
            "<p>A second example is <b>West Germany reunification</b>: outcome gdp, unit country, time year, treated unit West Germany, first treated period 1990; covariates such as trade and industry can be added (their missing years are skipped).</p>")),
        fam_guide_layout(),
        sec("Options worth knowing", "<ul><li>Use this version when the treated unit has a very different level from the donors but a similar trend.</li><li><b>Permutation inference</b> and the optional conformal inference are available.</li></ul>"))
    ,

    ascm = paste0(
        head("Augmented Synthetic Control", "Ben-Michael, Feller and Rothstein (2021): synthetic control weights plus a ridge regression correction. This analysis follows the R package <i>augsynth</i>."),
        sec("Sample data", paste0(
            fam_guide_table(S$prop99),
            "<p>For covariates open the <b>with covariates</b> dataset, or <b>West Germany reunification</b> (settings: gdp, country, year, West Germany, 1990; covariates infrate, trade, schooling, industry).</p>")),
        fam_guide_layout("<li><b>Covariates</b> are averaged over the pre-treatment periods (missing values omitted) and balanced together with the outcomes, or residualized if you tick that option.</li><li><b>Additional treated units</b> (comma-separated) are averaged with the main treated unit; they must be treated in the same period.</li>"),
        sec("Options worth knowing", "<ul><li><b>Ridge penalty</b>: leave blank for cross-validation.</li><li><b>Inference</b>: conformal (default), time jackknife+ or leave-one-unit-out jackknife.</li><li><b>Unit fixed effects</b> demeans each unit by its pre-treatment mean.</li></ul>"))
    ,

    enet = paste0(
        head("Elastic-Net Synthetic Control", "Doudchenko and Imbens (2016): unrestricted weights (may be negative, need not sum to one) with an intercept, fitted by elastic net."),
        sec("Sample data", paste0(fam_guide_table(S$prop99),
            "<p>A second example is <b>West Germany reunification</b>: outcome gdp, unit country, time year, treated unit West Germany, first treated period 1990. This analysis uses the outcome only.</p>")),
        fam_guide_layout("<li>This analysis uses the outcome only; it has no covariate option.</li>"),
        sec("Options worth knowing", "<ul><li><b>Choose alpha by cross-validation</b> (on by default) picks alpha and the penalty jointly, as in the paper; switch it off to fix <b>Mixing alpha</b> (0 = ridge, 1 = lasso). The selected alpha, penalty and intercept are reported.</li><li>With few pre-treatment periods the estimate can be unstable; compare with the simplex-weight methods.</li></ul>"))
    ,

    bscm = paste0(
        head("Bayesian Synthetic Control", "Kim, Lee and Gupta (2020): horseshoe-shrinkage regression of the treated outcome on the donors; the ATT comes with a 95% credible interval."),
        sec("Sample data", paste0(fam_guide_table(S$prop99),
            "<p>A second example is <b>West Germany reunification</b>: outcome gdp, unit country, time year, treated unit West Germany, first treated period 1990. This analysis uses the outcome only.</p>")),
        fam_guide_layout("<li>This analysis uses the outcome only; it has no covariate option.</li>"),
        sec("Options worth knowing", "<ul><li><b>Posterior draws</b>: more draws give smoother intervals but take longer.</li><li>The weights table shows posterior-mean coefficients; they need not be non-negative or sum to one.</li></ul>"))
    ,

    ife = paste0(
        head("Interactive Fixed Effects (GSC)", "Xu (2017): two-way fixed effects and latent factors estimated from the donor pool, with optional time-varying covariates."),
        sec("Sample data", paste0(fam_guide_table(S$prop99),
            "<p>Set <i>Number of factors (r)</i> to 1 or 2 and compare several values.</p>",
            "<p>A second example is <b>West Germany reunification</b>: gdp, country, year, West Germany, 1990, with no covariates (its covariates have missing years, so they cannot be used here).</p>")),
        fam_guide_layout("<li><b>Covariates</b> are time-varying regressors with one coefficient for all units. They must be observed in <b>every period for every unit</b> and must not be affected by the treatment. The sample datasets have no covariate that satisfies this (the gaps in income, beer and age, and the price response of retprice), so run the sample without covariates and use your own panel for covariates.</li>"),
        sec("Options worth knowing", "<ul><li>The number of factors cannot exceed the number of pre-treatment periods minus 2.</li><li>Conformal inference is available.</li></ul>"))
    ,

    mc = paste0(
        head("Matrix Completion", "Athey, Bayati, Doudchenko, Imbens and Khosravi (2021): the treated unit's post-treatment outcomes are imputed from a low-rank matrix plus unit and time effects."),
        sec("Sample data", paste0(fam_guide_table(S$prop99),
            "<p>A second example is <b>West Germany reunification</b>: outcome gdp, unit country, time year, treated unit West Germany, first treated period 1990. This analysis uses the outcome only.</p>")),
        fam_guide_layout("<li>This analysis uses the outcome only; it has no covariate option.</li>"),
        sec("Options worth knowing", "<ul><li><b>Nuclear-norm penalty</b> is relative to the largest singular value; try several values.</li><li>No donor weights are reported because the estimator is not a weighted average of donors.</li></ul>"))
    ,

    sdid = paste0(
        head("Synthetic Difference-in-Differences", "Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021). This analysis is a port of the R package <i>synthdid</i> and reproduces it."),
        sec("Sample data", paste0(fam_guide_table(S$prop99),
            "<p>Another example is <b>West Germany reunification</b>: gdp, country, year, West Germany, 1990.</p>")),
        fam_guide_layout("<li><b>Additional treated units</b>: units treated in the same first period, comma-separated. With one treated unit only the placebo standard error is available; the bootstrap and jackknife need several treated units.</li><li><b>Covariates</b> must be observed in every period for every unit and must not respond to the treatment. The sample datasets have none that qualify, so run them without covariates.</li>"),
        sec("Options worth knowing", "<ul><li><b>Estimator</b>: SDID, synthetic control or difference-in-differences within the same framework, for comparison.</li><li>Unit weights (omega) and time weights (lambda) are reported, with the effect curve for the post-treatment periods.</li><li>Standard errors can be slow (placebo with 200 replications takes about half a minute on the sample, much longer with covariates): try 20 replications first.</li></ul>"))
    ,

    stagascm = paste0(
        head("Staggered Augmented Synthetic Control", "Ben-Michael, Feller and Rothstein (2022): several units adopt the treatment at different times. This analysis follows the <i>multisynth</i> function of the R package <i>augsynth</i>."),
        sec("Sample data", paste0(
            fam_guide_table(list("Dataset" = "<b>Minimum wage and teen employment (staggered)</b>", "Outcome" = "lemp", "Unit" = "county", "Time" = "year",
                                 "First treated period" = "first_treat (blank for never-treated counties)",
                                 "Treatment indicator" = "leave empty (or use <i>treated</i> instead of first_treat)",
                                 "Covariates (optional)" = "lpop",
                                 "Unit fixed effects" = "<b>off</b> (2004 counties have a single pre-treatment year)",
                                 "Aggregate units by adoption time" = "<b>on</b> (otherwise 191 treated counties each get their own weights)",
                                 "Inference" = "wild bootstrap")))),
        sec("Data layout", paste0(
            "<p>Long format with one row per unit and period, as above, plus <b>either</b> a variable with each unit's first treated period (blank for never-treated units) <b>or</b> a 0/1 treatment indicator that switches from 0 to 1 at adoption and stays 1.</p><ul>",
            "<li>The panel must be balanced, with at least two never-treated donor units.</li>",
            "<li>Adoption times must be observed periods, and no unit may be treated in the first period.</li>",
            "<li>Periods after the last adoption are dropped from the pre-treatment window; units adopting in that last period count as treated in the final period.</li>",
            "<li>Covariates are averaged over the periods before the earliest adoption.</li></ul>")),
        sec("Options worth knowing", "<ul><li><b>Pooling parameter nu</b>: blank chooses it automatically; 0 fits each unit separately and 1 only balances the pooled average.</li><li>With many treated units, prefer the wild bootstrap; the jackknife refits the model once per unit.</li></ul>"))
    ,

    compare = paste0(
        head("Compare Estimators", "Runs several of the module's estimators on the same panel and shows their estimates, fit and counterfactual paths side by side. It uses pre-treatment outcomes only."),
        sec("Sample data", paste0(fam_guide_table(S$prop99),
            "<p>A second example is <b>West Germany reunification</b>: outcome gdp, unit country, time year, treated unit West Germany, first treated period 1990.</p>",
            "<p>Under <i>Methods to compare</i> tick the estimators to include. Synthetic Control, Demeaned, Augmented and SDID are ticked by default.</p>")),
        fam_guide_layout("<li>This analysis has no covariate option: every method uses the pre-treatment outcomes only. For covariates use the analysis of the method itself.</li>"),
        sec("Options worth knowing", paste0("<ul><li><b>Method settings</b> hold the one setting per method that matters most (penalty, number of factors, and so on); everything else follows the defaults of the method's own analysis.</li>",
            "<li><b>Add each method's own uncertainty measure</b> adds a permutation p-value, a conformal interval, a placebo standard error or a credible interval, depending on the method. This can take several minutes.</li>",
            "<li>The estimates equal those of the stand-alone analyses run with their defaults.</li></ul>"))),

    stop("Unknown guide: ", key))
    out
}

# Adds where to open the example data after the introduction.
fam_guide_html <- function(key) {
    x <- fam_guide_raw(key)
    sub("</p>", paste0("</p><p><i>Open the sample data from File &gt; Open &gt; Data Library &gt; Synthetic Control Family; ",
                       "the library lists each data set by the bold name below.</i></p>"), x, fixed = TRUE)
}
