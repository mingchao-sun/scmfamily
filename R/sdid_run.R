# Analysis body for the "sdid" jamovi analysis. Uses the synthdid port in sdid_engine.R.

sdid_info <- list(
    sdid = "Synthetic difference-in-differences (Algorithm 1 of Arkhangelsky et al., 2021): unit weights (omega) and time weights (lambda) are estimated by Frank-Wolfe with sparsification, and the effect is the weighted difference-in-differences.",
    sc   = "Synthetic control estimate in the synthdid framework: unit weights with an infinitesimal ridge penalty and no intercept, no time weights (the same specification as sc_estimate in the synthdid package).",
    did  = "Difference-in-differences estimate in the synthdid framework: uniform unit and time weights (the same specification as did_estimate in the synthdid package)."
)

sdid_run <- function(self, ck = NULL) {
    opts <- self$options
    notes <- self$results$notes
    if (isTRUE(opts$showGuide)) self$results$guide$setContent(fam_guide_html("sdid"))
    if (isTRUE(opts$showHelp)) self$results$helpText$setContent(fam_notes_html("sdid"))
    if (is.null(opts$outcome) || is.null(opts$unit) || is.null(opts$time) ||
        is.null(opts$treatedUnit) || !nzchar(opts$treatedUnit)) {
        notes$setContent("<p>Select an outcome, a unit identifier, a time variable and a treated unit to begin.</p>")
        return(invisible(NULL))
    }
    splitlist <- function(s) {
        v <- trimws(strsplit(if (is.null(s)) "" else s, ",", fixed = TRUE)[[1]])
        v[nzchar(v)]
    }
    reject <- function(msg) jmvcore::reject(msg)
    excl <- splitlist(opts$donorExclude)
    extra <- splitlist(opts$extraTreated)
    main <- opts$treatedUnit
    treated <- unique(c(main, extra))
    if (any(treated %in% excl)) reject("A treated unit cannot also be excluded from the donor pool.")
    covars <- as.character(fam_opt(opts, "predictors", character()))

    pan <- tryCatch(
        fam_panel(self$data, opts$outcome, opts$unit, opts$time, main, opts$treatmentTime, excl, covars),
        error = function(e) reject(conditionMessage(e)))
    bad <- setdiff(treated, colnames(pan$Y))
    if (length(bad) > 0) reject(paste0("Unknown treated unit(s): ", paste(bad, collapse = ", ")))
    controls <- setdiff(colnames(pan$Y), treated)
    N0 <- length(controls); N1 <- length(treated)
    T0 <- length(pan$pre); T1 <- length(pan$post)
    if (N0 < 2) reject("At least two control units are required.")
    ord <- c(controls, treated)
    Y <- t(pan$Y[, ord, drop = FALSE])                      # N x T, controls first
    X <- array(dim = c(dim(Y), 0))
    if (!is.null(pan$cov)) {
        Xc <- pan$cov[, ord, , drop = FALSE]
        if (anyNA(Xc)) {
            miss <- covars[vapply(seq_along(covars), function(k) anyNA(Xc[, , k]), NA)]
            reject(paste0("Covariates must be observed in every period for every unit. Missing values in: ",
                          paste(miss, collapse = ", "), "."))
        }
        X <- aperm(Xc, c(2, 1, 3))                           # N x T x K
    }
    times <- pan$times

    estimator <- fam_opt(opts, "estimator", "sdid")
    est <- tryCatch(switch(estimator,
        sdid = sdid_estimate(Y, N0, T0, X = X),
        sc   = sdid_sc_estimate(Y, N0, T0, X = X),
        did  = sdid_did_estimate(Y, N0, T0, X = X)),
        error = function(e) reject(conditionMessage(e)))
    tau <- as.numeric(est)
    w <- attr(est, "weights")

    # ---- standard error ----------------------------------------------------------------
    se <- NA_real_
    se_note <- ""
    method <- fam_opt(opts, "seMethod", "placebo")
    if (method != "none") {
        set.seed(as.integer(fam_opt(opts, "seed", 1)))
        se <- tryCatch(sdid_se(est, method, as.integer(fam_opt(opts, "seReps", 200)), checkpoint = ck),
                       error = function(e) { se_note <<- conditionMessage(e); NA_real_ })
        if (!is.finite(se) && !nzchar(se_note))
            se_note <- paste0("The ", method, " standard error is not available for this design ",
                              "(it needs more than one treated unit or more than one control with non-zero weight).")
    }

    # ---- counterfactual path and effect curve -------------------------------------------
    Xb <- sdid_contract3(X, w$beta)
    Yr <- Y - Xb
    sc_t <- as.vector(t(c(-w$omega, rep(1 / N1, N1))) %*% Yr)
    off <- sum(sc_t[seq_len(T0)] * w$lambda)
    gap <- sc_t - off
    treated_avg <- colMeans(Y[N0 + seq_len(N1), , drop = FALSE])
    path <- treated_avg - gap
    curve <- sdid_effect_curve(est)

    notes_txt <- paste0(
        "<h3>Synthetic difference-in-differences</h3><p>", sdid_info[[estimator]],
        " Control units: ", N0, "; treated units: ", N1, "; pre-treatment periods: ", T0,
        "; post-treatment periods: ", T1, ".</p>",
        if (length(covars) > 0) paste0("<p>Covariates (", paste(covars, collapse = ", "),
            ") enter as time-varying regressors with common coefficients, estimated jointly with the weights. ",
            "They must not be affected by the treatment.</p>") else "",
        if (nzchar(se_note)) paste0("<p>Standard error: ", se_note, "</p>") else "",
        "<p>Reference: Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021), American Economic Review 111(12), 4088-4118. ",
        "Estimation follows the R package synthdid.</p>")
    notes$setContent(notes_txt)

    tbl <- self$results$summary
    tbl$addRow(rowKey = "att", values = list(metric = "Average treatment effect on the treated (ATT)", value = tau))
    if (method != "none") {
        tbl$addRow(rowKey = "se", values = list(metric = paste0("Standard error (", method, ")"), value = se))
        if (is.finite(se)) {
            z <- stats::qnorm(0.975)
            tbl$addRow(rowKey = "lo", values = list(metric = "95% CI, lower (normal approximation)", value = tau - z * se))
            tbl$addRow(rowKey = "hi", values = list(metric = "95% CI, upper (normal approximation)", value = tau + z * se))
            tbl$addRow(rowKey = "p", values = list(metric = "Two-sided p-value", value = 2 * stats::pnorm(-abs(tau / se))))
        }
    }

    om <- stats::setNames(as.vector(w$omega), controls)
    for (u in names(sort(om, decreasing = TRUE))[sort(om, decreasing = TRUE) > 1e-8])
        self$results$weights$addRow(rowKey = u, values = list(unit = u, weight = om[[u]]))
    la <- as.vector(w$lambda)
    tp <- times[pan$pre]
    for (i in which(la > 1e-8))
        self$results$timeWeights$addRow(rowKey = as.character(tp[i]),
                                        values = list(time = as.character(tp[i]), weight = la[i]))
    if (!is.null(w$beta) && length(w$beta) == length(covars)) {
        for (k in seq_along(covars))
            self$results$covCoef$addRow(rowKey = covars[k],
                                        values = list(covariate = covars[k], estimate = w$beta[k]))
    }
    tq <- times[pan$post]
    for (i in seq_along(tq))
        self$results$effectCurve$addRow(rowKey = as.character(tq[i]),
                                        values = list(time = as.character(tq[i]), effect = curve[i]))
    if (isTRUE(opts$showPath)) {
        pt <- self$results$pathTable
        for (i in seq_along(times))
            pt$addRow(rowKey = as.character(times[i]),
                      values = list(time = as.character(times[i]), treated = treated_avg[i],
                                    synthetic = path[i], gap = gap[i]))
    }
    t0 <- times[pan$post[1]]
    self$results$pathsPlot$setState(list(times = times, treated = treated_avg, path = path,
                                         t0 = t0, label = opts$outcome))
    self$results$gapsPlot$setState(list(times = times, gap = gap, placebo = list(),
                                        t0 = t0, label = opts$outcome))
    invisible(NULL)
}
