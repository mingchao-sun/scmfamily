# Analysis body for the "ascm" jamovi analysis. Uses the augsynth port in augsynth_engine.R.

ascm_run <- function(self, ck = NULL) {
    opts <- self$options
    notes <- self$results$notes
    if (isTRUE(opts$showGuide)) self$results$guide$setContent(fam_guide_html("ascm"))
    if (isTRUE(opts$showHelp)) self$results$helpText$setContent(fam_notes_html("ascm"))
    if (is.null(opts$outcome) || is.null(opts$unit) || is.null(opts$time) ||
        is.null(opts$treatedUnit) || !nzchar(opts$treatedUnit)) {
        notes$setContent("<p>Select an outcome, a unit identifier, a time variable and a treated unit to begin.</p>")
        return(invisible(NULL))
    }
    reject <- function(msg) jmvcore::reject(msg)
    splitlist <- function(s) {
        v <- trimws(strsplit(if (is.null(s)) "" else s, ",", fixed = TRUE)[[1]])
        v[nzchar(v)]
    }
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
    times <- pan$times
    t_int <- opts$treatmentTime
    t0 <- length(pan$pre)
    tpost <- length(pan$post)
    t_final <- t0 + tpost

    lam <- NULL
    if (nzchar(trimws(opts$lambdaText))) {
        lam <- suppressWarnings(as.numeric(trimws(opts$lambdaText)))
        if (is.na(lam) || lam < 0) reject("The ridge penalty must be a non-negative number, or blank for cross-validation.")
    }
    progfunc <- opts$progfunc
    scm <- isTRUE(opts$scm)
    if (progfunc == "none" && !scm) reject("Without ridge augmentation the synthetic control weights are required.")

    Z <- NULL
    if (!is.null(pan$cov))
        Z <- tryCatch(asyn_covariates(pan$cov, times, t_int, colnames(pan$Y)),
                      error = function(e) reject(conditionMessage(e)))

    args <- list(Y = pan$Y, times = times, treated = treated, t_int = t_int, Z = Z,
                 progfunc = progfunc, scm = scm, fixedeff = isTRUE(opts$fixedeff))
    if (progfunc == "ridge") {
        args$min_1se <- isTRUE(opts$min1se)
        if (!is.null(lam)) args$lambda <- lam
        if (!is.null(Z)) args$residualize <- isTRUE(opts$residualize)
    }
    fit <- tryCatch(do.call(asyn_fit, args), error = function(e) reject(conditionMessage(e)))

    att <- asyn_predict(fit, att = TRUE)
    y_tr <- colMeans(cbind(fit$data$X, fit$data$y)[fit$data$trt == 1, , drop = FALSE])
    path <- asyn_predict(fit, att = FALSE)
    post_idx <- (t0 + 1):t_final
    att_avg <- mean(att[post_idx])

    # ---- inference ---------------------------------------------------------------------
    inf <- fam_opt(opts, "infType", "conformal")
    alpha <- fam_opt(opts, "alpha", 0.05)
    set.seed(as.integer(fam_opt(opts, "seed", 1)))
    se_vec <- lb <- ub <- pv <- rep(NA_real_, t_final + 1)
    inf_err <- ""
    if (inf != "none") {
        res <- tryCatch(switch(inf,
            conformal = asyn_conformal_inf(fit, alpha = alpha, type = opts$conformalType,
                                           ns = as.integer(opts$nPerm), grid_size = as.integer(opts$gridSize),
                                           checkpoint = ck),
            "jackknife+" = asyn_time_jackknife_plus(fit, alpha = alpha, checkpoint = ck),
            jackknife = asyn_jackknife_se(fit, checkpoint = ck)),
            error = function(e) { inf_err <<- conditionMessage(e); NULL })
        if (!is.null(res)) {
            if (!is.null(res$se)) se_vec <- res$se
            if (!is.null(res$lb)) lb <- res$lb
            if (!is.null(res$ub)) ub <- res$ub
            if (!is.null(res$p_val)) pv <- res$p_val
        }
    }
    nm_inf <- c(conformal = "conformal", "jackknife+" = "time jackknife+", jackknife = "jackknife", none = "none")[[inf]]

    notes_txt <- paste0(
        "<h3>Augmented synthetic control</h3><p>",
        if (progfunc == "ridge") "Simplex synthetic control weights (fitted by quadratic programming with osqp) are bias-corrected by a ridge regression of post-treatment outcomes on pre-treatment outcomes"
        else "Simplex synthetic control weights fitted by quadratic programming with osqp, with no outcome-model correction",
        if (progfunc == "ridge") sprintf("; the ridge penalty is %s (lambda = %s).",
            if (is.null(lam)) "chosen by leave-out cross-validation over pre-treatment periods" else "fixed by the user",
            format(signif(fit$lambda, 6))) else ".",
        if (!scm) " Uniform weights replace the synthetic control weights." else "",
        if (isTRUE(opts$fixedeff)) " Unit fixed effects: outcomes are demeaned by their pre-treatment means." else "",
        " Donors: ", sum(fit$data$trt == 0), "; treated units: ", sum(fit$data$trt == 1),
        "; pre-treatment periods: ", t0, "; post-treatment periods: ", tpost, ".</p>",
        if (length(covars) > 0) paste0("<p>Covariates (", paste(covars, collapse = ", "),
            "): pre-treatment means, ",
            if (isTRUE(opts$residualize) && progfunc == "ridge") "regressed out of the pre-treatment outcomes (residualized)."
            else "standardised to the scale of the outcomes and balanced together with them.", "</p>") else "",
        if (nzchar(inf_err)) paste0("<p>Inference failed: ", inf_err, "</p>") else "",
        "<p>Reference: Ben-Michael, Feller and Rothstein (2021), Journal of the American Statistical Association 116(536), 1789-1803. ",
        "Estimation follows the R package augsynth.</p>")
    notes$setContent(notes_txt)

    tbl <- self$results$summary
    tbl$addRow(rowKey = "att", values = list(metric = "Average post-treatment effect (ATT)", value = att_avg))
    if (inf != "none" && !nzchar(inf_err)) {
        if (inf == "jackknife")
            tbl$addRow(rowKey = "se", values = list(metric = "Standard error (jackknife)", value = se_vec[t_final + 1]))
        if (inf %in% c("conformal", "jackknife+") && is.finite(lb[t_final + 1]) && is.finite(ub[t_final + 1])) {
            lvl <- sprintf("%g%%", 100 * (1 - alpha))
            tbl$addRow(rowKey = "lo", values = list(metric = paste0(lvl, " interval, lower (", nm_inf, ")"), value = lb[t_final + 1]))
            tbl$addRow(rowKey = "hi", values = list(metric = paste0(lvl, " interval, upper (", nm_inf, ")"), value = ub[t_final + 1]))
        }
        if (inf == "conformal")
            tbl$addRow(rowKey = "p", values = list(metric = "p-value, average effect = 0 (conformal)", value = pv[t_final + 1]))
    }
    pre_gap <- att[seq_len(t0)]
    tbl$addRow(rowKey = "pre", values = list(metric = "RMSPE, pre-treatment", value = sqrt(mean(pre_gap^2))))
    tbl$addRow(rowKey = "post", values = list(metric = "RMSPE, post-treatment", value = sqrt(mean(att[post_idx]^2))))
    if (progfunc == "ridge")
        tbl$addRow(rowKey = "lambda", values = list(metric = "Ridge penalty (lambda)", value = fit$lambda))
    tbl$addRow(rowKey = "l2", values = list(metric = "Pre-treatment L2 imbalance", value = fit$l2_imbalance))
    tbl$addRow(rowKey = "sl2", values = list(metric = "L2 imbalance, scaled by uniform weights", value = fit$scaled_l2_imbalance))

    w <- fit$weights[, 1]
    for (u in names(sort(w, decreasing = TRUE))[sort(w, decreasing = TRUE) != 0 & abs(sort(w, decreasing = TRUE)) > 1e-8])
        self$results$weights$addRow(rowKey = u, values = list(unit = u, weight = w[[u]]))

    if (!is.null(Z)) {
        tr <- fit$data$trt == 1
        zt <- colMeans(Z[tr, , drop = FALSE])
        zs <- as.vector(t(Z[!tr, , drop = FALSE]) %*% fit$weights)
        for (k in seq_along(covars))
            self$results$balance$addRow(rowKey = covars[k],
                values = list(predictor = covars[k], treated = zt[k], synthetic = zs[k]))
    }

    tq <- times[pan$post]
    for (j in seq_along(tq)) {
        i <- t0 + j
        self$results$effects$addRow(rowKey = as.character(tq[j]), values = list(
            time = as.character(tq[j]), estimate = att[i], se = se_vec[i], lower = lb[i], upper = ub[i], p = pv[i]))
    }
    if (isTRUE(opts$showPath)) {
        pt <- self$results$pathTable
        for (i in seq_along(times))
            pt$addRow(rowKey = as.character(times[i]),
                      values = list(time = as.character(times[i]), treated = y_tr[i],
                                    synthetic = path[i], gap = att[i]))
    }
    t0_time <- times[pan$post[1]]
    self$results$pathsPlot$setState(list(times = times, treated = unname(y_tr), path = unname(path),
                                         t0 = t0_time, label = opts$outcome))
    self$results$gapsPlot$setState(list(times = times, gap = unname(att), placebo = list(),
                                        t0 = t0_time, label = opts$outcome))
    invisible(NULL)
}
