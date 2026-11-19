# Analysis body for the "stagascm" jamovi analysis (multisynth port in multisynth_engine.R).

stagascm_draw_event <- function(st) {
    if (is.null(st)) return(FALSE)
    ok <- !is.na(st$event) & !is.na(st$est)
    x <- st$event[ok]
    yl <- range(c(st$est[ok], st$lower[ok], st$upper[ok]), na.rm = TRUE)
    graphics::plot(x, st$est[ok], type = "n", ylim = yl, xlab = "Periods since adoption",
                   ylab = paste("Effect on", st$label), main = "Average effect by event time")
    if (any(!is.na(st$lower[ok])))
        graphics::polygon(c(x, rev(x)), c(st$lower[ok], rev(st$upper[ok])),
                          col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
    graphics::lines(x, st$est[ok], lwd = 2)
    graphics::abline(h = 0, lty = 3)
    graphics::abline(v = -0.5, lty = 3)
    TRUE
}

stagascm_run <- function(self, ck = NULL) {
    opts <- self$options
    notes <- self$results$notes
    if (isTRUE(opts$showGuide)) self$results$guide$setContent(fam_guide_html("stagascm"))
    if (isTRUE(opts$showHelp)) self$results$helpText$setContent(fam_notes_html("stagascm"))
    reject <- function(msg) jmvcore::reject(msg)
    if (is.null(opts$outcome) || is.null(opts$unit) || is.null(opts$time) ||
        (is.null(opts$treatInd) && is.null(opts$adoption))) {
        notes$setContent("<p>Select an outcome, unit, time, and either a treatment indicator or a first-treated-period variable to begin.</p>")
        return(invisible(NULL))
    }
    data <- self$data
    covars <- as.character(fam_opt(opts, "predictors", character()))

    # ---- build the panel ---------------------------------------------------------------
    tt <- suppressWarnings(as.numeric(as.character(data[[opts$time]])))
    if (anyNA(tt)) reject("The time variable must be numeric (for example a year).")
    d <- data.frame(u = as.character(data[[opts$unit]]), t = tt, y = as.numeric(data[[opts$outcome]]),
                    stringsAsFactors = FALSE)
    if (!is.null(opts$treatInd))
        d$g <- as.character(data[[opts$treatInd]]) %in% c("1", "TRUE", "true", "yes", "Yes")
    else d$g <- FALSE
    d$a <- if (!is.null(opts$adoption)) suppressWarnings(as.numeric(as.character(data[[opts$adoption]]))) else NA_real_
    cm <- matrix(NA_real_, nrow(data), length(covars), dimnames = list(NULL, covars))
    for (k in covars) cm[, k] <- suppressWarnings(as.numeric(as.character(data[[k]])))
    cc <- stats::complete.cases(d[, c("u", "t", "y")])
    d <- d[cc, , drop = FALSE]
    cm <- cm[cc, , drop = FALSE]
    if (anyDuplicated(d[, c("u", "t")])) reject("Each unit must have one row per period (duplicate unit-time rows found).")
    times <- sort(unique(d$t))
    units <- sort(unique(d$u))
    Y <- matrix(NA_real_, length(times), length(units), dimnames = list(NULL, units))
    ix <- cbind(match(d$t, times), match(d$u, units))
    Y[ix] <- d$y
    if (anyNA(Y)) reject("The panel is unbalanced: every unit needs a value in every period.")
    cov <- NULL
    if (length(covars) > 0) {
        cov <- array(NA_real_, dim = c(length(times), length(units), length(covars)),
                     dimnames = list(NULL, units, covars))
        for (k in seq_along(covars)) cov[cbind(ix, k)] <- cm[, k]
    }
    adopt <- vapply(units, function(u) {
        if (!is.null(opts$treatInd)) {
            x <- d$t[d$u == u & d$g]
            if (length(x)) min(x) else NA_real_
        } else {
            x <- d$a[d$u == u]
            x <- x[!is.na(x)]
            if (length(x)) x[1] else NA_real_
        }
    }, 0)
    names(adopt) <- units
    adopt[!is.na(adopt) & adopt > max(times)] <- NA
    if (all(is.na(adopt))) reject("No treated units were found (adoption time inside the sample).")
    if (any(!is.na(adopt) & !(adopt %in% times)))
        reject("Every adoption time must be one of the observed time periods.")
    early <- units[!is.na(adopt) & adopt <= times[1]]
    if (length(early))
        reject(paste0("These units are treated in the first period and have no pre-treatment data: ",
                      paste(early, collapse = ", "),
                      ". If the treatment indicator is already 1 in the first periods, use an indicator that switches from 0 to 1 at adoption."))

    # ---- options -----------------------------------------------------------------------
    num_or_null <- function(s, what, lo = -Inf) {
        s <- trimws(s)
        if (!nzchar(s)) return(NULL)
        v <- suppressWarnings(as.numeric(s))
        if (is.na(v) || v < lo) reject(paste0(what, " must be a number", if (is.finite(lo)) paste0(" >= ", lo) else "", ", or blank."))
        v
    }
    nu <- num_or_null(opts$nuText, "The pooling parameter nu", 0)
    if (!is.null(nu) && nu > 1) reject("The pooling parameter nu must lie between 0 and 1.")
    n_leads <- num_or_null(opts$nLeads, "The number of post-treatment periods", 1)
    n_lags <- num_or_null(opts$nLags, "The number of pre-treatment periods", 1)
    fixedeff <- isTRUE(opts$fixedeff)
    tcoh <- isTRUE(opts$timeCohort)

    Z <- NULL
    if (!is.null(cov)) Z <- tryCatch(ms_covariates(cov, times, adopt, sort(units)),
                                     error = function(e) reject(conditionMessage(e)))
    fit <- tryCatch(ms_fit(Y, times, adopt, Z = Z, n_leads = n_leads, n_lags = n_lags, nu = nu,
                           lambda = opts$lambda, fixedeff = fixedeff, time_cohort = tcoh),
                    error = function(e) reject(conditionMessage(e)))

    inf <- fam_opt(opts, "infType", "bootstrap")
    alpha <- fam_opt(opts, "alpha", 0.05)
    set.seed(as.integer(fam_opt(opts, "seed", 1)))
    inf_err <- ""
    sm <- tryCatch(suppressWarnings(
        if (inf == "bootstrap") ms_summary(fit, "bootstrap", n_boot = as.integer(opts$nBoot), alpha = alpha, checkpoint = ck)
        else if (inf == "jackknife") ms_summary(fit, "jackknife", alpha = alpha, checkpoint = ck)
        else ms_summary(fit, "none")),
        error = function(e) { inf_err <<- conditionMessage(e); NULL })
    if (is.null(sm)) {
        sm <- ms_summary(fit, "none")
    }

    nr <- nrow(sm$att)                       # last row = post-treatment average
    treated_names <- sm$groups[-1]
    ntr <- length(treated_names)
    levl <- sprintf("%g%%", 100 * (1 - alpha))
    notes$setContent(paste0(
        "<h3>Staggered-adoption augmented synthetic control (multisynth)</h3><p>",
        "Balancing weights for each treated ", if (tcoh) "adoption cohort" else "unit",
        " are estimated jointly by quadratic programming (osqp). They minimise a mix, set by the pooling parameter nu (",
        format(signif(fit$nu, 4)), if (is.null(nu)) ", chosen automatically from the separate fits" else ", user-set",
        "), of the imbalance of each treated ", if (tcoh) "cohort" else "unit",
        " and of the pooled average, in the residuals after removing ",
        if (fixedeff) "time effects (from never-treated units) and unit pre-treatment means"
        else "time effects (from never-treated units)", ". Effects are estimated up to ",
        fit$n_leads, " periods after adoption, balancing up to ", fit$n_lags, " pre-treatment periods. ",
        "Never-treated units: ", sum(!is.finite(fit$data$trt)), "; treated units: ", sum(is.finite(fit$data$trt)), ".</p>",
        if (length(covars) > 0) paste0("<p>Covariates (", paste(covars, collapse = ", "),
            "): means over the periods before the earliest adoption, balanced together with the outcomes.</p>") else "",
        if (inf == "bootstrap") paste0("<p>Inference: wild bootstrap (", opts$nBoot, " replications); ", levl, " intervals use bootstrap percentiles.</p>")
        else if (inf == "jackknife") paste0("<p>Inference: jackknife (leave one unit out); ", levl, " intervals use the normal approximation.</p>")
        else "",
        if (nzchar(inf_err)) paste0("<p>Inference failed: ", inf_err, if (inf == "jackknife") " Try the wild bootstrap instead." else "", "</p>") else "",
        "<p>Reference: Ben-Michael, Feller and Rothstein (2022), Journal of the Royal Statistical Society Series B 84(2), 351-381. ",
        "Estimation follows the R package augsynth (multisynth).</p>"))

    tbl <- self$results$summary
    tbl$addRow(rowKey = "att", values = list(metric = "Average post-treatment effect (ATT)", value = sm$att[nr, 1]))
    if (inf != "none" && !nzchar(inf_err)) {
        tbl$addRow(rowKey = "se", values = list(metric = paste0("Standard error (", inf, ")"), value = sm$se[nr, 1]))
        tbl$addRow(rowKey = "lo", values = list(metric = paste0(levl, " interval, lower"), value = sm$lower[nr, 1]))
        tbl$addRow(rowKey = "hi", values = list(metric = paste0(levl, " interval, upper"), value = sm$upper[nr, 1]))
    }
    tbl$addRow(rowKey = "ntr", values = list(metric = "Treated units", value = sum(is.finite(fit$data$trt))))
    tbl$addRow(rowKey = "nctl", values = list(metric = "Never-treated donor units", value = sum(!is.finite(fit$data$trt))))
    tbl$addRow(rowKey = "nu", values = list(metric = "Pooling parameter nu", value = fit$nu))
    tbl$addRow(rowKey = "gl2", values = list(metric = "Pooled pre-treatment imbalance (global L2)", value = fit$global_l2))
    tbl$addRow(rowKey = "sgl2", values = list(metric = "Global L2, scaled by uniform weights", value = fit$scaled_global_l2))
    tbl$addRow(rowKey = "il2", values = list(metric = "Individual imbalance (L2)", value = fit$ind_l2))
    tbl$addRow(rowKey = "sil2", values = list(metric = "Individual L2, scaled by uniform weights", value = fit$scaled_ind_l2))

    for (j in seq_along(treated_names))
        self$results$units$addRow(rowKey = treated_names[j], values = list(
            unit = treated_names[j], att = sm$att[nr, j + 1], se = sm$se[nr, j + 1],
            lower = sm$lower[nr, j + 1], upper = sm$upper[nr, j + 1]))

    for (i in seq_len(nr - 1))
        self$results$event$addRow(rowKey = as.character(sm$event[i]), values = list(
            k = as.integer(sm$event[i]), estimate = sm$att[i, 1], se = sm$se[i, 1],
            lower = sm$lower[i, 1], upper = sm$upper[i, 1]))

    if (!is.null(fit$imbalance_aux)) {
        for (k in rownames(fit$imbalance_aux))
            self$results$balance$addRow(rowKey = k, values = list(covariate = k, average = fit$imbalance_aux[k, 1]))
    }
    if (isTRUE(opts$showWeights)) {
        W <- fit$weights
        gl <- treated_names
        for (j in seq_along(gl)) {
            w <- W[, j]
            nz <- which(abs(w) > 1e-6)
            for (i in nz[order(-w[nz])])
                self$results$weights$addRow(rowKey = paste(gl[j], rownames(W)[i], sep = "|"),
                    values = list(treated = gl[j], donor = rownames(W)[i], weight = w[i]))
        }
    }
    self$results$eventPlot$setState(list(event = sm$event[-nr], est = sm$att[-nr, 1],
                                         lower = sm$lower[-nr, 1], upper = sm$upper[-nr, 1],
                                         label = opts$outcome))
    invisible(NULL)
}
