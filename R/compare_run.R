# Analysis body for the "compare" jamovi analysis: runs several of the module's estimators on
# the same panel (outcomes only, no covariates) and puts the results side by side.
# Each estimate is computed by the same code as the stand-alone analysis with its default
# settings (see tests/test_methods.R, which checks that the numbers are identical).

compare_methods <- list(
    scm  = list(opt = "cmpScm",  label = "Synthetic Control",          col = "#0072B2"),
    pen  = list(opt = "cmpPen",  label = "Penalized",                  col = "#E69F00"),
    dem  = list(opt = "cmpDem",  label = "Demeaned",                   col = "#009E73"),
    ascm = list(opt = "cmpAscm", label = "Augmented (ridge)",          col = "#D55E00"),
    enet = list(opt = "cmpEnet", label = "Elastic-Net",                col = "#CC79A7"),
    bscm = list(opt = "cmpBscm", label = "Bayesian",                   col = "#56B4E9"),
    ife  = list(opt = "cmpIfe",  label = "Interactive Fixed Effects",  col = "#6A3D9A"),
    mc   = list(opt = "cmpMc",   label = "Matrix Completion",          col = "#8C510A"),
    sdid = list(opt = "cmpSdid", label = "SDID",                       col = "#666666")
)

# One estimator. Returns path (counterfactual, one value per period), gap, ATT, fit summary and,
# when requested, the method's own uncertainty measure.
compare_one <- function(key, ctx, opts, inference, ck) {
    Y <- ctx$Y; tr <- ctx$tr; dn <- ctx$dn; pre <- ctx$pre; post <- ctx$post
    out <- list(key = key, se = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_, basis = "")
    summarise <- function(path, att = NULL) {
        gap <- as.vector(Y[, tr]) - path
        out$path <<- path
        out$gap <<- gap
        out$att <<- if (is.null(att)) mean(gap[post]) else att
        out$rmspe_pre <<- sqrt(mean(gap[pre]^2))
        out$rmspe_post <<- sqrt(mean(gap[post]^2))
        out$ratio <<- out$rmspe_post / out$rmspe_pre
    }
    if (key == "scm") {
        fit <- scm_fit(Y, treat = tr, pre = pre, post = post, X = NULL,
                       starts = opts$optStarts, seed = opts$seed)
        summarise(as.vector(fit$synth)); out$nw <- length(fit$weights)
        if (inference) {
            pl <- scm_placebo(Y, treat = tr, pre = pre, post = post, X = NULL, fit_true = fit,
                              checkpoint = ck, starts = opts$optStarts, seed = opts$seed)
            out$p <- pl$p_value
            out$basis <- "permutation, post/pre RMSPE ratio (all donors)"
        }
    } else if (key %in% c("pen", "dem", "enet", "bscm", "ife", "mc")) {
        method <- c(pen = "penalized", dem = "demeaned", enet = "enet", bscm = "bscm", ife = "ife", mc = "mc")[[key]]
        fitf <- function(tr_, dn_)
            fam_fit(method, Y, tr_, dn_, pre, post, opts$lambda, opts$nFactors, 0.5,
                    opts$bayesDraws, opts$mcLambda, NULL, 1, TRUE, NULL)
        fit <- fitf(tr, dn)
        summarise(fit$path, fit$att)
        out$nw <- if (is.null(fit$w)) NA_integer_ else sum(abs(fit$w) > 1e-6)
        if (key == "bscm") {
            out$lower <- fit$cri[1]; out$upper <- fit$cri[2]
            out$basis <- "posterior 95% credible interval"
        } else if (inference) {
            ratios <- vapply(dn, function(j) {
                ck()
                f <- tryCatch(fitf(j, setdiff(dn, j)), error = function(e) NULL)
                if (is.null(f)) NA_real_ else f$ratio
            }, 0)
            ratios <- ratios[is.finite(ratios)]
            out$p <- mean(c(fit$ratio, ratios) >= fit$ratio)
            out$basis <- "permutation, post/pre RMSPE ratio (all donors)"
        }
    } else if (key == "ascm") {
        t0n <- length(pre); t_final <- t0n + length(post)
        fit <- asyn_fit(Y = Y, times = ctx$times, treated = tr, t_int = ctx$t0,
                        Z = NULL, progfunc = "ridge", scm = TRUE, fixedeff = FALSE, min_1se = TRUE)
        att <- asyn_predict(fit, att = TRUE)
        path <- as.vector(Y[, tr]) - att[seq_len(t_final)]
        summarise(path, mean(att[(t0n + 1):t_final]))
        out$nw <- sum(abs(fit$weights[, 1]) > 1e-8)
        if (inference) {
            set.seed(as.integer(opts$seed))
            res <- asyn_conformal_inf(fit, alpha = 0.05, type = "iid", ns = 1000L, grid_size = 50L, checkpoint = ck)
            out$p <- res$p_val[t_final + 1]
            out$basis <- "conformal p-value, average effect = 0 (no interval for the average)"
        }
    } else if (key == "sdid") {
        N0 <- length(dn); T0 <- length(pre)
        Ys <- t(Y[, c(dn, tr), drop = FALSE])
        est <- sdid_estimate(Ys, N0, T0)
        w <- attr(est, "weights")
        sc_t <- as.vector(t(c(-w$omega, 1)) %*% Ys)
        gap <- sc_t - sum(sc_t[seq_len(T0)] * w$lambda)
        summarise(as.vector(Y[, tr]) - gap, as.numeric(est))
        out$nw <- sum(as.vector(w$omega) > 1e-8)
        if (inference) {
            set.seed(as.integer(opts$seed))
            se <- sdid_se(est, "placebo", 200L, checkpoint = ck)
            if (is.finite(se)) {
                z <- stats::qnorm(0.975)
                out$se <- se; out$lower <- out$att - z * se; out$upper <- out$att + z * se
                out$p <- 2 * stats::pnorm(-abs(out$att / se))
                out$basis <- "placebo standard error, normal approximation"
            }
        }
    }
    out
}

compare_run <- function(self, ck = function() NULL) {
    opts <- self$options
    notes <- self$results$notes
    if (isTRUE(opts$showGuide)) self$results$guide$setContent(fam_guide_html("compare"))
    if (isTRUE(opts$showHelp)) self$results$helpText$setContent(fam_notes_html("compare"))
    if (is.null(opts$outcome) || is.null(opts$unit) || is.null(opts$time) ||
        is.null(opts$treatedUnit) || !nzchar(opts$treatedUnit)) {
        notes$setContent("<p>Select an outcome, a unit identifier, a time variable and a treated unit to begin.</p>")
        return(invisible(NULL))
    }
    excl <- trimws(strsplit(opts$donorExclude, ",", fixed = TRUE)[[1]])
    excl <- excl[nzchar(excl)]
    pan <- tryCatch(
        fam_panel(self$data, opts$outcome, opts$unit, opts$time, opts$treatedUnit,
                  opts$treatmentTime, excl, character()),
        error = function(e) jmvcore::reject(conditionMessage(e)))
    ctx <- list(Y = pan$Y, tr = pan$treated, dn = pan$donors, pre = pan$pre, post = pan$post,
                times = pan$times, t0 = opts$treatmentTime)
    keys <- names(compare_methods)[vapply(compare_methods, function(m) isTRUE(opts[[m$opt]]), NA)]
    if (length(keys) == 0) {
        notes$setContent("<p>Tick at least one method under <i>Methods to compare</i>.</p>")
        return(invisible(NULL))
    }
    o <- list(lambda = opts$lambda, nFactors = opts$nFactors, mcLambda = opts$mcLambda,
              bayesDraws = opts$bayesDraws, optStarts = opts$optStarts, seed = opts$seed)
    inference <- isTRUE(opts$inference)

    res <- list(); failed <- character()
    for (k in keys) {
        ck()
        r <- tryCatch(compare_one(k, ctx, o, inference, ck), error = function(e) {
            if (inherits(e, "jmvcoreCancel") || grepl("interrupt|cancel", conditionMessage(e), ignore.case = TRUE)) stop(e)
            failed <<- c(failed, paste0(compare_methods[[k]]$label, ": ", conditionMessage(e)))
            NULL
        })
        if (!is.null(r)) res[[k]] <- r
    }
    if (length(res) == 0) jmvcore::reject(paste("No method could be estimated.", paste(failed, collapse = " ")))

    notes$setContent(paste0(
        "<h3>Comparison of synthetic control estimators</h3><p>", length(res), " method(s) applied to the same panel, using ",
        "pre-treatment outcomes only (no covariates). Donor pool: ", length(ctx$dn), " units. Pre-treatment periods: ",
        length(ctx$pre), "; post-treatment periods: ", length(ctx$post), ". Each method uses the default settings of its own analysis ",
        "except where set under <i>Method settings</i>. Report all results, not only the most favourable one.</p>",
        if (length(failed) > 0) paste0("<p><b>Not estimated:</b> ", paste(failed, collapse = "; "), "</p>") else ""))

    tbl <- self$results$summary
    for (k in names(res)) {
        r <- res[[k]]
        tbl$addRow(rowKey = k, values = list(method = compare_methods[[k]]$label, att = r$att,
            rmspePre = r$rmspe_pre, rmspePost = r$rmspe_post, ratio = r$ratio,
            nw = if (is.na(r$nw)) NA_integer_ else as.integer(r$nw),
            se = r$se, lower = r$lower, upper = r$upper, p = r$p, basis = r$basis))
    }
    if (isTRUE(opts$showPath)) {
        pt <- self$results$pathTable
        for (i in seq_along(ctx$times)) {
            vals <- list(time = as.character(ctx$times[i]), treated = ctx$Y[i, ctx$tr])
            for (k in names(res)) vals[[k]] <- res[[k]]$path[i]
            pt$addRow(rowKey = as.character(ctx$times[i]), values = vals)
        }
    }
    cols <- vapply(compare_methods[names(res)], function(m) m$col, "")
    labs <- vapply(compare_methods[names(res)], function(m) m$label, "")
    self$results$pathsPlot$setState(list(times = ctx$times, treated = as.vector(ctx$Y[, ctx$tr]),
        paths = lapply(res, function(r) r$path), cols = cols, labels = labs,
        t0 = ctx$times[ctx$post[1]], label = opts$outcome))
    self$results$gapsPlot$setState(list(times = ctx$times, paths = lapply(res, function(r) r$gap),
        cols = cols, labels = labs, t0 = ctx$times[ctx$post[1]], label = opts$outcome))
    self$results$attPlot$setState(list(labels = labs, cols = cols,
        att = vapply(res, function(r) r$att, 0), lower = vapply(res, function(r) r$lower, 0),
        upper = vapply(res, function(r) r$upper, 0)))
    invisible(NULL)
}

# ---- figures ---------------------------------------------------------------------------
compare_legend <- function(labels, cols, extra = NULL) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    lt <- rep(c(2, 3, 4, 5, 6), length.out = length(labels))
    leg <- c(extra$label, labels)
    graphics::legend("center", legend = leg, col = c(extra$col, cols), lty = c(extra$lty, lt),
                     lwd = 2, ncol = min(3, length(leg)), bty = "n", cex = 0.9)
}

compare_draw_paths <- function(st) {
    if (is.null(st)) return(FALSE)
    graphics::layout(matrix(1:2, 2), heights = c(5, 1.2))
    on.exit(graphics::layout(1))
    graphics::par(mar = c(4.1, 4.1, 3.1, 1))
    yl <- range(c(st$treated, unlist(st$paths)), na.rm = TRUE)
    graphics::plot(st$times, st$treated, type = "l", lwd = 3, ylim = yl, xlab = "Time", ylab = st$label,
                   main = "Treated unit vs counterfactuals")
    lt <- rep(c(2, 3, 4, 5, 6), length.out = length(st$paths))
    for (i in seq_along(st$paths)) graphics::lines(st$times, st$paths[[i]], col = st$cols[i], lty = lt[i], lwd = 2)
    graphics::abline(v = st$t0, lty = 3)
    compare_legend(st$labels, st$cols, list(label = "Treated", col = "black", lty = 1))
    TRUE
}

compare_draw_gaps <- function(st) {
    if (is.null(st)) return(FALSE)
    graphics::layout(matrix(1:2, 2), heights = c(5, 1.2))
    on.exit(graphics::layout(1))
    graphics::par(mar = c(4.1, 4.1, 3.1, 1))
    rng <- range(unlist(st$paths), na.rm = TRUE)
    graphics::plot(st$times, st$paths[[1]], type = "n", ylim = rng, xlab = "Time",
                   ylab = paste("Gap in", st$label), main = "Gap: treated minus counterfactual")
    lt <- rep(c(2, 3, 4, 5, 6), length.out = length(st$paths))
    for (i in seq_along(st$paths)) graphics::lines(st$times, st$paths[[i]], col = st$cols[i], lty = lt[i], lwd = 2)
    graphics::abline(h = 0, lty = 3)
    graphics::abline(v = st$t0, lty = 3)
    compare_legend(st$labels, st$cols)
    TRUE
}

compare_draw_att <- function(st) {
    if (is.null(st)) return(FALSE)
    n <- length(st$att)
    lo <- ifelse(is.finite(st$lower), st$lower, st$att)
    hi <- ifelse(is.finite(st$upper), st$upper, st$att)
    graphics::par(mar = c(9, 12, 3.1, 1), xpd = FALSE)
    graphics::plot(st$att, rev(seq_len(n)), type = "n", xlim = range(c(lo, hi, st$att, 0)), ylim = c(0.5, n + 0.5),
                   yaxt = "n", ylab = "", xlab = "Average post-treatment effect (ATT)", main = "ATT by method")
    graphics::axis(2, at = rev(seq_len(n)), labels = st$labels, las = 1)
    graphics::abline(v = 0, lty = 3)
    for (i in seq_len(n)) {
        y <- n + 1 - i
        if (is.finite(st$lower[i]) && is.finite(st$upper[i]))
            graphics::segments(lo[i], y, hi[i], y, lwd = 2, col = st$cols[i])
        graphics::points(st$att[i], y, pch = 19, col = st$cols[i], cex = 1.3)
    }
    graphics::par(xpd = NA)
    graphics::legend("bottom", inset = c(0, -0.37), legend = c(st$labels, "95% interval (if available)"),
                     col = c(st$cols, "grey30"), pch = c(rep(19, n), NA), lty = c(rep(NA, n), 1), lwd = 2,
                     ncol = 3, bty = "n", cex = 0.8)
    TRUE
}
