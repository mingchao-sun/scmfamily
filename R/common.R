# Shared engine for the synthetic-control family (module "scmfamily").
# Estimation is pure R (needs quadprog only), so it can be tested without jamovi.
# Each analysis calls fam_run(self, method) with one of:
#   "penalized", "demeaned", "ife", "enet", "bscm", "mc".
# (The "sdid", "ascm" and "stagascm" analyses have their own engines: sdid_engine.R,
# augsynth_engine.R and multisynth_engine.R.)

fam_info <- list(
    penalized = list(
        title = "Penalized synthetic control",
        cite  = "Abadie and L'Hour (2021), Journal of the American Statistical Association",
        note  = "Donor weights minimise the squared pre-treatment predictor gap plus a penalty on the distance between the treated unit and each donor. Predictors are the pre-treatment outcomes, standardised by their cross-unit SD; lambda is in those units."),
    demeaned = list(
        title = "Demeaned synthetic control",
        cite  = "Ferman and Pinto (2021), Quantitative Economics",
        note  = "Each unit's pre-treatment mean is removed before the simplex weights are fitted, so donors can match trends when levels differ. The counterfactual adds back the treated unit's pre-treatment mean."),
    ife = list(
        title = "Interactive fixed effects (generalized synthetic control)",
        cite  = "Xu (2017), Political Analysis",
        note  = "Two-way fixed effects and r latent factors are estimated from the donor pool only; the treated unit's factor loadings come from a pre-treatment regression. The number of factors is user-set; compare several values as a robustness check. This is a control-only estimator, not the full gsynth routine with cross-validated r."),
    enet = list(
        title = "Elastic-net synthetic control",
        cite  = "Doudchenko and Imbens (2016), NBER Working Paper 22791",
        note  = "Weights are unrestricted (they may be negative and need not sum to one) and an intercept is included. The penalty mixes ridge and lasso through alpha. As in the paper, lambda (and alpha, if chosen by cross-validation) minimises the squared prediction error in the post-treatment periods when each donor in turn plays the treated unit."),
    bscm = list(
        title = "Bayesian synthetic control (horseshoe shrinkage)",
        cite  = "Kim, Lee and Gupta (2020), Journal of Marketing Research 57(5)",
        note  = "Donor coefficients carry horseshoe shrinkage priors and are not constrained to the simplex; the intercept is free. The ATT is the posterior mean with a 95% credible interval that includes the sampling error of the post-treatment outcomes. The weights table shows posterior-mean coefficients, which need not be non-negative or sum to one. This regression form follows Kim et al. (2020) in spirit; check the exact prior specification against the paper before publication."),
    mc = list(
        title = "Matrix completion (nuclear-norm regularised)",
        cite  = "Athey, Bayati, Doudchenko, Imbens and Khosravi (2021), Journal of the American Statistical Association 116(536)",
        note  = "The treated unit's post-treatment outcomes are treated as missing and imputed from a low-rank matrix plus unit and time fixed effects. Unit and time effects are estimated jointly with the low-rank matrix by coordinate descent, as in the authors' MCPanel package. The nuclear-norm penalty is relative to the smallest penalty that sets the low-rank part to zero; check a few values. Weights are not reported because the estimator is not a weighted average of donors.")
)

# min 0.5 w'Qw + cvec'w  subject to  sum(w) = 1, w >= 0
fam_qp_simplex <- function(Q, cvec) {
    if (!requireNamespace("quadprog", quietly = TRUE))
        stop("The R package 'quadprog' is required. Install it with install.packages('quadprog').")
    cvec <- as.vector(cvec)
    n <- length(cvec)
    Q <- (Q + t(Q)) / 2
    Q <- Q + diag(max(1e-10, 1e-10 * mean(abs(diag(Q)))), n)
    sol <- quadprog::solve.QP(Dmat = Q, dvec = -cvec,
                              Amat = cbind(rep(1, n), diag(n)),
                              bvec = c(1, rep(0, n)), meq = 1)
    w <- pmax(sol$solution, 0)
    w / sum(w)
}

# Build the balanced T x N outcome matrix and the pre/post indices. Optional covariates
# are returned as a T x N x K array (NA allowed; each estimator decides how to use them).
fam_panel <- function(df, y, unit, time, treated, t0, exclude = character(),
                      covars = character()) {
    covars <- as.character(covars)
    tt <- df[[time]]
    if (is.factor(tt)) tt <- as.character(tt)
    tt <- suppressWarnings(as.numeric(as.character(tt)))
    if (anyNA(tt)) stop("The time variable must be numeric (for example a year).")
    d <- data.frame(u = as.character(df[[unit]]), t = tt,
                    y = as.numeric(df[[y]]), stringsAsFactors = FALSE)
    cm <- matrix(NA_real_, nrow(df), length(covars), dimnames = list(NULL, covars))
    for (k in covars) cm[, k] <- suppressWarnings(as.numeric(as.character(df[[k]])))
    cc <- stats::complete.cases(d)
    d <- d[cc, , drop = FALSE]
    cm <- cm[cc, , drop = FALSE]
    if (!(treated %in% d$u)) stop("The treated unit was not found in the data.")
    unknown <- setdiff(exclude, unique(d$u))
    if (length(unknown) > 0)
        stop(paste0("Unknown unit(s) in 'Exclude donors': ", paste(unknown, collapse = ", ")))
    keep <- !(d$u %in% setdiff(exclude, treated))
    d <- d[keep, , drop = FALSE]
    cm <- cm[keep, , drop = FALSE]
    if (anyDuplicated(d[, c("u", "t")]))
        stop("Each unit must have one row per period (duplicate unit-time rows found).")
    times <- sort(unique(d$t))
    if (!(t0 %in% times)) stop("The first treated period must be one of the observed time periods.")
    units <- sort(unique(d$u))
    Y <- matrix(NA_real_, nrow = length(times), ncol = length(units),
                dimnames = list(NULL, units))
    ix <- cbind(match(d$t, times), match(d$u, units))
    Y[ix] <- d$y
    if (anyNA(Y)) stop("The panel is unbalanced: every unit needs a value in every period.")
    cov <- NULL
    if (length(covars) > 0) {
        cov <- array(NA_real_, dim = c(length(times), length(units), length(covars)),
                     dimnames = list(NULL, units, covars))
        for (k in seq_along(covars))
            cov[cbind(ix, k)] <- cm[, k]
    }
    pre <- which(times < t0)
    post <- which(times >= t0)
    if (length(pre) < 2) stop("At least two pre-treatment periods are required.")
    if (length(post) < 1) stop("No post-treatment periods: check the first treated period.")
    donors <- setdiff(units, treated)
    if (length(donors) < 2) stop("At least two donor units are required.")
    list(Y = Y, times = times, pre = pre, post = post, treated = treated, donors = donors,
         cov = cov)
}

# Covariate predictor rows (K x units): pre-treatment covariate means, z-scored across the
# units in `cols`, then scaled so that the covariate block carries `covw` times the weight
# of the outcome block. `scale` is the SD scale of one outcome row (1 when outcomes are
# already standardised).
fam_cov_rows <- function(cov, pre, cols, covw = 1, scale = 1) {
    K <- dim(cov)[3]
    nm <- dimnames(cov)[[3]]
    M <- matrix(apply(cov[pre, cols, , drop = FALSE], c(2, 3), mean, na.rm = TRUE),
                length(cols), K, dimnames = list(cols, nm))
    bad <- !is.finite(M)
    if (any(bad)) {
        k <- which(colSums(bad) > 0)[1]
        stop(sprintf("Covariate '%s' has no pre-treatment values for unit '%s'.",
                     nm[k], cols[which(bad[, k])[1]]), call. = FALSE)
    }
    sdv <- apply(M, 2, stats::sd)
    sdv[!is.finite(sdv) | sdv <= 0] <- 1
    Z <- sweep(sweep(M, 2, colMeans(M)), 2, sdv, "/")
    t(Z) * (scale * sqrt(max(covw, 0) * length(pre) / K))
}

# Penalized SC: min ||x1 - X0 w||^2 + lambda * sum_j w_j ||x1 - x_j||^2 on the simplex.
fam_pen <- function(Y, tr, dn, pre, lambda, cov = NULL, covw = 1) {
    X <- Y[pre, , drop = FALSE]
    s <- apply(X[, dn, drop = FALSE], 1, stats::sd)
    s[!is.finite(s) | s <= 0] <- 1
    Xs <- X / s
    if (!is.null(cov)) {
        cols <- c(tr, dn)
        Xs <- rbind(Xs[, cols, drop = FALSE], fam_cov_rows(cov, pre, cols, covw, 1))
    }
    x1 <- Xs[, tr]
    X0 <- Xs[, dn, drop = FALSE]
    d <- colSums((X0 - x1)^2)
    w <- fam_qp_simplex(2 * crossprod(X0), -2 * crossprod(X0, x1) + lambda * d)
    list(path = as.vector(Y[, dn, drop = FALSE] %*% w), w = w)
}

# Demeaned SC: weights fitted to pre-period deviations from each unit's pre-mean.
fam_dsc <- function(Y, tr, dn, pre, cov = NULL, covw = 1) {
    mu <- colMeans(Y[pre, , drop = FALSE])
    Yc <- sweep(Y, 2, mu)
    B <- Yc[pre, dn, drop = FALSE]
    a <- Yc[pre, tr]
    if (!is.null(cov)) {
        sc <- mean(apply(Yc[pre, , drop = FALSE], 1, stats::sd))
        Cr <- fam_cov_rows(cov, pre, c(tr, dn), covw, sc)
        B <- rbind(B, Cr[, dn, drop = FALSE])
        a <- c(a, Cr[, tr])
    }
    w <- fam_qp_simplex(2 * crossprod(B), -2 * crossprod(B, a))
    list(path = as.vector(Yc[, dn, drop = FALSE] %*% w) + mu[[tr]], w = w)
}

# Interactive fixed effects / generalized SC (Xu 2017), control-only version.
# Two-way fixed effects and r factors come from the donor pool; the treated unit's
# loadings come from OLS on the pre-treatment periods.
fam_ife <- function(Y, tr, dn, pre, post, r = 1, cov = NULL) {
    C <- Y[, dn, drop = FALSE]
    Tn <- nrow(C)
    n0 <- ncol(C)
    r <- as.integer(r)
    if (is.na(r) || r < 1 || r > length(pre) - 2)
        stop("The number of factors must be at least 1 and at most (pre-treatment periods - 2).")
    twoway <- function(M) M - rowMeans(M) - matrix(colMeans(M), nrow(M), ncol(M), byrow = TRUE) + mean(M)
    beta <- NULL
    if (is.null(cov)) {
        mC <- mean(C)
        xi <- rowMeans(C) - mC
        R <- twoway(C)
        Fm <- svd(R)$u[, seq_len(r), drop = FALSE]
        z <- Y[pre, tr] - xi[pre]
        D <- cbind(1, Fm[pre, , drop = FALSE])
        coef <- stats::lm.fit(D, z)$coefficients
        path <- as.vector(xi + coef[1] + Fm %*% coef[-1])
    } else {
        # Xu (2017): common slopes on time-varying covariates, estimated jointly with the
        # factors from the donor pool (alternate between OLS and a rank-r SVD).
        K <- dim(cov)[3]
        nm <- dimnames(cov)[[3]]
        Xc <- cov[, dn, , drop = FALSE]
        Xt <- matrix(cov[, tr, ], Tn, K)
        if (anyNA(Xc) || anyNA(Xt)) {
            miss <- nm[vapply(seq_len(K), function(k) anyNA(Xc[, , k]) || anyNA(Xt[, k]), NA)]
            stop("Interactive fixed effects needs covariates observed in every period for every unit. ",
                 "Missing values in: ", paste(miss, collapse = ", "), ".", call. = FALSE)
        }
        Xmat <- vapply(seq_len(K), function(k) as.vector(twoway(Xc[, , k])), numeric(Tn * n0))
        Xmat <- matrix(Xmat, Tn * n0, K)
        XtX <- crossprod(Xmat)
        if (qr(XtX)$rank < K)
            stop("The covariates are collinear after removing unit and time effects.", call. = FALSE)
        Cw <- twoway(C)
        beta <- solve(XtX, crossprod(Xmat, as.vector(Cw)))
        for (it in 1:500) {
            fit_x <- Reduce(`+`, lapply(seq_len(K), function(k) Xc[, , k] * beta[k]))
            sv <- svd(twoway(C - fit_x))
            L <- sv$u[, seq_len(r), drop = FALSE] %*%
                (sv$d[seq_len(r)] * t(sv$v[, seq_len(r), drop = FALSE]))
            bn <- solve(XtX, crossprod(Xmat, as.vector(Cw - L)))
            done <- max(abs(bn - beta)) < 1e-8 * max(1, max(abs(beta)))
            beta <- bn
            if (done) break
        }
        beta <- as.vector(beta)
        names(beta) <- nm
        fit_x <- Reduce(`+`, lapply(seq_len(K), function(k) Xc[, , k] * beta[k]))
        Cres <- C - fit_x
        xi <- rowMeans(Cres) - mean(Cres)
        sv <- svd(twoway(Cres))
        Fm <- sv$u[, seq_len(r), drop = FALSE]
        xt <- as.vector(Xt %*% beta)
        z <- Y[pre, tr] - xt[pre] - xi[pre]
        D <- cbind(1, Fm[pre, , drop = FALSE])
        coef <- stats::lm.fit(D, z)$coefficients
        path <- as.vector(xt + xi + coef[1] + Fm %*% coef[-1])
    }
    list(path = path, w = NULL, beta = beta)
}

# Elastic-net synthetic control (Doudchenko and Imbens 2016, section 4.4). The intercept mu and the
# weights omega minimise  sum_t (Y_1t - mu - omega' Y_0t)^2 + lambda * ((1-alpha)/2 ||omega||_2^2 + alpha ||omega||_1)
# over the pre-treatment periods, with an unpenalised intercept, no standardisation of the donors
# and no constraints on the weights. The penalty is chosen by the paper's cross-validation: each
# donor in turn is the pseudo-treated unit, its pre-treatment outcomes are fitted on the other
# donors, and the squared prediction error in the post-treatment periods is averaged over donors.
# alpha is either fixed by the user or searched on {0.1, ..., 0.9}.
#
# glmnet is used as the solver. It standardises the response internally, which changes the ratio
# of the ridge to the lasso part of the penalty, so the outcome and the donors are first divided by
# the response's SD s (a common scale leaves the weights unchanged); glmnet's lambda is then
# lambda_paper / (2 * T0 * s^2) and its objective is exactly the paper's.
fam_enet_fit <- function(X, y, alpha, lambda_paper) {
    n <- length(y)
    sc <- sqrt(mean((y - mean(y))^2))
    if (!is.finite(sc) || sc <= 0) sc <- 1
    fit <- suppressWarnings(glmnet::glmnet(X / sc, y / sc, alpha = alpha, standardize = FALSE,
                                           lambda = lambda_paper / (2 * n * sc^2),
                                           thresh = 1e-12, maxit = 1e6))
    list(fit = fit, scale = sc, n = n)
}
fam_enet_coef <- function(ef, lambda_paper) {
    cf <- as.matrix(stats::coef(ef$fit, s = lambda_paper / (2 * ef$n * ef$scale^2)))
    list(w = cf[-1, 1], mu = cf[1, 1] * ef$scale)
}
fam_enet <- function(Y, tr, dn, pre, post, alpha = 0.5, cv_alpha = FALSE, fixed = NULL) {
    if (!requireNamespace("glmnet", quietly = TRUE))
        stop("The R package 'glmnet' is required. Install it with install.packages('glmnet').")
    if (!cv_alpha && (!is.finite(alpha) || alpha < 0 || alpha > 1))
        stop("The mixing parameter alpha must lie in [0, 1].")
    if (!is.null(fixed)) {                           # penalty given (conformal inference): no cross-validation
        ef <- fam_enet_fit(Y[pre, dn, drop = FALSE], Y[pre, tr], fixed$alpha, fixed$lambda)
        cf <- fam_enet_coef(ef, fixed$lambda)
        w <- cf$w; names(w) <- dn
        return(list(path = as.vector(cf$mu + Y[, dn, drop = FALSE] %*% w), w = w,
                    extra = list(alpha = fixed$alpha, lambda = fixed$lambda, cv_error = NA_real_, intercept = cf$mu)))
    }
    if (length(dn) < 3) stop("At least three donor units are required for the elastic-net cross-validation.")
    fit_pre <- pre
    hold <- setdiff(post, pre)                       # periods used to score the penalty
    if (length(hold) == 0) {                         # all periods are used for fitting (conformal inference)
        k <- max(1, ceiling(length(pre) / 5))
        hold <- utils::tail(pre, k)
        fit_pre <- utils::head(pre, length(pre) - k)
        if (length(fit_pre) < 3) stop("Too few periods for the elastic-net cross-validation.")
    }
    alphas <- if (cv_alpha) seq(0.1, 0.9, by = 0.1) else alpha
    best <- list(err = Inf)
    for (a in alphas) {
        # common grid of penalties, taken from the path of the treated unit's own fit
        nf <- length(fit_pre)
        scy <- sqrt(mean((Y[fit_pre, tr] - mean(Y[fit_pre, tr]))^2)); if (!is.finite(scy) || scy <= 0) scy <- 1
        path <- suppressWarnings(glmnet::glmnet(Y[fit_pre, dn, drop = FALSE] / scy, Y[fit_pre, tr] / scy, alpha = a,
                                                standardize = FALSE, thresh = 1e-12, maxit = 1e6))
        lams <- path$lambda * 2 * nf * scy^2         # paper scale
        errs <- matrix(NA_real_, length(dn), length(lams))
        for (jj in seq_along(dn)) {
            others <- dn[-jj]
            ef <- fam_enet_fit(Y[fit_pre, others, drop = FALSE], Y[fit_pre, dn[jj]], a, lams)
            pr <- stats::predict(ef$fit, newx = Y[hold, others, drop = FALSE] / ef$scale) * ef$scale
            errs[jj, seq_len(ncol(pr))] <- colMeans((Y[hold, dn[jj]] - pr)^2)
        }
        cvm <- colMeans(errs, na.rm = TRUE)
        k <- which.min(cvm)
        if (cvm[k] < best$err) best <- list(err = cvm[k], alpha = a, lambda = lams[k], lams = lams)
    }
    ef <- fam_enet_fit(Y[pre, dn, drop = FALSE], Y[pre, tr], best$alpha, best$lams)
    cf <- fam_enet_coef(ef, best$lambda)
    w <- cf$w
    names(w) <- dn
    list(path = as.vector(cf$mu + Y[, dn, drop = FALSE] %*% w), w = w,
         extra = list(alpha = best$alpha, lambda = best$lambda, cv_error = best$err, intercept = cf$mu))
}

# Bayesian synthetic control (Kim, Lee and Gupta 2020): donor coefficients with
# horseshoe shrinkage priors, no simplex constraint, free intercept. Gibbs sampler with
# the auxiliary-variable scheme of Makalic and Schmidt (2016).
fam_bscm <- function(Y, tr, dn, pre, post, ndraw = 1000, seed = 1) {
    burn <- ndraw %/% 2
    Z <- Y[pre, dn, drop = FALSE]
    y <- Y[pre, tr]
    zb <- colMeans(Z)
    yb <- mean(y)
    Zc <- sweep(Z, 2, zb)
    yc <- y - yb
    n <- length(yc)
    p <- ncol(Zc)
    # the sampler works on standardised data (unit-SD donors and outcome), as is usual for horseshoe
    # priors; coefficients are transformed back to the donors' own scale
    zs <- apply(Zc, 2, stats::sd); zs[!is.finite(zs) | zs <= 0] <- 1
    ys <- stats::sd(yc); if (!is.finite(ys) || ys <= 0) ys <- 1
    Zs <- sweep(Zc, 2, zs, "/")
    ysd <- yc / ys
    set.seed(seed)
    beta <- rep(0, p)
    sig2 <- stats::var(ysd)
    lam2 <- rep(1, p)
    nu <- rep(1, p)
    tau2 <- 1
    xi <- 1
    nk <- ndraw + burn
    Bdr <- matrix(NA_real_, nk, p)
    sdr <- numeric(nk)
    for (it in seq_len(nk)) {
        # beta | rest ~ N((X'X + D^-1)^-1 X'y, sig2 (X'X + D^-1)^-1), D = diag(lam2 * tau2): drawn with the
        # n x n system of Bhattacharya, Chakraborty and Mallick (2016), which stays stable when p > n
        d <- pmin(pmax(lam2 * tau2, 1e-12), 1e12)
        u <- sqrt(d) * stats::rnorm(p)
        v <- as.vector(Zs %*% u) + stats::rnorm(n)
        wv <- solve(Zs %*% (d * t(Zs)) + diag(n), ysd / sqrt(sig2) - v)
        beta <- sqrt(sig2) * (u + d * as.vector(crossprod(Zs, wv)))
        rss <- sum((ysd - Zs %*% beta)^2)
        sig2 <- 1 / stats::rgamma(1, shape = (n + p) / 2, rate = (rss + sum(beta^2 / d)) / 2)
        lam2 <- 1 / stats::rgamma(p, shape = 1, rate = 1 / nu + beta^2 / (2 * sig2 * tau2))
        nu <- 1 / stats::rgamma(p, shape = 1, rate = 1 + 1 / lam2)
        tau2 <- 1 / stats::rgamma(1, shape = (p + 1) / 2, rate = 1 / xi + sum(beta^2 / lam2) / (2 * sig2))
        xi <- 1 / stats::rgamma(1, shape = 1, rate = 1 + 1 / tau2)
        Bdr[it, ] <- beta * ys / zs
        sdr[it] <- sqrt(sig2) * ys
    }
    Bk <- Bdr[(burn + 1):nk, , drop = FALSE]
    bm <- colMeans(Bk)
    a <- yb - sum(zb * bm)
    Yd <- Y[, dn, drop = FALSE]
    a_k <- yb - as.vector(Bk %*% zb)
    cf_post_k <- a_k + as.vector(Bk %*% colMeans(Yd[post, , drop = FALSE]))
    # The effect is the constant shift of the treated outcome after the treatment: besides the uncertainty of the
    # coefficients, the posterior includes the sampling error of the average of the post-treatment errors, N(0, sigma^2 / T1).
    sd_k <- sdr[(burn + 1):nk]
    att_k <- mean(Y[post, tr]) - cf_post_k + stats::rnorm(length(sd_k), 0, sd_k / sqrt(length(post)))
    list(path = as.vector(a + Yd %*% bm), w = bm, tau = mean(att_k),
         cri = stats::quantile(att_k, c(0.025, 0.975), names = FALSE))
}

# Dispatch and summarise one fit. Gaps are treated minus counterfactual path.
fam_fit <- function(method, Y, tr, dn, pre, post, lambda = 0.01,
                    nfac = 1, alpha = 0.5, ndraw = 1000, mcl = 0.1,
                    cov = NULL, covw = 1, acv = FALSE, enfix = NULL) {
    out <- switch(method,
        penalized = fam_pen(Y, tr, dn, pre, lambda, cov, covw),
        demeaned  = fam_dsc(Y, tr, dn, pre, cov, covw),
        ife       = fam_ife(Y, tr, dn, pre, post, nfac, cov),
        enet      = fam_enet(Y, tr, dn, pre, post, alpha, acv, enfix),
        bscm      = fam_bscm(Y, tr, dn, pre, post, ndraw),
        mc        = fam_mc(Y, tr, dn, pre, post, mcl),
        stop("Unknown method: ", method))
    gap <- as.vector(Y[, tr]) - out$path
    if (!is.null(out$w)) names(out$w) <- dn
    pre_gap <- gap[pre]
    post_gap <- gap[post]
    rmspe_pre <- sqrt(mean(pre_gap^2))
    rmspe_post <- sqrt(mean(post_gap^2))
    att <- if (is.null(out$tau)) mean(post_gap) else out$tau
    list(method = method, path = out$path, gap = gap, w = out$w, att = att,
         rmspe_pre = rmspe_pre, rmspe_post = rmspe_post,
         ratio = rmspe_post / rmspe_pre, cri = out$cri, beta = out$beta, extra = out$extra)
}

# jmvcore's options$name errors when an option is not defined for this analysis
# (e.g. lambda in the ASCM analysis), so read it defensively.
fam_opt <- function(opts, name, default) {
    v <- tryCatch(opts[[name]], error = function(e) NULL)
    if (is.null(v)) default else v
}

# Shared analysis body used by all four jamovi analyses.
fam_run <- function(self, method, ck = function() NULL) {
    opts <- self$options
    info <- fam_info[[method]]
    notes <- self$results$notes
    if (isTRUE(opts$showGuide)) self$results$guide$setContent(fam_guide_html(method))
    if (isTRUE(opts$showHelp)) self$results$helpText$setContent(fam_notes_html(method))
    if (is.null(opts$outcome) || is.null(opts$unit) || is.null(opts$time) ||
        is.null(opts$treatedUnit) || !nzchar(opts$treatedUnit)) {
        notes$setContent("<p>Select an outcome, a unit identifier, a time variable and a treated unit to begin.</p>")
        return(invisible(NULL))
    }
    excl <- trimws(strsplit(opts$donorExclude, ",", fixed = TRUE)[[1]])
    excl <- excl[nzchar(excl)]
    covars <- as.character(fam_opt(opts, "predictors", character()))
    pan <- tryCatch(
        fam_panel(self$data, opts$outcome, opts$unit, opts$time, opts$treatedUnit,
                  opts$treatmentTime, excl, covars),
        error = function(e) jmvcore::reject(conditionMessage(e)))
    lam <- fam_opt(opts, "lambda", 0.01)
    nfc <- fam_opt(opts, "nFactors", 1)
    enA <- fam_opt(opts, "enetAlpha", 0.5)
    acv <- isTRUE(fam_opt(opts, "enetCV", FALSE))
    ndr <- fam_opt(opts, "bayesDraws", 1000)
    mcl <- fam_opt(opts, "mcLambda", 0.1)
    covw <- fam_opt(opts, "covWeight", 1)
    cov <- pan$cov
    Y <- pan$Y; tr <- pan$treated; dn <- pan$donors
    pre <- pan$pre; post <- pan$post; times <- pan$times
    t0 <- times[post[1]]

    fit <- tryCatch(fam_fit(method, Y, tr, dn, pre, post, lam, nfc, enA, ndr, mcl, cov, covw, acv),
                    error = function(e) jmvcore::reject(conditionMessage(e)))

    notes_txt <- paste0(
        "<h3>", info$title, "</h3><p>", info$note, "</p><p>Reference: ", info$cite,
        ". Donor pool: ", length(dn), " units. Pre-treatment periods: ", length(pre),
        "; post-treatment periods: ", length(post), ".</p>")
    if (!is.null(cov)) {
        notes_txt <- paste0(notes_txt, "<p>Covariates: ", paste(covars, collapse = ", "), ". ",
            if (method == "ife")
                "They enter as time-varying regressors with coefficients common to all units, estimated jointly with the factors. Covariates must not be affected by the treatment: post-treatment values that respond to it absorb part of the effect."
            else sprintf("Their pre-treatment means are standardised across units and matched together with the pre-treatment outcomes (covariate block weight %s relative to the outcome block).", format(covw)),
            "</p>")
    }
    notes$setContent(notes_txt)

    tbl <- self$results$summary
    tbl$addRow(rowKey = "att", values = list(metric = "Average treatment effect on the treated (ATT)", value = fit$att))
    tbl$addRow(rowKey = "pre", values = list(metric = "RMSPE, pre-treatment", value = fit$rmspe_pre))
    tbl$addRow(rowKey = "post", values = list(metric = "RMSPE, post-treatment", value = fit$rmspe_post))
    tbl$addRow(rowKey = "ratio", values = list(metric = "Post / pre RMSPE ratio", value = fit$ratio))
    if (!is.null(fit$cri)) {
        tbl$addRow(rowKey = "cri_lo", values = list(metric = "ATT 95% credible interval, lower", value = fit$cri[1]))
        tbl$addRow(rowKey = "cri_hi", values = list(metric = "ATT 95% credible interval, upper", value = fit$cri[2]))
    }
    if (method == "enet" && !is.null(fit$extra)) {
        tbl$addRow(rowKey = "en_alpha", values = list(metric = "Elastic-net mixing alpha (selected)", value = fit$extra$alpha))
        tbl$addRow(rowKey = "en_lambda", values = list(metric = "Elastic-net penalty lambda (selected, paper scale)", value = fit$extra$lambda))
        tbl$addRow(rowKey = "en_int", values = list(metric = "Intercept", value = fit$extra$intercept))
    }
    if (!is.null(fit$w)) {
        tbl$addRow(rowKey = "nw", values = list(metric = "Donors with nonzero weight", value = sum(abs(fit$w) > 1e-6)))
        w <- sort(fit$w, decreasing = TRUE)
        for (u in names(w))
            self$results$weights$addRow(rowKey = u, values = list(unit = u, weight = as.numeric(w[[u]])))
    }

    if (!is.null(cov) && method %in% c("penalized", "demeaned")) {
        M <- matrix(apply(cov[pre, c(tr, dn), , drop = FALSE], c(2, 3), mean, na.rm = TRUE),
                    length(dn) + 1, length(covars), dimnames = list(c(tr, dn), covars))
        bt <- self$results$balance
        for (k in covars)
            bt$addRow(rowKey = k, values = list(predictor = k, treated = M[tr, k],
                      synthetic = sum(fit$w[dn] * M[dn, k])))
    }
    if (!is.null(fit$beta)) {
        ct <- self$results$covCoef
        for (k in names(fit$beta))
            ct$addRow(rowKey = k, values = list(covariate = k, estimate = fit$beta[[k]]))
    }

    if (isTRUE(fam_opt(opts, "conformal", FALSE))) {
        if (method %in% c("penalized", "demeaned", "enet", "ife")) {
            pars <- list(lambda = lam, nfac = nfc, alpha = enA, ndraw = ndr, cov = cov, covw = covw, acv = acv,
                         enfix = if (method == "enet" && !is.null(fit$extra)) list(alpha = fit$extra$alpha, lambda = fit$extra$lambda) else NULL)
            cf <- tryCatch(fam_conformal(method, Y, tr, dn, pre, post, fit$att, fit$gap, pars, ck = ck),
                           error = function(e) NULL)
            if (!is.null(cf)) {
                tbl$addRow(rowKey = "conf_p", values = list(metric = "Conformal p-value, ATT = 0 (Chernozhukov et al., 2021)", value = cf$p0))
                tbl$addRow(rowKey = "conf_lo", values = list(metric = "Conformal 95% interval, lower", value = cf$lower))
                tbl$addRow(rowKey = "conf_hi", values = list(metric = "Conformal 95% interval, upper", value = cf$upper))
            }
        } else {
            notes$setContent(paste0(notes_txt, "<p>Conformal inference is available for penalized, demeaned, elastic-net and interactive fixed-effects estimators only.</p>"))
        }
    }
    if (isTRUE(opts$showPath)) {
        pt <- self$results$pathTable
        for (i in seq_along(times))
            pt$addRow(rowKey = as.character(times[i]),
                      values = list(time = as.character(times[i]), treated = Y[i, tr],
                                    synthetic = fit$path[i], gap = fit$gap[i]))
    }

    plac <- NULL
    if (isTRUE(opts$placebo)) {
        plac <- list()
        for (j in dn) {
            ck()
            f <- tryCatch(fam_fit(method, Y, j, setdiff(dn, j), pre, post, lam, nfc, enA, ndr, mcl, cov, covw, acv),
                          error = function(e) NULL)
            if (!is.null(f)) { f$unit <- j; plac[[j]] <- f }
        }
        if (length(plac) < 1) jmvcore::reject("No placebo fits could be computed.")
        all_ratio <- c(fit$ratio, vapply(plac, function(f) f$ratio, 0))
        all_att <- c(fit$att, vapply(plac, function(f) f$att, 0))
        all_pre <- c(fit$rmspe_pre, vapply(plac, function(f) f$rmspe_pre, 0))
        incl <- all_pre <= opts$placeboPreRatio * fit$rmspe_pre
        incl[1] <- TRUE
        p_ratio <- mean(all_ratio >= fit$ratio)
        p_filt <- mean(all_ratio[incl] >= fit$ratio)
        p_att <- mean(abs(all_att) >= abs(fit$att))

        inf <- self$results$inference
        inf$addRow(rowKey = "p_ratio", values = list(metric = "Permutation p, RMSPE ratio (all donors)", value = p_ratio))
        inf$addRow(rowKey = "p_filt", values = list(metric = sprintf("Permutation p, RMSPE ratio (pre-fit filter %.1f)", opts$placeboPreRatio), value = p_filt))
        inf$addRow(rowKey = "p_att", values = list(metric = "Permutation p, |ATT| (two-sided)", value = p_att))
        inf$addRow(rowKey = "n_kept", values = list(metric = "Placebo units kept by the pre-fit filter", value = sum(incl) - 1))

        pl <- self$results$placeboTable
        pl$addRow(rowKey = "__treated__", values = list(unit = paste0(tr, " (treated)"), ratio = fit$ratio,
                  preRmspe = fit$rmspe_pre, att = fit$att, included = "yes"))
        for (f in plac)
            pl$addRow(rowKey = f$unit, values = list(unit = f$unit, ratio = f$ratio, preRmspe = f$rmspe_pre,
                      att = f$att, included = if (f$rmspe_pre <= opts$placeboPreRatio * fit$rmspe_pre) "yes" else "no"))
    }

    self$results$pathsPlot$setState(list(times = times, treated = Y[, tr], path = fit$path,
                                         t0 = t0, label = opts$outcome))
    gaps <- if (is.null(plac)) list() else lapply(plac, function(f) f$gap)
    self$results$gapsPlot$setState(list(times = times, gap = fit$gap, placebo = gaps,
                                        t0 = t0, label = opts$outcome))
    invisible(NULL)
}

fam_draw_paths <- function(st) {
    if (is.null(st)) return(FALSE)
    yl <- range(c(st$treated, st$path), na.rm = TRUE)
    op <- graphics::par(mar = c(6.5, 4.1, 4.1, 2.1), xpd = FALSE)
    on.exit(graphics::par(op))
    graphics::plot(st$times, st$treated, type = "l", lwd = 2, ylim = yl,
                   xlab = "Time", ylab = st$label,
                   main = "Treated unit vs counterfactual")
    graphics::lines(st$times, st$path, lty = 2, lwd = 2, col = "steelblue")
    graphics::abline(v = st$t0, lty = 3)
    graphics::par(xpd = NA)
    graphics::legend("bottom", inset = c(0, -0.33), horiz = TRUE,
                     c("Treated", "Counterfactual"), lty = c(1, 2), lwd = 2,
                     col = c("black", "steelblue"), bty = "n")
    TRUE
}

fam_draw_gaps <- function(st) {
    if (is.null(st)) return(FALSE)
    rng <- range(c(st$gap, unlist(st$placebo)), na.rm = TRUE)
    graphics::plot(st$times, st$gap, type = "n", ylim = rng, xlab = "Time",
                   ylab = paste("Gap in", st$label), main = "Gap: treated minus counterfactual")
    for (g in st$placebo) graphics::lines(st$times, g, col = "grey70")
    graphics::lines(st$times, st$gap, lwd = 2)
    graphics::abline(h = 0, lty = 3)
    graphics::abline(v = st$t0, lty = 3)
    TRUE
}


# Matrix completion with nuclear-norm regularisation (Athey, Bayati, Doudchenko, Imbens and
# Khosravi 2021, MC-NNM). Minimises, over L, unit effects u and time effects v,
#   (1/|O|) sum_{(i,t) in O} (Y_it - L_it - u_i - v_t)^2 + lambda_L ||L||_*
# by cyclic coordinate descent (u, then v, then a singular-value-thresholding step for L),
# the algorithm of the authors' MCPanel package. O is the set of observed cells: the treated
# unit's post-treatment cells are missing. lam_rel = lambda_L / lambda_L,max, where lambda_L,max
# is the smallest penalty that makes L = 0, so lam_rel lies in (0, 1].
fam_mc <- function(Y, tr, dn, pre, post, lam_rel = 0.1, maxit = 20000, tol = 1e-12) {
    O <- Y[, c(tr, dn), drop = FALSE]
    O[post, 1] <- NA
    obs <- !is.na(O)
    M <- ifelse(obs, O, 0)
    nobs <- sum(obs)
    cnt_u <- colSums(obs)                              # observed cells per unit
    cnt_v <- rowSums(obs)                              # observed cells per period
    objective <- function(L, u, v, sing_sum, lam) {
        est <- L + matrix(u, nrow(M), ncol(M), byrow = TRUE) + matrix(v, nrow(M), ncol(M))
        sum(((est - M) * obs)^2) / nobs + lam * sing_sum
    }
    upd_u <- function(L, v) {
        r <- colSums(((L + v - M) * obs))
        ifelse(cnt_u > 0, -r / pmax(cnt_u, 1), 0)
    }
    upd_v <- function(L, u) {
        r <- rowSums(((L + matrix(u, nrow(M), ncol(M), byrow = TRUE) - M) * obs))
        ifelse(cnt_v > 0, -r / pmax(cnt_v, 1), 0)
    }
    # effects with L = 0, and the penalty that sets L to zero
    L <- matrix(0, nrow(M), ncol(M))
    u <- numeric(ncol(M)); v <- numeric(nrow(M))
    obj <- objective(L, u, v, 0, 0)
    for (it in seq_len(1000)) {
        u <- upd_u(L, v)
        v <- upd_v(L, u)
        new <- objective(L, u, v, 0, 0)
        rel <- (new - obj) / obj
        if (is.finite(rel) && rel < tol && rel >= 0) break
        obj <- new
    }
    E <- L + matrix(u, nrow(M), ncol(M), byrow = TRUE) + matrix(v, nrow(M), ncol(M))
    lambda_max <- 2 * max(svd(((M - E) * obs))$d) / nobs
    lambda_L <- lam_rel * lambda_max
    if (!is.finite(lambda_L) || lambda_L <= 0) lambda_L <- 1e-8
    thr <- lambda_L * nobs / 2
    sv <- svd(L)
    obj <- objective(L, u, v, sum(sv$d), lambda_L)
    for (it in seq_len(maxit)) {
        u <- upd_u(L, v)
        v <- upd_v(L, u)
        E <- L + matrix(u, nrow(M), ncol(M), byrow = TRUE) + matrix(v, nrow(M), ncol(M))
        sv <- svd(((M - E) * obs) + L)
        d <- pmax(sv$d - thr, 0)
        L <- sv$u %*% (d * t(sv$v))
        new <- objective(L, u, v, sum(d), lambda_L)
        rel <- (obj - new) / obj
        if (new < 1e-8) break
        if (is.finite(rel) && rel < tol && rel >= 0) break
        obj <- new
    }
    list(path = as.vector(L[, 1] + u[1] + v), w = NULL, tau = NULL, L = L, u = u, v = v)
}

# Conformal inference for the ATT (Chernozhukov, Wuthrich and Zhu 2021). Under the sharp
# null tau0 the counterfactual is re-estimated on all periods; p-values compare the
# post-period residual size with cyclic rotations of the residual series. The interval
# is obtained by inverting the test over a grid of tau0.
fam_conformal <- function(method, Y, tr, dn, pre, post, att, gap, pars,
                          alpha = 0.05, ngrid = 101, ck = function() NULL) {
    Tn <- nrow(Y)
    allp <- seq_len(Tn)
    Tp <- length(post)
    pval <- function(tau0) {
        ck()
        Y2 <- Y
        Y2[post, tr] <- Y[post, tr] - tau0
        f <- fam_fit(method, Y2, tr, dn, allp, post, pars$lambda,
                     pars$nfac, pars$alpha, pars$ndraw,
                     cov = pars$cov, covw = if (is.null(pars$covw)) 1 else pars$covw, acv = isTRUE(pars$acv), enfix = pars$enfix)
        u <- as.vector(Y2[, tr]) - f$path
        S <- vapply(0:(Tn - 1), function(k) {
            v <- u[((seq_len(Tn) - 1 + k) %% Tn) + 1]
            mean(abs(v[(Tn - Tp + 1):Tn]))
        }, 0)
        mean(S >= S[1])
    }
    p0 <- pval(0)
    w <- max(3 * stats::sd(gap), 2 * abs(att))
    if (!is.finite(w) || w <= 0) w <- 1
    grid <- seq(att - w, att + w, length.out = ngrid)
    pv <- vapply(grid, pval, 0)
    keep <- grid[pv >= alpha]
    list(p0 = p0,
         lower = if (length(keep)) min(keep) else NA_real_,
         upper = if (length(keep)) max(keep) else NA_real_,
         grid = grid, pv = pv)
}
