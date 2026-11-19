# Ridge-augmented synthetic control engine for a single treatment period.
#
# This file is a port of the single-period part of the R package 'augsynth'
# (https://github.com/ebenmichael/augsynth, Eli Ben-Michael, MIT licence):
# fit_synth.R, ridge.R, ridge_lambda.R, parts of augsynth.R, format.R and inference.R.
# The algorithms are kept line for line so that results match the package: osqp-based
# simplex weights, ridge augmentation with cross-validated lambda, auxiliary covariates,
# unit fixed effects, and conformal / jackknife+ / jackknife inference. Only the data
# plumbing (dplyr, Formula) is replaced by base R, and function names are prefixed `asyn_`.
# Reference: Ben-Michael, Feller and Rothstein (2021), JASA 116(536), 1789-1803.

# ---- synthetic control weights (fit_synth.R) ------------------------------------------

asyn_make_V <- function(t0, V) {
    if (is.null(V)) {
        V <- diag(rep(1, t0))
    } else if (is.vector(V)) {
        if (length(V) != t0) stop(paste("`V` must be a vector with", t0, "elements or a", t0, "x", t0, "matrix"))
        V <- diag(V)
    } else if (ncol(V) == 1 & nrow(V) == t0) {
        V <- diag(c(V))
    } else if (ncol(V) == t0 & nrow(V) == 1) {
        V <- diag(c(V))
    } else if (nrow(V) == t0) {
    } else {
        stop(paste("`V` must be a vector with", t0, "elements or a", t0, "x", t0, "matrix"))
    }
    V
}

asyn_synth_qp <- function(X1, X0, V) {
    if (!requireNamespace("osqp", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE))
        stop("The R packages 'osqp' and 'Matrix' are required. Install them with install.packages(c('osqp', 'Matrix')).")
    n0 <- nrow(X0)
    t0 <- ncol(X0)
    Pmat <- Matrix::bdiag(Matrix::Matrix(0, n0, n0, sparse = TRUE), V)
    qvec <- c(numeric(n0), -t(V) %*% X1)
    Amat <- Matrix::rbind2(
        Matrix::rbind2(
            Matrix::Matrix(c(rep(1, n0), rep(0, t0)), nrow = 1, sparse = TRUE),
            Matrix::cbind2(Matrix::Diagonal(n0), Matrix::Matrix(0, n0, t0, sparse = TRUE))
        ),
        Matrix::cbind2(-t(X0), Matrix::Diagonal(t0))
    )
    l <- c(1, numeric(n0), numeric(t0))
    u <- c(1, rep(1, n0), numeric(t0))
    settings <- osqp::osqpSettings(verbose = FALSE, eps_rel = 1e-8, eps_abs = 1e-8)
    sol <- osqp::solve_osqp(P = Pmat, q = qvec, A = Amat, l = l, u = u, pars = settings)
    sol$x[1:n0]
}

asyn_fit_synth <- function(synth_data, V = NULL) {
    t0 <- dim(synth_data$Z0)[1]
    V <- asyn_make_V(t0, V)
    weights <- asyn_synth_qp(synth_data$X1, t(synth_data$X0), V)
    l2_imbalance <- sqrt(sum((synth_data$Z0 %*% weights - synth_data$Z1)^2))
    uni_w <- matrix(1 / ncol(synth_data$Z0), nrow = ncol(synth_data$Z0), ncol = 1)
    unif_l2_imbalance <- sqrt(sum((synth_data$Z0 %*% uni_w - synth_data$Z1)^2))
    list(weights = weights, l2_imbalance = l2_imbalance,
         scaled_l2_imbalance = l2_imbalance / unif_l2_imbalance)
}

# ---- choosing the ridge penalty (ridge_lambda.R) -----------------------------------------

asyn_get_lambda_errors <- function(lambdas, X_c, X_t, synth_data, trt, holdout_length = 1, scm = TRUE) {
    errors <- matrix(0, nrow = ncol(X_c) - holdout_length, ncol = length(lambdas))
    for (i in 1:(ncol(X_c) - holdout_length)) {
        X_0 <- X_c[, -(i:(i + holdout_length - 1))]
        X_1 <- matrix(X_t[-(i:(i + holdout_length - 1))])
        X_0v <- X_c[, i:(i + holdout_length - 1)]
        X_1v <- matrix(X_t[i:(i + holdout_length - 1)], ncol = 1)
        new_synth_data <- synth_data
        new_synth_data$Z1 <- X_1
        new_synth_data$X1 <- X_1
        new_synth_data$Z0 <- t(X_0)
        new_synth_data$X0 <- t(X_0)
        if (scm) syn <- asyn_fit_synth(new_synth_data)$weights
        else syn <- rep(1 / sum(trt == 0), sum(trt == 0))
        for (j in 1:length(lambdas)) {
            ridge_weights <- t(X_1 - t(X_0) %*% syn) %*% solve(t(X_0) %*% X_0 + lambdas[j] * diag(ncol(X_0))) %*% t(X_0)
            aug_weights <- syn + t(ridge_weights)
            error <- X_1v - t(X_0v) %*% aug_weights
            errors[i, j] <- sum(error^2)
        }
    }
    list(lambda_errors = apply(errors, 2, mean),
         lambda_errors_se = apply(errors, 2, function(x) stats::sd(x) / sqrt(length(x))))
}

asyn_get_lambda_max <- function(X_c) svd(X_c)$d[1]^2

asyn_create_lambda_list <- function(lambda_max, lambda_min_ratio, n_lambda) {
    scaler <- (lambda_min_ratio)^(1 / n_lambda)
    lambda_max * (scaler^(seq(0:n_lambda) - 1))
}

asyn_choose_lambda <- function(lambdas, lambda_errors, lambda_errors_se, min_1se) {
    min_idx <- which.min(lambda_errors)
    min_error <- lambda_errors[min_idx]
    min_se <- lambda_errors_se[min_idx]
    lambda_min <- lambdas[min_idx]
    lambda_1se <- max(lambdas[lambda_errors <= min_error + min_se])
    if (min_1se) lambda_1se else lambda_min
}

asyn_cv_lambda <- function(X_c, X_1, synth_data, trt, holdout_length, scm,
                           lambda_max, lambda_min_ratio, n_lambda, min_1se) {
    if (is.null(lambda_max)) lambda_max <- asyn_get_lambda_max(X_c)
    lambdas <- asyn_create_lambda_list(lambda_max, lambda_min_ratio, n_lambda)
    lambda_out <- asyn_get_lambda_errors(lambdas, X_c, X_1, synth_data, trt, holdout_length, scm)
    lambda <- asyn_choose_lambda(lambdas, lambda_out$lambda_errors, lambda_out$lambda_errors_se, min_1se)
    list(lambda = lambda, lambda_errors = lambda_out$lambda_errors,
         lambda_errors_se = lambda_out$lambda_errors_se, lambdas = lambdas)
}

# ---- ridge-augmented weights (ridge.R) -----------------------------------------------------

asyn_fit_ridgeaug_inner <- function(X_c, X_1, trt, synth_data, lambda, ridge, scm,
                                    lambda_min_ratio, n_lambda, lambda_max, holdout_length, min_1se) {
    lambda_errors <- NULL
    lambda_errors_se <- NULL
    lambdas <- NULL
    if (scm) syn <- asyn_fit_synth(synth_data)$weights
    else syn <- rep(1 / sum(trt == 0), sum(trt == 0))
    if (ridge) {
        if (is.null(lambda)) {
            cv_out <- asyn_cv_lambda(X_c, X_1, synth_data, trt, holdout_length, scm,
                                     lambda_max, lambda_min_ratio, n_lambda, min_1se)
            lambda <- cv_out$lambda
            lambda_errors <- cv_out$lambda_errors
            lambda_errors_se <- cv_out$lambda_errors_se
            lambdas <- cv_out$lambdas
        }
        ridge_w <- t(t(X_1) - t(X_c) %*% syn) %*% solve(t(X_c) %*% X_c + lambda * diag(ncol(X_c))) %*% t(X_c)
    } else {
        ridge_w <- matrix(0, ncol = sum(trt == 0), nrow = 1)
    }
    weights <- syn + t(ridge_w)
    list(weights = weights, synw = syn, lambda = lambda, lambdas = lambdas,
         lambda_errors = lambda_errors, lambda_errors_se = lambda_errors_se)
}

asyn_fit_ridgeaug <- function(wide_data, synth_data, Z = NULL, lambda = NULL, ridge = TRUE, scm = TRUE,
                              lambda_min_ratio = 1e-8, n_lambda = 20, lambda_max = NULL,
                              holdout_length = 1, min_1se = TRUE, V = NULL,
                              residualize = FALSE, ...) {
    extra_params <- list(...)
    if (length(extra_params) > 0)
        warning("Unused parameters in using ridge augmented weights: ", paste(names(extra_params), collapse = ", "))
    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt

    X_cent <- apply(X, 2, function(x) x - mean(x[trt == 0]))
    X_c <- X_cent[trt == 0, , drop = FALSE]
    X_1 <- matrix(colMeans(X_cent[trt == 1, , drop = FALSE]), nrow = 1)
    y_cent <- apply(y, 2, function(x) x - mean(x[trt == 0]))
    y_c <- y_cent[trt == 0, , drop = FALSE]

    t0 <- ncol(X_c)
    V <- asyn_make_V(t0, V)
    X_c <- X_c %*% V
    X_1 <- X_1 %*% V
    new_synth_data <- synth_data

    if (!is.null(Z)) {
        Z_cent <- apply(Z, 2, function(x) x - mean(x[trt == 0]))
        Z_c <- Z_cent[trt == 0, , drop = FALSE]
        Z_1 <- matrix(colMeans(Z_cent[trt == 1, , drop = FALSE]), nrow = 1)
        if (residualize) {
            Xc_hat <- Z_c %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c) %*% X_c
            X1_hat <- Z_1 %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c) %*% X_c
            res_t <- X_1 - X1_hat
            res_c <- X_c - Xc_hat
            X_c <- res_c
            X_1 <- res_t
            X_cent[trt == 0, ] <- res_c
            X_cent[trt == 1, ] <- res_t
            new_synth_data$Z1 <- t(res_t)
            new_synth_data$X1 <- t(res_t)
            new_synth_data$Z0 <- t(res_c)
            new_synth_data$X0 <- t(res_c)
        } else {
            sdz <- apply(Z_c, 2, stats::sd)
            sdx <- stats::sd(X_c)
            Z_c <- sdx * t(t(Z_c) / sdz)
            Z_1 <- sdx * Z_1 / sdz
            X_c <- cbind(X_c, Z_c)
            X_1 <- cbind(X_1, Z_1)
            new_synth_data$Z1 <- t(X_1)
            new_synth_data$X1 <- t(X_1)
            new_synth_data$Z0 <- t(X_c)
            new_synth_data$X0 <- t(X_c)
            V <- diag(ncol(X_c))
        }
    } else {
        new_synth_data$Z1 <- t(X_1)
        new_synth_data$X1 <- t(X_1)
        new_synth_data$Z0 <- t(X_c)
        new_synth_data$X0 <- t(X_c)
    }
    out <- asyn_fit_ridgeaug_inner(X_c, X_1, trt, new_synth_data, lambda, ridge, scm,
                                   lambda_min_ratio, n_lambda, lambda_max, holdout_length, min_1se)
    weights <- out$weights
    synw <- out$synw
    lambda <- out$lambda
    lambdas <- out$lambdas
    lambda_errors <- out$lambda_errors
    lambda_errors_se <- out$lambda_errors_se

    if (!is.null(Z)) {
        if (residualize) {
            no_cov_weights <- weights
            ridge_w <- t(t(Z_1) - t(Z_c) %*% weights) %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c)
            weights <- weights + t(ridge_w)
        } else {
            no_cov_weights <- NULL
        }
    }

    l2_imbalance <- sqrt(sum((synth_data$X0 %*% weights - synth_data$X1)^2))
    uni_w <- matrix(1 / ncol(synth_data$X0), nrow = ncol(synth_data$X0), ncol = 1)
    unif_l2_imbalance <- sqrt(sum((synth_data$X0 %*% uni_w - synth_data$X1)^2))
    scaled_l2_imabalance <- l2_imbalance / unif_l2_imbalance

    mhat <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    ridge_mhat <- mhat
    if (!is.null(Z)) {
        if (residualize) {
            ridge_mhat <- ridge_mhat + Z_cent %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c) %*% y_c
            yc_hat <- ridge_mhat[trt == 0, , drop = FALSE]
            y_c <- y_c - yc_hat
        } else {
            X_cent <- cbind(X_cent, Z_cent)
        }
    }
    if (ridge) {
        ridge_mhat <- ridge_mhat + X_cent %*% solve(t(X_c) %*% X_c + lambda * diag(ncol(X_c))) %*% t(X_c) %*% y_c
    }
    output <- list(weights = weights, l2_imbalance = l2_imbalance, scaled_l2_imbalance = scaled_l2_imabalance,
                   mhat = mhat, lambda = lambda, ridge_mhat = ridge_mhat, synw = synw, lambdas = lambdas,
                   lambda_errors = lambda_errors, lambda_errors_se = lambda_errors_se)
    if (!is.null(Z)) {
        output$no_cov_weights <- no_cov_weights
        z_l2_imbalance <- sqrt(sum((t(Z_c) %*% weights - t(Z_1))^2))
        z_unif_l2_imbalance <- sqrt(sum((t(Z_c) %*% uni_w - t(Z_1))^2))
        output$covariate_l2_imbalance <- z_l2_imbalance
        output$scaled_covariate_l2_imbalance <- z_l2_imbalance / z_unif_l2_imbalance
    }
    output
}

# ---- data formatting (format.R, base R) ----------------------------------------------------

# Y: T x N outcome matrix with unit names as colnames; units are sorted by name, as in augsynth.
asyn_format <- function(Y, times, treated, t_int) {
    units <- sort(colnames(Y))
    Y <- Y[, units, drop = FALSE]
    pre <- which(times < t_int)
    post <- which(times >= t_int)
    X <- t(Y[pre, , drop = FALSE]); colnames(X) <- as.character(times[pre])
    y <- t(Y[post, , drop = FALSE]); colnames(y) <- as.character(times[post])
    trt <- as.numeric(units %in% treated)
    list(X = X, trt = trt, y = y, units = units)
}

asyn_format_synth <- function(X, trt, y) {
    sd <- list()
    sd$Z0 <- t(X[trt == 0, , drop = FALSE])
    sd$Z1 <- as.matrix(colMeans(X[trt == 1, , drop = FALSE]), ncol = 1)
    sd$Y0plot <- t(cbind(X[trt == 0, , drop = FALSE], y[trt == 0, , drop = FALSE]))
    sd$Y1plot <- as.matrix(colMeans(cbind(X[trt == 1, , drop = FALSE], y[trt == 1, , drop = FALSE])), ncol = 1)
    sd$X0 <- sd$Z0
    sd$X1 <- sd$Z1
    sd
}

asyn_demean <- function(wide_data, synth_data) {
    means <- rowMeans(wide_data$X)
    new_wide_data <- list()
    new_X <- wide_data$X - means
    trt <- wide_data$trt
    new_wide_data$X <- new_X
    new_wide_data$y <- wide_data$y - means
    new_wide_data$trt <- trt
    new_synth_data <- list()
    new_synth_data$X0 <- t(new_X[trt == 0, , drop = FALSE])
    new_synth_data$Z0 <- new_synth_data$X0
    new_synth_data$X1 <- as.matrix((colMeans(new_X[trt == 1, , drop = FALSE])), ncol = 1)
    new_synth_data$Z1 <- new_synth_data$X1
    mhat <- replicate(ncol(wide_data$X) + ncol(wide_data$y), means)
    list(wide = new_wide_data, synth_data = new_synth_data, mhat = mhat)
}

# ---- fitting and prediction (augsynth.R) -----------------------------------------------------

asyn_fit_internal <- function(wide, synth_data, Z, progfunc, scm, fixedeff, V = NULL, ...) {
    n <- nrow(wide$X)
    t0 <- ncol(wide$X)
    ttot <- t0 + ncol(wide$y)
    if (fixedeff) {
        demeaned <- asyn_demean(wide, synth_data)
        fit_wide <- demeaned$wide
        fit_synth_data <- demeaned$synth_data
        mhat <- demeaned$mhat
    } else {
        fit_wide <- wide
        fit_synth_data <- synth_data
        mhat <- matrix(0, n, ttot)
    }
    if (is.null(progfunc)) progfunc <- "none"
    progfunc <- tolower(progfunc)
    if (progfunc == "ridge") {
        augsynth <- do.call(asyn_fit_ridgeaug,
                            list(wide_data = fit_wide, synth_data = fit_synth_data,
                                 Z = Z, V = V, ridge = TRUE, scm = scm, ...))
    } else if (progfunc == "none") {
        if (!scm) stop("Cannot run with progfunc='none' AND no SCM weights")
        augsynth <- do.call(asyn_fit_ridgeaug,
                            c(list(wide_data = fit_wide, synth_data = fit_synth_data,
                                   Z = Z, V = V, ridge = FALSE, scm = scm, ...)))
    } else {
        stop("progfunc must be one of 'ridge', 'none'")
    }
    augsynth$mhat <- mhat + cbind(matrix(0, nrow = n, ncol = t0), augsynth$mhat)
    augsynth$data <- wide
    augsynth$data$Z <- Z
    augsynth$data$synth_data <- synth_data
    augsynth$progfunc <- progfunc
    augsynth$scm <- scm
    augsynth$fixedeff <- fixedeff
    augsynth$extra_args <- list(...)
    if (progfunc == "ridge") augsynth$extra_args$lambda <- augsynth$lambda
    class(augsynth) <- "asyn_fit"
    augsynth
}

# Fit from a T x N panel: the one entry point used by the jamovi analysis.
asyn_fit <- function(Y, times, treated, t_int, Z = NULL, progfunc = "ridge", scm = TRUE,
                     fixedeff = FALSE, ...) {
    wide <- asyn_format(Y, times, treated, t_int)
    synth_data <- asyn_format_synth(wide$X, wide$trt, wide$y)
    if (!is.null(Z)) Z <- Z[wide$units, , drop = FALSE]
    fit <- asyn_fit_internal(wide, synth_data, Z, tolower(progfunc), scm, fixedeff, ...)
    fit$data$time <- times
    fit$units <- wide$units
    fit$weights <- matrix(fit$weights)
    rownames(fit$weights) <- wide$units[wide$trt == 0]
    fit$t_int <- t_int
    fit
}

# covariates: mean over the pre-treatment periods (omitting NAs), as augsynth does by default
asyn_covariates <- function(cov, times, t_int, units) {
    pre <- which(times < t_int)
    K <- dim(cov)[3]
    Z <- matrix(apply(cov[pre, units, , drop = FALSE], c(2, 3), mean, na.rm = TRUE),
                length(units), K, dimnames = list(units, dimnames(cov)[[3]]))
    if (any(!is.finite(Z))) stop("Some units are missing all covariate data in the pre-treatment periods.")
    zs <- apply(Z, 2, stats::sd)
    if (any(zs == 0))
        stop(paste("The following covariates have no variation across units:",
                   paste(colnames(Z)[zs == 0], collapse = ", ")))
    Z
}

asyn_predict <- function(object, att = FALSE) {
    X <- object$data$X
    y <- object$data$y
    comb <- cbind(X, y)
    trt <- object$data$trt
    mhat <- object$mhat
    m1 <- colMeans(mhat[trt == 1, , drop = FALSE])
    resid <- (comb[trt == 0, , drop = FALSE] - mhat[trt == 0, drop = FALSE])
    y0 <- m1 + t(resid) %*% object$weights
    if (att) {
        colMeans(comb[trt == 1, , drop = FALSE]) - c(y0)
    } else {
        rnames <- rownames(y0)
        y0_vec <- c(y0)
        names(y0_vec) <- rnames
        y0_vec
    }
}

# ---- inference (inference.R) -------------------------------------------------------------------

asyn_drop_time_t <- function(wide_data, Z, t_drop) {
    new_wide_data <- list()
    new_wide_data$trt <- wide_data$trt
    new_wide_data$X <- wide_data$X[, -t_drop, drop = FALSE]
    new_wide_data$y <- cbind(wide_data$X[, t_drop, drop = FALSE], wide_data$y)
    X0 <- new_wide_data$X[new_wide_data$trt == 0, , drop = FALSE]
    x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1, , drop = FALSE]), ncol = 1)
    new_synth_data <- list()
    new_synth_data$Z0 <- t(X0)
    new_synth_data$X0 <- t(X0)
    new_synth_data$Z1 <- x1
    new_synth_data$X1 <- x1
    list(wide_data = new_wide_data, synth_data = new_synth_data, Z = Z)
}

asyn_drop_unit_i <- function(wide_data, Z, i) {
    new_wide_data <- list()
    new_wide_data$trt <- wide_data$trt[-i]
    new_wide_data$X <- wide_data$X[-i, , drop = FALSE]
    new_wide_data$y <- wide_data$y[-i, , drop = FALSE]
    X0 <- new_wide_data$X[new_wide_data$trt == 0, , drop = FALSE]
    x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1, , drop = FALSE]), ncol = 1)
    new_synth_data <- list()
    new_synth_data$Z0 <- t(X0)
    new_synth_data$X0 <- t(X0)
    new_synth_data$Z1 <- x1
    new_synth_data$X1 <- x1
    new_Z <- if (!is.null(Z)) Z[-i, , drop = FALSE] else NULL
    list(wide_data = new_wide_data, synth_data = new_synth_data, Z = new_Z)
}

asyn_refit <- function(ascm, wide, synth_data, Z) {
    do.call(asyn_fit_internal,
            c(list(wide = wide, synth_data = synth_data, Z = Z, progfunc = ascm$progfunc,
                   scm = ascm$scm, fixedeff = ascm$fixedeff), ascm$extra_args))
}

asyn_time_jackknife_plus <- function(ascm, alpha = 0.05, conservative = FALSE, checkpoint = NULL) {
    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    Z <- wide_data$Z
    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    jack_ests <- lapply(1:t0, function(tdrop) {
        if (!is.null(checkpoint)) checkpoint()
        new_data <- asyn_drop_time_t(wide_data, Z, tdrop)
        new_ascm <- asyn_refit(ascm, new_data$wide_data, new_data$synth_data, new_data$Z)
        est <- asyn_predict(new_ascm, att = FALSE)[(t0 + 1):t_final]
        est <- c(est, mean(est))
        err <- c(colMeans(wide_data$X[wide_data$trt == 1, tdrop, drop = FALSE]) -
                     asyn_predict(new_ascm, att = FALSE)[t0])
        list(err, rbind(est + abs(err), est - abs(err), est + err, est))
    })
    held_out_errs <- vapply(jack_ests, `[[`, numeric(1), 1)
    jack_dist <- vapply(jack_ests, `[[`, matrix(0, nrow = 4, ncol = tpost + 1), 2)
    out <- list()
    att <- asyn_predict(ascm, att = TRUE)
    out$att <- c(att, mean(att[(t0 + 1):t_final]))
    out$heldout_att <- c(held_out_errs, att[(t0 + 1):t_final], mean(att[(t0 + 1):t_final]))
    if (conservative) {
        qerr <- stats::quantile(abs(held_out_errs), 1 - alpha)
        out$lb <- c(rep(NA, t0), apply(jack_dist[4, , ], 1, min) - qerr)
        out$ub <- c(rep(NA, t0), apply(jack_dist[4, , ], 1, max) + qerr)
    } else {
        out$lb <- c(rep(NA, t0), apply(jack_dist[2, , ], 1, stats::quantile, alpha / 2))
        out$ub <- c(rep(NA, t0), apply(jack_dist[1, , ], 1, stats::quantile, 1 - alpha / 2))
    }
    y1 <- asyn_predict(ascm, att = FALSE) + att
    y1 <- c(y1, mean(y1[(t0 + 1):t_final]))
    shifted_lb <- y1 - out$ub
    shifted_ub <- y1 - out$lb
    out$lb <- shifted_lb
    out$ub <- shifted_ub
    out$alpha <- alpha
    out
}

asyn_jackknife_se <- function(ascm, checkpoint = NULL) {
    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    Z <- wide_data$Z
    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    nnz_weights <- numeric(n)
    nnz_weights[wide_data$trt == 0] <- round(ascm$weights, 3) != 0
    if (sum(wide_data$trt) > 1) nnz_weights[wide_data$trt == 1] <- 1
    trt_idxs <- (1:n)[as.logical(nnz_weights)]
    ests <- vapply(trt_idxs, function(i) {
        if (!is.null(checkpoint)) checkpoint()
        new_data <- asyn_drop_unit_i(wide_data, Z, i)
        new_ascm <- asyn_refit(ascm, new_data$wide_data, new_data$synth_data, new_data$Z)
        est <- asyn_predict(new_ascm, att = TRUE)[(t0 + 1):t_final]
        c(est, mean(est))
    }, numeric(tpost + 1))
    ests <- matrix(ests, nrow = tpost + 1, ncol = length(trt_idxs))
    se2 <- apply(ests, 1, function(x) (n - 1) / n * sum((x - mean(x, na.rm = TRUE))^2))
    se <- sqrt(se2)
    out <- list()
    att <- asyn_predict(ascm, att = TRUE)
    out$att <- c(att, mean(att[(t0 + 1):t_final]))
    out$se <- c(rep(NA, t0), se)
    out
}

asyn_permute_test_stats <- function(wide_data, ascm, h0, post_length, type, q, ns, stat_func) {
    new_wide_data <- wide_data
    t0 <- ncol(wide_data$X) - post_length
    tpost <- t0 + post_length
    new_wide_data$X[wide_data$trt == 1, (t0 + 1):tpost] <- new_wide_data$X[wide_data$trt == 1, (t0 + 1):tpost] - h0
    X0 <- new_wide_data$X[new_wide_data$trt == 0, , drop = FALSE]
    x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1, , drop = FALSE]), ncol = 1)
    new_synth_data <- list()
    new_synth_data$Z0 <- t(X0)
    new_synth_data$X0 <- t(X0)
    new_synth_data$Z1 <- x1
    new_synth_data$X1 <- x1
    new_ascm <- asyn_refit(ascm, new_wide_data, new_synth_data, wide_data$Z)
    resids <- asyn_predict(new_ascm, att = TRUE)[1:tpost]
    if (is.null(stat_func)) stat_func <- function(x) (sum(abs(x)^q) / sqrt(length(x)))^(1 / q)
    if (type == "iid") {
        test_stats <- sapply(1:ns, function(x) {
            reorder <- sample(resids)
            stat_func(reorder[(t0 + 1):tpost])
        })
    } else {
        test_stats <- sapply(1:tpost, function(j) {
            reorder <- resids[(0:tpost - 1 + j) %% tpost + 1]
            stat_func(reorder[(t0 + 1):tpost])
        })
    }
    list(resids = resids, test_stats = test_stats, stat_func = stat_func)
}

asyn_permute_pval <- function(wide_data, ascm, h0, post_length, type, q, ns, stat_func) {
    t0 <- ncol(wide_data$X) - post_length
    tpost <- t0 + post_length
    out <- asyn_permute_test_stats(wide_data, ascm, h0, post_length, type, q, ns, stat_func)
    mean(out$stat_func(out$resids[(t0 + 1):tpost]) <= out$test_stats)
}

asyn_permute_ci <- function(wide_data, ascm, grid, post_length, alpha, type, q, ns, stat_func, checkpoint = NULL) {
    grid <- c(grid, 0)
    ps <- sapply(grid, function(x) {
        if (!is.null(checkpoint)) checkpoint()
        asyn_permute_pval(wide_data, ascm, x, post_length, type, q, ns, stat_func)
    })
    c(min(grid[ps >= alpha]), max(grid[ps >= alpha]), ps[grid == 0])
}

asyn_conformal_inf <- function(ascm, alpha = 0.05, stat_func = NULL, type = "iid",
                               q = 1, ns = 1000, grid_size = 50, checkpoint = NULL) {
    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    att <- asyn_predict(ascm, att = TRUE)
    post_att <- att[(t0 + 1):t_final]
    post_sd <- sqrt(mean(post_att^2))
    cis <- vapply(1:tpost, function(j) {
        new_wide_data <- wide_data
        new_wide_data$X <- cbind(wide_data$X, wide_data$y[, j, drop = TRUE])
        if (tpost > 1) new_wide_data$y <- wide_data$y[, -j, drop = FALSE]
        else new_wide_data$y <- matrix(1, nrow = n, ncol = 1)
        grid <- seq(att[t0 + j] - 2 * post_sd, att[t0 + j] + 2 * post_sd, length.out = grid_size)
        asyn_permute_ci(new_wide_data, ascm, grid, 1, alpha, type, q, ns, stat_func, checkpoint)
    }, numeric(3))
    new_wide_data <- wide_data
    new_wide_data$X <- cbind(wide_data$X, wide_data$y)
    new_wide_data$y <- matrix(1, nrow = n, ncol = 1)
    null_p <- asyn_permute_pval(new_wide_data, ascm, 0, ncol(wide_data$y), type, q, ns, stat_func)
    out <- list()
    att <- asyn_predict(ascm, att = TRUE)
    out$att <- c(att, mean(att[(t0 + 1):t_final]))
    out$lb <- c(rep(NA, t0), cis[1, ], NA)
    out$ub <- c(rep(NA, t0), cis[2, ], NA)
    out$p_val <- c(rep(NA, t0), cis[3, ], null_p)
    out$alpha <- alpha
    out
}
