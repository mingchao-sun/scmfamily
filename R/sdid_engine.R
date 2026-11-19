# Synthetic difference-in-differences engine.
#
# This file is a port of the R package 'synthdid' (https://github.com/synth-inference/synthdid),
# Copyright 2019, Stanford University, distributed under GPL (>= 2) | BSD-3-Clause.
# Authors: Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, Stefan Wager.
# The algorithms (Frank-Wolfe weight solver, sparsification, covariate solver, collapsed form,
# placebo / bootstrap / jackknife standard errors) are kept line for line so that results
# match the package; only the function names are prefixed with `sdid_` and the S4 plumbing,
# plotting and data-reshaping helpers are left out. Reference: Arkhangelsky, Athey, Hirshberg,
# Imbens and Wager (2021), American Economic Review 111(12), 4088-4118.

sdid_contract3 <- function(X, v) {
    stopifnot(length(dim(X)) == 3, dim(X)[3] == length(v))
    out <- array(0, dim = dim(X)[1:2])
    if (length(v) == 0) return(out)
    for (ii in 1:length(v)) out <- out + v[ii] * X[, , ii]
    out
}

# a Frank-Wolfe step for ||Ax - b||^2 + eta * ||x||^2 with x in the unit simplex
sdid_fw_step <- function(A, x, b, eta, alpha = NULL) {
    Ax <- A %*% x
    half.grad <- t(Ax - b) %*% A + eta * x
    i <- which.min(half.grad)
    if (!is.null(alpha)) {
        x <- x * (1 - alpha)
        x[i] <- x[i] + alpha
        return(x)
    } else {
        d.x <- -x; d.x[i] <- 1 - x[i]
        if (all(d.x == 0)) return(x)
        d.err <- A[, i] - Ax
        step <- -t(c(half.grad)) %*% d.x / (sum(d.err^2) + eta * sum(d.x^2))
        constrained.step <- min(1, max(0, step))
        return(x + constrained.step * d.x)
    }
}

# a Frank-Wolfe solver for synthetic control weights using exact line search
sdid_sc_weight_fw <- function(Y, zeta, intercept = TRUE, lambda = NULL, min.decrease = 1e-3, max.iter = 1000) {
    T0 <- ncol(Y) - 1
    N0 <- nrow(Y)
    if (is.null(lambda)) lambda <- rep(1 / T0, T0)
    if (intercept) Y <- apply(Y, 2, function(col) col - mean(col))
    t <- 0
    vals <- rep(NA, max.iter)
    A <- Y[, 1:T0]
    b <- Y[, T0 + 1]
    eta <- N0 * Re(zeta^2)
    while (t < max.iter && (t < 2 || vals[t - 1] - vals[t] > min.decrease^2)) {
        t <- t + 1
        lambda.p <- sdid_fw_step(A, lambda, b, eta)
        lambda <- lambda.p
        err <- Y[1:N0, ] %*% c(lambda, -1)
        vals[t] <- Re(zeta^2) * sum(lambda^2) + sum(err^2) / N0
    }
    list(lambda = lambda, vals = vals)
}

# Frank-Wolfe + gradient solver for lambda, omega and beta when there are covariates
sdid_sc_weight_fw_covariates <- function(Y, X = array(0, dim = c(dim(Y), 0)), zeta.lambda = 0, zeta.omega = 0,
                                         lambda.intercept = TRUE, omega.intercept = TRUE,
                                         min.decrease = 1e-3, max.iter = 1000,
                                         lambda = NULL, omega = NULL, beta = NULL,
                                         update.lambda = TRUE, update.omega = TRUE) {
    stopifnot(length(dim(Y)) == 2, length(dim(X)) == 3, all(dim(Y) == dim(X)[1:2]),
              all(is.finite(Y)), all(is.finite(X)))
    T0 <- ncol(Y) - 1
    N0 <- nrow(Y) - 1
    if (length(dim(X)) == 2) dim(X) <- c(dim(X), 1)
    if (is.null(lambda)) lambda <- rep(1 / T0, T0)
    if (is.null(omega)) omega <- rep(1 / N0, N0)
    if (is.null(beta)) beta <- rep(0, dim(X)[3])

    update.weights <- function(Y, lambda, omega) {
        Y.lambda <- if (lambda.intercept) apply(Y[1:N0, ], 2, function(row) row - mean(row)) else Y[1:N0, ]
        if (update.lambda) lambda <- sdid_fw_step(Y.lambda[, 1:T0], lambda, Y.lambda[, T0 + 1], N0 * Re(zeta.lambda^2))
        err.lambda <- Y.lambda %*% c(lambda, -1)

        Y.omega <- if (omega.intercept) apply(t(Y[, 1:T0]), 2, function(row) row - mean(row)) else t(Y[, 1:T0])
        if (update.omega) omega <- sdid_fw_step(Y.omega[, 1:N0], omega, Y.omega[, N0 + 1], T0 * Re(zeta.omega^2))
        err.omega <- Y.omega %*% c(omega, -1)

        val <- Re(zeta.omega^2) * sum(omega^2) + Re(zeta.lambda^2) * sum(lambda^2) +
            sum(err.omega^2) / T0 + sum(err.lambda^2) / N0
        list(val = val, lambda = lambda, omega = omega, err.lambda = err.lambda, err.omega = err.omega)
    }

    vals <- rep(NA, max.iter)
    t <- 0
    Y.beta <- Y - sdid_contract3(X, beta)
    weights <- update.weights(Y.beta, lambda, omega)
    while (t < max.iter && (t < 2 || abs(vals[t - 1] - vals[t]) > min.decrease^2)) {
        t <- t + 1
        grad.beta <- -if (dim(X)[3] == 0) c() else {
            apply(X, 3, function(Xi) {
                t(weights$err.lambda) %*% Xi[1:N0, ] %*% c(weights$lambda, -1) / N0 +
                    t(weights$err.omega) %*% t(Xi[, 1:T0]) %*% c(weights$omega, -1) / T0
            })
        }
        alpha <- 1 / t
        beta <- beta - alpha * grad.beta
        Y.beta <- Y - sdid_contract3(X, beta)
        weights <- update.weights(Y.beta, weights$lambda, weights$omega)
        vals[t] <- weights$val
    }
    list(lambda = weights$lambda, omega = weights$omega, beta = beta, vals = vals)
}

sdid_sparsify <- function(v) { v[v <= max(v) / 4] <- 0; v / sum(v) }

# collapse Y to an (N0+1) x (T0+1) matrix by averaging the treated rows and post columns
sdid_collapsed_form <- function(Y, N0, T0) {
    N <- nrow(Y); T <- ncol(Y)
    rbind(cbind(Y[1:N0, 1:T0, drop = FALSE], rowMeans(Y[1:N0, (T0 + 1):T, drop = FALSE])),
          cbind(t(colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE])), mean(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE])))
}

sdid_pairwise_sum_decreasing <- function(x, y) {
    na.x <- is.na(x)
    na.y <- is.na(y)
    x[is.na(x)] <- min(x[!na.x])
    y[is.na(y)] <- min(y[!na.y])
    pairwise.sum <- x + y
    pairwise.sum[na.x & na.y] <- NA
    pairwise.sum
}

# Y: N x T with the N0 control rows first and the treated rows last; first T0 columns pre-treatment
sdid_estimate <- function(Y, N0, T0, X = array(dim = c(dim(Y), 0)),
                          noise.level = stats::sd(apply(Y[1:N0, 1:T0], 1, diff)),
                          eta.omega = ((nrow(Y) - N0) * (ncol(Y) - T0))^(1/4), eta.lambda = 1e-6,
                          zeta.omega = eta.omega * noise.level, zeta.lambda = eta.lambda * noise.level,
                          omega.intercept = TRUE, lambda.intercept = TRUE,
                          weights = list(omega = NULL, lambda = NULL),
                          update.omega = is.null(weights$omega), update.lambda = is.null(weights$lambda),
                          min.decrease = 1e-5 * noise.level, max.iter = 1e4,
                          sparsify = sdid_sparsify, max.iter.pre.sparsify = 100) {
    stopifnot(nrow(Y) > N0, ncol(Y) > T0, length(dim(X)) %in% c(2, 3), dim(X)[1:2] == dim(Y), is.list(weights),
              is.null(weights$lambda) || length(weights$lambda) == T0,
              is.null(weights$omega) || length(weights$omega) == N0,
              !is.null(weights$lambda) || update.lambda, !is.null(weights$omega) || update.omega)
    if (length(dim(X)) == 2) dim(X) <- c(dim(X), 1)
    if (is.null(sparsify)) max.iter.pre.sparsify <- max.iter
    N1 <- nrow(Y) - N0
    T1 <- ncol(Y) - T0

    if (dim(X)[3] == 0) {
        weights$vals <- NULL
        weights$lambda.vals <- NULL
        weights$omega.vals <- NULL
        if (update.lambda) {
            Yc <- sdid_collapsed_form(Y, N0, T0)
            lambda.opt <- sdid_sc_weight_fw(Yc[1:N0, ], zeta = zeta.lambda, intercept = lambda.intercept,
                                            lambda = weights$lambda, min.decrease = min.decrease,
                                            max.iter = max.iter.pre.sparsify)
            if (!is.null(sparsify)) {
                lambda.opt <- sdid_sc_weight_fw(Yc[1:N0, ], zeta = zeta.lambda, intercept = lambda.intercept,
                                                lambda = sparsify(lambda.opt$lambda),
                                                min.decrease = min.decrease, max.iter = max.iter)
            }
            weights$lambda <- lambda.opt$lambda
            weights$lambda.vals <- lambda.opt$vals
            weights$vals <- lambda.opt$vals
        }
        if (update.omega) {
            Yc <- sdid_collapsed_form(Y, N0, T0)
            omega.opt <- sdid_sc_weight_fw(t(Yc[, 1:T0]), zeta = zeta.omega, intercept = omega.intercept,
                                           lambda = weights$omega, min.decrease = min.decrease,
                                           max.iter = max.iter.pre.sparsify)
            if (!is.null(sparsify)) {
                omega.opt <- sdid_sc_weight_fw(t(Yc[, 1:T0]), zeta = zeta.omega, intercept = omega.intercept,
                                               lambda = sparsify(omega.opt$lambda),
                                               min.decrease = min.decrease, max.iter = max.iter)
            }
            weights$omega <- omega.opt$lambda
            weights$omega.vals <- omega.opt$vals
            if (is.null(weights$vals)) weights$vals <- omega.opt$vals
            else weights$vals <- sdid_pairwise_sum_decreasing(weights$vals, omega.opt$vals)
        }
    } else {
        Yc <- sdid_collapsed_form(Y, N0, T0)
        Xc <- apply(X, 3, function(Xi) sdid_collapsed_form(Xi, N0, T0))
        dim(Xc) <- c(dim(Yc), dim(X)[3])
        weights <- sdid_sc_weight_fw_covariates(Yc, Xc, zeta.lambda = zeta.lambda, zeta.omega = zeta.omega,
                                                lambda.intercept = lambda.intercept, omega.intercept = omega.intercept,
                                                min.decrease = min.decrease, max.iter = max.iter,
                                                lambda = weights$lambda, omega = weights$omega,
                                                update.lambda = update.lambda, update.omega = update.omega)
    }

    X.beta <- sdid_contract3(X, weights$beta)
    estimate <- t(c(-weights$omega, rep(1 / N1, N1))) %*% (Y - X.beta) %*% c(-weights$lambda, rep(1 / T1, T1))

    class(estimate) <- "sdid_estimate"
    attr(estimate, "estimator") <- "sdid_estimate"
    attr(estimate, "weights") <- weights
    attr(estimate, "setup") <- list(Y = Y, X = X, N0 = N0, T0 = T0)
    attr(estimate, "opts") <- list(zeta.omega = zeta.omega, zeta.lambda = zeta.lambda,
                                   omega.intercept = omega.intercept, lambda.intercept = lambda.intercept,
                                   update.omega = update.omega, update.lambda = update.lambda,
                                   min.decrease = min.decrease, max.iter = max.iter)
    estimate
}

# synthetic control estimator (infinitesimal ridge, no time weights)
sdid_sc_estimate <- function(Y, N0, T0, eta.omega = 1e-6, ...) {
    estimate <- sdid_estimate(Y, N0, T0, eta.omega = eta.omega,
                              weights = list(lambda = rep(0, T0)), omega.intercept = FALSE, ...)
    attr(estimate, "estimator") <- "sc_estimate"
    estimate
}

# difference-in-differences estimator (uniform weights)
sdid_did_estimate <- function(Y, N0, T0, ...) {
    estimate <- sdid_estimate(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0), omega = rep(1 / N0, N0)), ...)
    attr(estimate, "estimator") <- "did_estimate"
    estimate
}

# the effect curve that is averaged to produce the estimate (one value per post period)
sdid_effect_curve <- function(estimate) {
    setup <- attr(estimate, "setup")
    weights <- attr(estimate, "weights")
    X.beta <- sdid_contract3(setup$X, weights$beta)
    N1 <- nrow(setup$Y) - setup$N0
    T1 <- ncol(setup$Y) - setup$T0
    tau.sc <- t(c(-weights$omega, rep(1 / N1, N1))) %*% (setup$Y - X.beta)
    tau.sc[setup$T0 + (1:T1)] - c(tau.sc[1:setup$T0] %*% weights$lambda)
}

sdid_sum_normalize <- function(x) {
    if (sum(x) != 0) x / sum(x) else rep(1 / length(x), length(x))
}

# ---- standard errors (Algorithms 2-4 of Arkhangelsky et al.) ----------------------------

sdid_bootstrap_sample <- function(estimate, replications, checkpoint = NULL) {
    setup <- attr(estimate, "setup")
    opts <- attr(estimate, "opts")
    weights <- attr(estimate, "weights")
    if (setup$N0 == nrow(setup$Y) - 1) return(NA)
    theta <- function(ind) {
        if (all(ind <= setup$N0) || all(ind > setup$N0)) NA
        else {
            weights.boot <- weights
            weights.boot$omega <- sdid_sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])
            do.call(sdid_estimate, c(list(Y = setup$Y[sort(ind), ], N0 = sum(ind <= setup$N0), T0 = setup$T0,
                                          X = setup$X[sort(ind), , ], weights = weights.boot), opts))
        }
    }
    bootstrap.estimates <- rep(NA, replications)
    count <- 0
    while (count < replications) {
        if (!is.null(checkpoint)) checkpoint()
        bootstrap.estimates[count + 1] <- theta(sample(1:nrow(setup$Y), replace = TRUE))
        if (!is.na(bootstrap.estimates[count + 1])) count <- count + 1
    }
    bootstrap.estimates
}

sdid_bootstrap_se <- function(estimate, replications, checkpoint = NULL)
    sqrt((replications - 1) / replications) * stats::sd(sdid_bootstrap_sample(estimate, replications, checkpoint))

sdid_jackknife <- function(x, theta, checkpoint = NULL) {
    n <- length(x)
    u <- rep(0, n)
    for (i in 1:n) {
        if (!is.null(checkpoint)) checkpoint()
        u[i] <- theta(x[-i])
    }
    sqrt(((n - 1) / n) * (n - 1) * stats::var(u))
}

sdid_jackknife_se <- function(estimate, weights = attr(estimate, "weights"), checkpoint = NULL) {
    setup <- attr(estimate, "setup")
    opts <- attr(estimate, "opts")
    if (!is.null(weights)) opts$update.omega <- opts$update.lambda <- FALSE
    if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) return(NA)
    theta <- function(ind) {
        weights.jk <- weights
        if (!is.null(weights)) weights.jk$omega <- sdid_sum_normalize(weights$omega[ind[ind <= setup$N0]])
        do.call(sdid_estimate, c(list(Y = setup$Y[ind, ], N0 = sum(ind <= setup$N0), T0 = setup$T0,
                                      X = setup$X[ind, , ], weights = weights.jk), opts))
    }
    sdid_jackknife(1:nrow(setup$Y), theta, checkpoint)
}

sdid_placebo_se <- function(estimate, replications, checkpoint = NULL) {
    setup <- attr(estimate, "setup")
    opts <- attr(estimate, "opts")
    weights <- attr(estimate, "weights")
    N1 <- nrow(setup$Y) - setup$N0
    if (setup$N0 <= N1) stop("must have more controls than treated units to use the placebo se")
    theta <- function(ind) {
        N0 <- length(ind) - N1
        weights.boot <- weights
        weights.boot$omega <- sdid_sum_normalize(weights$omega[ind[1:N0]])
        do.call(sdid_estimate, c(list(Y = setup$Y[ind, ], N0 = N0, T0 = setup$T0,
                                      X = setup$X[ind, , ], weights = weights.boot), opts))
    }
    reps <- vapply(seq_len(replications), function(i) {
        if (!is.null(checkpoint)) checkpoint()
        as.numeric(theta(sample(1:setup$N0)))
    }, numeric(1))
    sqrt((replications - 1) / replications) * stats::sd(reps)
}

sdid_se <- function(estimate, method = c("placebo", "bootstrap", "jackknife"), replications = 200,
                    checkpoint = NULL) {
    method <- match.arg(method)
    switch(method,
           bootstrap = sdid_bootstrap_se(estimate, replications, checkpoint),
           jackknife = sdid_jackknife_se(estimate, checkpoint = checkpoint),
           placebo = sdid_placebo_se(estimate, replications, checkpoint))
}
