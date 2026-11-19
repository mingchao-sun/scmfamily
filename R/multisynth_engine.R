# Staggered-adoption augmented synthetic control ("multisynth").
#
# This file is a port of the multisynth part of the R package 'augsynth'
# (https://github.com/ebenmichael/augsynth, Eli Ben-Michael, MIT licence): multisynth_class.R,
# multi_synth_qp.R, outcome_multi.R, eligible_donors.R, format.R (format_data_stag,
# check_data_stag, extract_covariates) and inference.R (jackknife_se_multi,
# weighted_bootstrap_multi). The algorithms are kept line for line; dplyr / purrr / tidyr
# plumbing is replaced by base R and function names are prefixed `ms_`.
# Supported: balancing weights with partial pooling (nu), unit and time fixed effects or time
# effects only, balanced auxiliary covariates, time cohorts, jackknife and wild-bootstrap
# inference. Not ported: matching on covariates (exact / knn), multiple outcomes, and the
# time-regression outcome model (time_w), none of which the default `multisynth()` call uses.
# Reference: Ben-Michael, Feller and Rothstein (2022), JRSS-B 84(2), 351-381.

# ---- data (format.R) -------------------------------------------------------------------------

# Y: T x N outcomes (colnames = units), times sorted; adopt: named vector of first treated
# period per unit (NA for never treated), names = units.
ms_format <- function(Y, times, adopt) {
    units <- sort(colnames(Y))
    Y <- Y[, units, drop = FALSE]
    adopt <- adopt[units]
    trt_time <- ifelse(is.na(adopt), Inf, adopt)
    t_int <- max(trt_time[is.finite(trt_time)])
    t_vec <- sort(unique(times))
    pre_times <- t_vec[t_vec < t_int]
    post_times <- t_vec[t_vec >= t_int]
    Xy <- t(Y)
    colnames(Xy) <- as.character(t_vec)
    X <- Xy[, as.character(pre_times), drop = FALSE]
    y <- Xy[, as.character(post_times), drop = FALSE]
    trt <- vapply(trt_time, function(x) if (is.finite(x)) which(t_vec == x) - 1 else Inf, numeric(1))
    names(trt) <- NULL
    fin <- is.finite(trt_time)
    mask <- t(vapply(which(fin), function(i) as.numeric(t_vec[t_vec < t_int] < trt_time[i]), numeric(length(pre_times))))
    if (sum(fin) == 1) mask <- matrix(mask, nrow = 1)
    list(X = X, trt = trt, y = y, mask = mask, time = t_vec, units = units)
}

ms_check_data <- function(wide, fixedeff) {
    less_5 <- wide$units[wide$trt < 5]
    warn_msg <- ""
    if (length(less_5) != 0)
        warn_msg <- paste0("The following units have less than 5 pre-treatment outcomes: (",
                           paste(less_5, collapse = ","), "). Be cautious!")
    always_trt <- wide$units[wide$trt == 0]
    n1 <- wide$units[wide$trt == 1]
    err_msg <- ""
    if (length(always_trt) != 0)
        err_msg <- paste0(err_msg, "The following units are always treated and should be removed: (",
                          paste(always_trt, collapse = ","), ")\n")
    if (length(n1) != 0 & fixedeff) {
        if (nchar(err_msg) > 0) err_msg <- paste0(err_msg, "  Also: ")
        err_msg <- paste0(err_msg, "You are including a fixed effect, but the following units only have ",
                          "one pre-treatment outcome: (", paste(n1, collapse = ","),
                          "). Either remove these units or switch off unit fixed effects.\n")
    }
    if (max(wide$trt) < ncol(wide$X) + ncol(wide$y)) {
        if (nchar(err_msg) > 0) err_msg <- paste0(err_msg, "  Also: ")
        err_msg <- paste0(err_msg, "All units are eventually treated. The last treatment time is ",
                          wide$time[max(wide$trt)], ". Remove all periods after this time, or ",
                          "include never-treated units.\n")
    }
    if (nchar(warn_msg) > 0) warning(warn_msg)
    if (nchar(err_msg) > 0) stop(err_msg)
    invisible(NULL)
}

# covariates: mean over the periods before the earliest adoption (omitting NAs)
ms_covariates <- function(cov, times, adopt, units) {
    t_first <- min(adopt, na.rm = TRUE)
    pre <- which(times < t_first)
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

# ---- eligible donors (eligible_donors.R) -----------------------------------------------------

ms_get_eligible_donors <- function(trt, time_cohort, n_leads) {
    grps <- if (time_cohort) unique(trt[is.finite(trt)]) else trt[is.finite(trt)]
    J <- length(grps)
    lapply(1:J, function(j) trt > n_leads + grps[j])
}

ms_get_nona_donors <- function(X, y, trt, n_lags, n_leads, time_cohort) {
    n <- length(trt)
    fulldat <- cbind(X, y)
    is_na <- is.na(fulldat[is.finite(trt), , drop = FALSE])
    if (time_cohort) {
        grps <- unique(trt[is.finite(trt)])
        finite_trt <- trt[is.finite(trt)]
        is_na <- t(sapply(grps, function(tj) apply(is_na[finite_trt == tj, , drop = FALSE], 2, all)))
    } else {
        grps <- trt[is.finite(trt)]
    }
    not_na <- !is.na(fulldat)
    J <- length(grps)
    lapply(1:J, function(j) {
        idxs <- max(grps[j] - n_lags + 1, 1):min(grps[j] + n_leads, ncol(fulldat))
        isna_j <- is_na[j, idxs]
        apply(not_na[, idxs, drop = FALSE][, !isna_j, drop = FALSE], 1, all)
    })
}

ms_get_donors <- function(X, y, trt, time_cohort, n_lags, n_leads) {
    donors <- ms_get_eligible_donors(trt, time_cohort, n_leads)
    nona_donors <- ms_get_nona_donors(X, y, trt, n_lags, n_leads, time_cohort)
    lapply(1:length(donors), function(j) donors[[j]] & nona_donors[[j]])
}

# ---- outcome model: fixed effects (outcome_multi.R) --------------------------------------------

ms_fit_feff <- function(X, trt, mask, force, time_cohort) {
    ttot <- ncol(X)
    n <- nrow(X)
    grps <- unique(trt[is.finite(trt)])
    J <- length(grps)
    which_t <- (1:n)[is.finite(trt)]
    if (force %in% c(2, 3)) {
        time_eff <- matrix(colMeans(X[!is.finite(trt), , drop = FALSE], na.rm = TRUE),
                           nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
    } else {
        time_eff <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    }
    residuals <- X - time_eff
    y0hat <- time_eff
    if (force %in% c(1, 3)) {
        unit_eff <- lapply(grps, function(tj) matrix(rowMeans(residuals[, 1:tj, drop = FALSE], na.rm = TRUE),
                                                     nrow = nrow(X), ncol = ncol(X)))
        residuals <- lapply(1:J, function(j) residuals - unit_eff[[j]])
        y0hat <- unit_eff
    }
    if (force == 3) y0hat <- lapply(unit_eff, function(ufj) time_eff + ufj)
    if (force %in% c(1, 3) & !time_cohort) {
        names(residuals) <- as.character(grps)
        residuals <- residuals[as.character(trt[is.finite(trt)])]
        names(y0hat) <- as.character(grps)
        y0hat <- y0hat[as.character(trt[is.finite(trt)])]
    }
    list(y0hat = y0hat, residuals = residuals)
}

# ---- the QP for the weights (multi_synth_qp.R) ----------------------------------------------------

ms_qp <- function(X, trt, mask, Z = NULL, n_leads = NULL, n_lags = NULL,
                  relative = TRUE, nu = 0, lambda = 0, V = NULL, time_cohort = FALSE,
                  donors = NULL, norm_pool = NULL, norm_sep = NULL,
                  verbose = FALSE, eps_rel = 1e-4, eps_abs = 1e-4) {
    if (!requireNamespace("osqp", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE))
        stop("The R packages 'osqp' and 'Matrix' are required. Install them with install.packages(c('osqp', 'Matrix')).")
    if (!is.null(Z)) if (ncol(Z) == 0) Z <- NULL
    n <- if (typeof(X) == "list") dim(X[[1]])[1] else dim(X)[1]
    d <- if (typeof(X) == "list") dim(X[[1]])[2] else dim(X)[2]
    if (is.null(n_leads)) n_leads <- d + 1 else if (n_leads > d) n_leads <- d + 1
    if (is.null(n_lags)) n_lags <- d else if (n_lags > d) n_lags <- d
    V <- asyn_make_V(n_lags, V)
    if (time_cohort) {
        grps <- unique(trt[is.finite(trt)])
        which_t <- lapply(grps, function(tj) (1:n)[trt == tj])
        mask <- unique(mask)
    } else {
        grps <- trt[is.finite(trt)]
        which_t <- (1:n)[is.finite(trt)]
    }
    J <- length(grps)
    if (is.null(norm_sep)) norm_sep <- 1
    if (is.null(norm_pool)) norm_pool <- 1
    n1 <- sapply(1:J, function(j) length(which_t[[j]]))
    if (is.null(donors)) donors <- ms_get_eligible_donors(trt, time_cohort, n_leads)

    if (typeof(X) == "list") {
        x_t <- lapply(1:J, function(j) colSums(X[[j]][which_t[[j]], mask[j, ] == 1, drop = FALSE]))
        Xc <- lapply(1:nrow(mask), function(j) X[[j]][donors[[j]], mask[j, ] == 1, drop = FALSE])
        sdx <- stats::sd(X[[1]][is.finite(trt)])
    } else {
        x_t <- lapply(1:J, function(j) colSums(X[which_t[[j]], mask[j, ] == 1, drop = FALSE]))
        Xc <- lapply(1:nrow(mask), function(j) X[donors[[j]], mask[j, ] == 1, drop = FALSE])
        sdx <- stats::sd(X[is.finite(trt)])
    }
    if (!is.null(Z)) {
        Z_scale <- sdx * apply(Z, 2, function(z) (z - mean(z[!is.finite(trt)])) / stats::sd(z[!is.finite(trt)]))
        z_t <- lapply(1:J, function(j) colSums(Z_scale[which_t[[j]], , drop = FALSE]))
        Zc <- lapply(1:J, function(j) Z_scale[donors[[j]], , drop = FALSE])
    } else {
        z_t <- lapply(1:J, function(j) c(0))
        Zc <- lapply(1:J, function(j) Matrix::Matrix(0, nrow = sum(donors[[j]]), ncol = 1))
    }
    dz <- ncol(Zc[[1]])

    x_t <- lapply(x_t, function(xtk) { xtk[is.na(xtk)] <- 0; xtk })
    Xc <- lapply(Xc, function(xck) { xck[is.na(xck)] <- 0; xck })
    n0s <- sapply(Xc, nrow)
    if (any(n0s == 0)) stop("Some treated units have no possible donor units!")
    n0 <- sum(n0s)

    const_mats <- ms_make_constraint_mats(trt, grps, n_leads, n_lags, Xc, Zc, d, n1)
    Amat <- const_mats$Amat
    lvec <- const_mats$lvec
    uvec <- const_mats$uvec
    qvec <- ms_make_qvec(Xc, x_t, z_t, nu, n_lags, d, V, norm_pool, norm_sep)
    Pmat <- ms_make_Pmat(Xc, x_t, dz, nu, n_lags, lambda, d, V, norm_pool, norm_sep)

    settings <- do.call(osqp::osqpSettings, c(list(verbose = verbose, eps_rel = eps_rel, eps_abs = eps_abs)))
    out <- osqp::solve_osqp(Pmat, qvec, Amat, lvec, uvec, pars = settings)

    total_ctrls <- n0 * J
    weights <- matrix(out$x[1:total_ctrls], nrow = n0)
    nj0 <- as.numeric(lapply(Xc, nrow))
    nj0cumsum <- c(0, cumsum(nj0))
    imbalance <- vapply(1:J, function(j) {
        dj <- length(x_t[[j]])
        ndim <- min(dj, n_lags)
        c(numeric(d - ndim),
          x_t[[j]][(dj - ndim + 1):dj] -
              t(Xc[[j]][, (dj - ndim + 1):dj, drop = FALSE]) %*% out$x[(nj0cumsum[j] + 1):nj0cumsum[j + 1]])
    }, numeric(d))
    avg_imbal <- rowMeans(t(t(imbalance)))

    Vsq <- t(V) %*% V
    global_l2 <- c(sqrt(t(avg_imbal[(d - n_lags + 1):d]) %*% Vsq %*% avg_imbal[(d - n_lags + 1):d])) / sqrt(d)
    avg_l2 <- mean(apply(imbalance, 2, function(x) c(sqrt(t(x[(d - n_lags + 1):d]) %*% Vsq %*% x[(d - n_lags + 1):d]))))
    ind_l2 <- sqrt(mean(apply(imbalance, 2, function(x)
        c(x[(d - n_lags + 1):d] %*% Vsq %*% x[(d - n_lags + 1):d]) / sum(x[(d - n_lags + 1):d] != 0))))
    weights <- vapply(1:J, function(j) {
        weightj <- numeric(n)
        weightj[donors[[j]]] <- out$x[(nj0cumsum[j] + 1):nj0cumsum[j + 1]]
        weightj
    }, numeric(n))
    weights <- t(t(weights) / n1)
    weights <- pmax(weights, 0)
    col_sums <- colSums(weights)
    weights <- sweep(weights, 2, ifelse(col_sums > 0, col_sums, 1), "/")

    output <- list(weights = weights, imbalance = cbind(avg_imbal, imbalance),
                   global_l2 = global_l2, ind_l2 = ind_l2, avg_l2 = avg_l2, V = V)
    if (!is.null(Z)) {
        z_t <- sapply(1:J, function(j) colMeans(Z[which_t[[j]], , drop = FALSE]))
        imbal_z <- z_t - t(Z) %*% weights
        avg_imbal_z <- rowSums(t(t(imbal_z) * n1)) / sum(n1)
        global_l2_z <- sqrt(sum(avg_imbal_z^2))
        ind_l2_z <- sum(apply(imbal_z, 2, function(x) sqrt(sum(x^2))))
        imbal_z <- cbind(avg_imbal_z, imbal_z)
        rownames(imbal_z) <- colnames(Z)
        output$imbalance_aux <- imbal_z
        output$global_l2_aux <- global_l2_z
        output$ind_l2_aux <- ind_l2_z
    }
    output
}

ms_make_constraint_mats <- function(trt, grps, n_leads, n_lags, Xc, Zc, d, n1) {
    J <- length(grps)
    A1 <- Matrix::bdiag(lapply(1:J, function(j) rep(1, nrow(Xc[[j]]))))
    Amat <- as.matrix(Matrix::t(A1))
    Amat <- Matrix::rbind2(Matrix::t(A1), Matrix::Diagonal(nrow(A1)))
    dz <- ncol(Zc[[1]])
    A_trans1 <- do.call(Matrix::bdiag, lapply(1:J, function(j) {
        dj <- ncol(Xc[[j]])
        ndim <- min(dj, n_lags)
        max_dim <- min(d, n_lags)
        mat <- Xc[[j]][, (dj - ndim + 1):dj, drop = FALSE]
        n0 <- nrow(mat)
        zero_mat <- Matrix::Matrix(0, n0, max_dim - ndim)
        Matrix::t(cbind(zero_mat, mat))
    }))
    sum_tj <- min(d, n_lags) * J
    A_trans2 <- -Matrix::Diagonal(sum_tj)
    A_trans <- Matrix::cbind2(Matrix::cbind2(A_trans1, A_trans2),
                              Matrix::Matrix(0, nrow = nrow(A_trans1), ncol = dz * J))
    A_transz <- Matrix::t(Matrix::bdiag(Zc))
    A_transz <- Matrix::cbind2(Matrix::cbind2(A_transz, Matrix::Matrix(0, nrow = nrow(A_transz), ncol = sum_tj)),
                               -Matrix::Diagonal(dz * J))
    Amat <- Matrix::cbind2(Amat, Matrix::Matrix(0, nrow = nrow(Amat), ncol = sum_tj + dz * J))
    Amat <- Matrix::rbind2(Matrix::rbind2(Amat, A_trans), A_transz)
    lvec <- c(n1, numeric(nrow(A1)), numeric(sum_tj), numeric(dz * J))
    uvec <- c(n1, rep(Inf, nrow(A1)), numeric(sum_tj), numeric(dz * J))
    list(Amat = Amat, lvec = lvec, uvec = uvec)
}

ms_make_qvec <- function(Xc, x_t, z_t, nu, n_lags, d, V, norm_pool, norm_sep) {
    J <- length(x_t)
    Vsq <- t(V) %*% V
    qvec <- lapply(1:J, function(j) {
        dj <- length(x_t[[j]])
        ndim <- min(dj, n_lags)
        max_dim <- min(d, n_lags)
        vec <- x_t[[j]][(dj - ndim + 1):dj] / ndim
        Vsq %*% c(numeric(max_dim - ndim), vec)
    })
    avg_target_vec <- Reduce(`+`, lapply(x_t, function(xtk) {
        dk <- length(xtk)
        ndim <- min(dk, n_lags)
        max_dim <- min(d, n_lags)
        c(numeric(max_dim - ndim), xtk[(dk - ndim + 1):dk])
    })) %*% Vsq
    qvec_avg <- rep(avg_target_vec, J)
    qvec <- -(nu * qvec_avg / (norm_pool * n_lags * J^2) + (1 - nu) * Reduce(c, qvec) / (norm_sep * J))
    qvec_avg_z <- Reduce(`+`, z_t)
    qvec_avg_z <- rep(qvec_avg_z, J)
    qvec_z <- -(nu * qvec_avg_z / (norm_pool * J^2) + (1 - nu) * Reduce(c, z_t) / (norm_sep * J)) / length(z_t[[1]])
    total_ctrls <- sum(unlist(lapply(Xc, nrow)))
    c(numeric(total_ctrls), qvec, qvec_z)
}

ms_make_Pmat <- function(Xc, x_t, dz, nu, n_lags, lambda, d, V, norm_pool, norm_sep) {
    J <- length(x_t)
    Vsq <- t(V) %*% V
    ndims <- vapply(1:J, function(j) min(length(x_t[[j]]), n_lags), numeric(1))
    max_dim <- min(d, n_lags)
    total_dim <- max_dim * J
    V1 <- Matrix::bdiag(lapply(ndims, function(ndim) Vsq / ndim))
    tile_sparse <- function(j) kronecker(Matrix::Matrix(1, nrow = j, ncol = j), Vsq)
    tile_sparse_cov <- function(d, j) kronecker(Matrix::Matrix(1, nrow = j, ncol = j), Matrix::Diagonal(d))
    V2 <- tile_sparse(J) / n_lags
    Pmat <- nu * V2 / (norm_pool * J^2) + (1 - nu) * V1 / (norm_sep * J)
    V1_z <- Matrix::Diagonal(dz * J, 1 / dz)
    V2_z <- tile_sparse_cov(dz, J) / dz
    Pmat_z <- nu * V2_z / (norm_pool * J^2) + (1 - nu) * V1_z / (norm_sep * J)
    total_ctrls <- sum(unlist(lapply(Xc, nrow)))
    Pmat <- Matrix::bdiag(Matrix::Matrix(0, nrow = total_ctrls, ncol = total_ctrls), Pmat, Pmat_z)
    I0 <- Matrix::bdiag(Matrix::Diagonal(total_ctrls),
                        Matrix::Matrix(0, nrow = total_dim + dz * J, ncol = total_dim + dz * J))
    Pmat + lambda * I0
}

# ---- fitting (multisynth_class.R) --------------------------------------------------------------------

ms_fit_formatted <- function(wide, n_leads, n_lags, nu, lambda, V, force, scm, time_cohort,
                             eps_abs, eps_rel, verbose) {
    relative <- TRUE
    grps <- if (time_cohort) unique(wide$trt[is.finite(wide$trt)]) else wide$trt[is.finite(wide$trt)]
    J <- length(grps)
    if (force == 0) {
        pure_ctrl <- cbind(wide$X, wide$y)[!is.finite(wide$trt), , drop = FALSE]
        y0hat <- matrix(colMeans(pure_ctrl, na.rm = TRUE), nrow = nrow(wide$X), ncol = ncol(pure_ctrl), byrow = TRUE)
        residuals <- cbind(wide$X, wide$y) - y0hat
        params <- NULL
    } else {
        fullmask <- cbind(wide$mask, matrix(0, nrow = nrow(wide$mask), ncol = ncol(wide$y)))
        out <- ms_fit_feff(cbind(wide$X, wide$y), wide$trt, fullmask, force, time_cohort)
        y0hat <- out$y0hat
        residuals <- out$residuals
        params <- NULL
    }
    if (typeof(residuals) == "list") bal_mat <- lapply(residuals, function(x) x[, 1:ncol(wide$X)])
    else bal_mat <- residuals[, 1:ncol(wide$X)]

    Zw <- wide$Z
    if (scm) {
        donors <- ms_get_donors(wide$X, wide$y, wide$trt, time_cohort, n_lags, n_leads)
        sep_fit <- ms_qp(X = bal_mat, trt = wide$trt, mask = wide$mask, Z = Zw,
                         n_leads = n_leads, n_lags = n_lags, relative = relative, nu = 0, lambda = lambda,
                         V = V, time_cohort = time_cohort, donors = donors,
                         eps_rel = eps_rel, eps_abs = eps_abs, verbose = verbose)
        if (is.null(nu)) {
            glbl <- sep_fit$global_l2 * sqrt(nrow(sep_fit$imbalance))
            ind <- sep_fit$avg_l2
            nu <- glbl / ind
        }
        msynth <- ms_qp(X = bal_mat, trt = wide$trt, mask = wide$mask, Z = Zw,
                        n_leads = n_leads, n_lags = n_lags, relative = relative, nu = nu, lambda = lambda,
                        V = V, time_cohort = time_cohort, donors = donors,
                        norm_pool = sep_fit$global_l2^2, norm_sep = sep_fit$ind_l2^2,
                        eps_rel = eps_rel, eps_abs = eps_abs, verbose = verbose)
    } else {
        donors <- NULL
        msynth <- list(weights = matrix(0, nrow = nrow(wide$X), ncol = J), imbalance = NA,
                       global_l2 = NA, ind_l2 = NA)
    }
    msynth$data <- wide
    msynth$relative <- relative
    msynth$n_leads <- n_leads
    msynth$n_lags <- n_lags
    msynth$nu <- nu
    msynth$lambda <- lambda
    msynth$scm <- scm
    msynth$time_cohort <- time_cohort
    msynth$grps <- grps
    msynth$y0hat <- y0hat
    msynth$residuals <- residuals
    msynth$force <- force
    msynth$params <- params
    msynth$extra_pars <- list(eps_abs = eps_abs, eps_rel = eps_rel, verbose = verbose)
    msynth$donors <- donors
    class(msynth) <- "ms_fit"
    msynth
}

# Y: T x N panel (colnames = units); adopt: named vector of first treated periods (NA = never)
ms_fit <- function(Y, times, adopt, Z = NULL, n_leads = NULL, n_lags = NULL, nu = NULL, lambda = 0,
                   V = NULL, fixedeff = TRUE, scm = TRUE, time_cohort = FALSE,
                   eps_abs = 1e-4, eps_rel = 1e-4, verbose = FALSE) {
    wide <- ms_format(Y, times, adopt)
    ms_check_data(wide, fixedeff)
    force <- if (fixedeff) 3 else 2
    if (!is.null(Z)) Z <- Z[wide$units, , drop = FALSE]
    wide$Z <- Z
    if (is.null(n_leads)) {
        n_leads <- ncol(wide$y)
    } else if (n_leads > max(apply(1 - wide$mask, 1, sum, na.rm = TRUE)) + ncol(wide$y)) {
        n_leads <- max(apply(1 - wide$mask, 1, sum, na.rm = TRUE)) + ncol(wide$y)
    }
    if (is.null(n_lags)) n_lags <- ncol(wide$X) else if (n_lags > ncol(wide$X)) n_lags <- ncol(wide$X)

    msynth <- ms_fit_formatted(wide, n_leads, n_lags, nu, lambda, V, force, scm, time_cohort,
                               eps_abs, eps_rel, verbose)
    rownames(msynth$weights) <- wide$units
    if (scm) {
        unif <- ms_qp(X = wide$X, trt = wide$trt, mask = wide$mask, Z = Z,
                      n_leads = n_leads, n_lags = n_lags, relative = TRUE, nu = 0, lambda = 1e10, V = V,
                      time_cohort = time_cohort, donors = msynth$donors, eps_rel = eps_rel,
                      eps_abs = eps_abs, verbose = verbose)
        msynth$scaled_global_l2 <- msynth$global_l2 / unif$global_l2
        msynth$scaled_ind_l2 <- msynth$ind_l2 / unif$ind_l2
    }
    msynth$trt_unit <- msynth$data$units[msynth$data$trt < Inf]
    msynth
}

# ---- prediction (relative event time only) -------------------------------------------------------------------

ms_predict <- function(object, att = FALSE, att_weight = NULL, bs_weight = NULL) {
    multisynth <- object
    time_cohort <- multisynth$time_cohort
    n_leads <- multisynth$n_leads
    d <- ncol(multisynth$data$X)
    n <- nrow(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)
    grps <- multisynth$grps
    J <- length(grps)
    if (is.null(bs_weight)) bs_weight <- rep(1, n)
    if (time_cohort) {
        which_t <- lapply(grps, function(tj) (1:n)[multisynth$data$trt == tj])
        mask <- unique(multisynth$data$mask)
    } else {
        which_t <- (1:n)[is.finite(multisynth$data$trt)]
        mask <- multisynth$data$mask
    }
    n1 <- sapply(1:J, function(j) length(which_t[[j]]))
    fullmask <- cbind(mask, matrix(0, nrow = J, ncol = (ttot - d)))

    mu1hat <- vapply(1:J, function(j) colMeans((bs_weight * fulldat)[which_t[[j]], , drop = FALSE]), numeric(ttot))
    if (typeof(multisynth$y0hat) == "list") {
        mu0hat <- vapply(1:J, function(j) {
            y0hat <- colMeans((bs_weight * multisynth$y0hat[[j]])[which_t[[j]], , drop = FALSE])
            weightsj <- multisynth$weights[, j] * bs_weight
            resj <- multisynth$residuals[[j]][weightsj != 0, , drop = FALSE]
            y0hat + t(resj) %*% weightsj[weightsj != 0]
        }, numeric(ttot))
    } else {
        mu0hat <- vapply(1:J, function(j) {
            y0hat <- colMeans((bs_weight * multisynth$y0hat)[which_t[[j]], , drop = FALSE])
            weightsj <- multisynth$weights[, j] * bs_weight
            resj <- multisynth$residuals[weightsj != 0, , drop = FALSE]
            y0hat + t(resj) %*% weightsj[weightsj != 0]
        }, numeric(ttot))
    }
    tauhat <- mu1hat - mu0hat
    if (is.null(att_weight)) att_weight <- rep(1, J)

    total_len <- min(d + n_leads, ttot + d - min(grps))
    mu0hat <- vapply(1:J, function(j) {
        vec <- c(rep(NA, d - grps[j]), mu0hat[1:grps[j], j],
                 mu0hat[(grps[j] + 1):(min(grps[j] + n_leads, ttot)), j])
        c(vec, rep(NA, total_len - length(vec)), mean(mu0hat[(grps[j] + 1):(min(grps[j] + n_leads, ttot)), j]))
    }, numeric(total_len + 1))
    tauhat <- vapply(1:J, function(j) {
        vec <- c(rep(NA, d - grps[j]), tauhat[1:grps[j], j],
                 tauhat[(grps[j] + 1):(min(grps[j] + n_leads, ttot)), j])
        c(vec, rep(NA, total_len - length(vec)), mean(tauhat[(grps[j] + 1):(min(grps[j] + n_leads, ttot)), j]))
    }, numeric(total_len + 1))
    if (is.null(dim(att_weight))) {
        if (J == 1) att_weight <- matrix(replicate(total_len + 1, att_weight), ncol = 1)
        else att_weight <- t(replicate(total_len + 1, att_weight))
    }
    att_weight_new <- vapply(1:J, function(j) {
        vec <- c(rep(NA, d - grps[j]), att_weight[1:grps[j], j],
                 att_weight[(grps[j] + 1):(min(grps[j] + n_leads, ttot)), j])
        c(vec, rep(NA, total_len - length(vec)), mean(att_weight[(grps[j] + 1):(min(grps[j] + n_leads, ttot)), j]))
    }, numeric(total_len + 1))

    avg <- sapply(1:nrow(mu0hat), function(k)
        sum(n1 * mu0hat[k, ] * att_weight_new[k, ], na.rm = TRUE) /
            sum(n1 * (!is.na(mu0hat[k, ])) * att_weight_new[k, ], na.rm = TRUE))
    mu0hat <- cbind(avg, mu0hat)
    avg <- sapply(1:nrow(mu0hat), function(k)
        sum(n1 * tauhat[k, ] * att_weight_new[k, ], na.rm = TRUE) /
            sum(n1 * (!is.na(tauhat[k, ])) * att_weight_new[k, ], na.rm = TRUE))
    tauhat <- cbind(avg, tauhat)
    if (att) tauhat else mu0hat
}

# ---- inference (inference.R) ------------------------------------------------------------------------------------------

ms_drop_unit_i <- function(msyn, i) {
    n <- nrow(msyn$data$X)
    time_cohort <- msyn$time_cohort
    which_t <- (1:n)[is.finite(msyn$data$trt)]
    not_miss_j <- which_t %in% setdiff(which_t, i)
    drop_i <- msyn$data
    drop_i$X <- msyn$data$X[-i, , drop = FALSE]
    drop_i$y <- msyn$data$y[-i, , drop = FALSE]
    drop_i$trt <- msyn$data$trt[-i]
    drop_i$mask <- msyn$data$mask[not_miss_j, , drop = FALSE]
    if (!is.null(msyn$data$Z)) drop_i$Z <- msyn$data$Z[-i, , drop = FALSE] else drop_i$Z <- NULL
    args_list <- list(wide = drop_i, n_leads = msyn$n_leads, n_lags = msyn$n_lags,
                      nu = msyn$nu, lambda = msyn$lambda, V = msyn$V, force = msyn$force,
                      scm = msyn$scm, time_cohort = msyn$time_cohort)
    msyn_i <- do.call(ms_fit_formatted, c(args_list, msyn$extra_pars))
    if (time_cohort) dropped <- which(!msyn$grps %in% msyn_i$grps) else dropped <- which(!not_miss_j)
    list(msyn_i, dropped)
}

ms_jackknife_se <- function(multisynth, alpha = 0.05, att_weight = NULL, checkpoint = NULL) {
    n <- nrow(multisynth$data$X)
    att <- ms_predict(multisynth, att = TRUE, att_weight = att_weight)
    outddim <- nrow(att)
    J <- length(multisynth$grps)
    jack_est <- vapply(1:n, function(i) {
        if (!is.null(checkpoint)) checkpoint()
        msyn_i <- ms_drop_unit_i(multisynth, i)
        pred <- ms_predict(msyn_i[[1]], att = TRUE, att_weight = att_weight)
        if (nrow(pred) < outddim) {
            pred <- rbind(pred[1:(nrow(pred) - 1), ],
                          matrix(NA, nrow = outddim - nrow(pred), ncol = ncol(pred)),
                          pred[nrow(pred), ])
        }
        if (length(msyn_i[[2]]) != 0) {
            out <- matrix(NA, nrow = nrow(pred), ncol = (J + 1))
            out[, -(msyn_i[[2]] + 1)] <- pred
        } else {
            out <- pred
        }
        out
    }, matrix(0, nrow = outddim, ncol = (J + 1)))
    se2 <- apply(jack_est, c(1, 2), function(x) (n - 1) / n * sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
    list(att = att, se = sqrt(se2),
         lower_bound = att - stats::qnorm(1 - alpha / 2) * sqrt(se2),
         upper_bound = att + stats::qnorm(1 - alpha / 2) * sqrt(se2))
}

ms_rwild <- function(n) {
    sample(c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2), n, replace = TRUE,
           prob = c((sqrt(5) + 1) / (2 * sqrt(5)), (sqrt(5) - 1) / (2 * sqrt(5))))
}

ms_bootstrap <- function(multisynth, rweight = ms_rwild, n_boot = 1000, alpha = 0.05, att_weight = NULL,
                         checkpoint = NULL) {
    n <- nrow(multisynth$data$X)
    att <- ms_predict(multisynth, att = TRUE, att_weight = att_weight)
    outddim <- nrow(att)
    n1 <- sum(is.finite(multisynth$data$trt))
    J <- length(multisynth$grps)
    bs_est <- vapply(1:n_boot, function(i) {
        if (!is.null(checkpoint)) checkpoint()
        Z <- rweight(n)
        ms_predict(multisynth, att = TRUE, att_weight = att_weight, bs_weight = Z) - sum(Z) / n1 * att
    }, matrix(0, nrow = outddim, ncol = (J + 1)))
    se2 <- apply(bs_est, c(1, 2), function(x) mean((x - mean(x))^2, na.rm = TRUE))
    bias <- apply(bs_est, c(1, 2), function(x) mean(x, na.rm = TRUE))
    upper_bound <- att - apply(bs_est, c(1, 2), function(x) stats::quantile(x, alpha / 2, na.rm = TRUE))
    lower_bound <- att - apply(bs_est, c(1, 2), function(x) stats::quantile(x, 1 - alpha / 2, na.rm = TRUE))
    list(att = att, bias = bias, se = sqrt(se2), upper_bound = upper_bound, lower_bound = lower_bound)
}

# summary: event-time table with Average and per-group columns (as summary.multisynth)
ms_summary <- function(object, inf_type = c("bootstrap", "jackknife", "none"), att_weight = NULL, ...) {
    inf_type <- match.arg(inf_type)
    multisynth <- object
    n_leads <- multisynth$n_leads
    d <- ncol(multisynth$data$X)
    n <- nrow(multisynth$data$X)
    ttot <- d + ncol(multisynth$data$y)
    trt <- multisynth$data$trt
    time_cohort <- multisynth$time_cohort
    if (time_cohort) {
        grps <- unique(trt[is.finite(trt)])
    } else {
        grps <- trt[is.finite(trt)]
    }
    which_t <- (1:n)[is.finite(trt)]
    times <- multisynth$data$time
    if (inf_type == "jackknife") {
        attse <- ms_jackknife_se(multisynth, att_weight = att_weight, ...)
    } else if (inf_type == "bootstrap") {
        attse <- ms_bootstrap(multisynth, att_weight = att_weight, ...)
    } else {
        att <- ms_predict(multisynth, att = TRUE, att_weight = att_weight)
        attse <- list(att = att, se = matrix(NA, nrow(att), ncol(att)),
                      upper_bound = matrix(NA, nrow(att), ncol(att)),
                      lower_bound = matrix(NA, nrow(att), ncol(att)))
    }
    event <- c(-(d - 1):min(n_leads, ttot - min(grps)), NA) - 1
    col_names <- if (time_cohort) c("Average", as.character(times[grps + 1]))
                 else c("Average", as.character(multisynth$data$units[which_t]))
    colnames(attse$att) <- colnames(attse$se) <- colnames(attse$lower_bound) <- colnames(attse$upper_bound) <- col_names
    list(event = event, att = attse$att, se = attse$se, lower = attse$lower_bound, upper = attse$upper_bound,
         groups = col_names, grps = grps, n_leads = n_leads, n_lags = multisynth$n_lags,
         global_l2 = multisynth$global_l2, scaled_global_l2 = multisynth$scaled_global_l2,
         ind_l2 = multisynth$ind_l2, scaled_ind_l2 = multisynth$scaled_ind_l2, nu = multisynth$nu)
}
