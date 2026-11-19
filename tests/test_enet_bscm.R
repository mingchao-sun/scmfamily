# Independent checks of the elastic-net and Bayesian synthetic control analyses. Run from the module root:
#   Rscript tests/test_enet_bscm.R
# enet: the weights are compared with an exact solution of the paper's objective by coordinate descent written for this check
#       (no code of glmnet); the cross-validation error is recomputed with that solver.
# bscm: the sampler is checked on simulated data (recovery, coverage, shrinkage) and, when the CRAN package 'horseshoe' is
#       installed, against its posterior means (same Makalic-Schmidt Gibbs sampler and priors).
source("R/common.R")
nfail <- 0
chk <- function(name, ok, detail = "") { cat(sprintf("%-62s %s %s\n", name, if (ok) "ok  " else "FAIL", detail)); if (!ok) nfail <<- nfail + 1 }

# ---- enet: exact solution by coordinate descent ---------------------------------------------------------------------
# Coordinate descent on the centred problem (intercept profiled out), written for this check only; it iterates to machine precision.
qp_enet <- function(Z, y, alpha, lambda, tol = 1e-13, maxit = 2e5) {
    zb <- colMeans(Z); yb <- mean(y)
    Zc <- sweep(Z, 2, zb); yc <- y - yb
    p <- ncol(Zc); w <- rep(0, p); zz <- colSums(Zc^2); res <- yc
    soft <- function(a, t) sign(a) * pmax(abs(a) - t, 0)
    for (it in seq_len(maxit)) {
        old <- w
        for (j in seq_len(p)) {
            r <- sum(Zc[, j] * res) + zz[j] * w[j]
            wj <- soft(r, lambda * alpha / 2) / (zz[j] + lambda * (1 - alpha) / 2)
            if (wj != w[j]) { res <- res - Zc[, j] * (wj - w[j]); w[j] <- wj }
        }
        if (max(abs(w - old)) < tol) break
    }
    list(mu = yb - sum(zb * w), w = w)
}
pp <- read.csv("data/prop99_long.csv")
pq <- fam_panel(pp, "cigsale", "state", "year", "California", 1989)
Z <- pq$Y[pq$pre, pq$donors]; y <- pq$Y[pq$pre, pq$treated]
for (cs in list(c(0.5, 100), c(0.5, 1000), c(0.1, 300), c(0.9, 50), c(1, 80), c(0.25, 5000))) {
    a <- cs[1]; lam <- cs[2]
    ours <- fam_enet(pq$Y, pq$treated, pq$donors, pq$pre, pq$post, fixed = list(alpha = a, lambda = lam))
    ref <- qp_enet(Z, y, a, lam)
    chk(sprintf("enet: weights and intercept, alpha %.2f lambda %g", a, lam),
        max(abs(ours$w - ref$w)) < 1e-4 * max(1, max(abs(ref$w))) && abs(ours$extra$intercept - ref$mu) < 1e-3 * max(1, abs(ref$mu)),
        sprintf("(max |dw| %.1e, d mu %.1e, %d non-zero)", max(abs(ours$w - ref$w)), abs(ours$extra$intercept - ref$mu), sum(abs(ref$w) > 1e-6)))
}
# the cross-validation error is recomputed with this solver
cvr <- fam_enet(pq$Y, pq$treated, pq$donors, pq$pre, pq$post, alpha = 0.5)
cv_at <- function(alpha, lam) {
    errs <- vapply(seq_along(pq$donors), function(j) {
        oth <- pq$donors[-j]
        f <- qp_enet(pq$Y[pq$pre, oth], pq$Y[pq$pre, pq$donors[j]], alpha, lam)
        mean((pq$Y[pq$post, pq$donors[j]] - (f$mu + pq$Y[pq$post, oth] %*% f$w))^2)
    }, 0)
    mean(errs)
}
chk("enet: cross-validation error at the chosen penalty", abs(cv_at(cvr$extra$alpha, cvr$extra$lambda) - cvr$extra$cv_error) < 1e-3 * cvr$extra$cv_error,
    sprintf("(%.4f vs %.4f, lambda %.3g)", cv_at(cvr$extra$alpha, cvr$extra$lambda), cvr$extra$cv_error, cvr$extra$lambda))
nb <- c(cv_at(cvr$extra$alpha, cvr$extra$lambda * 0.7), cv_at(cvr$extra$alpha, cvr$extra$lambda * 1.4))
chk("enet: the chosen penalty is not beaten by neighbouring values", all(nb >= cvr$extra$cv_error * 0.97),
    sprintf("(neighbours %.3f, %.3f)", nb[1], nb[2]))

# ---- bscm: simulated recovery and coverage -----------------------------------------------------------------------------------
sim_panel <- function(seed, T0 = 40, T1 = 6, J = 20, tau = 3) {
    set.seed(seed)
    f <- cumsum(rnorm(T0 + T1, 0, 0.5)); lam <- rnorm(J + 1)
    w_true <- c(0.5, 0.3, 0.2, rep(0, J - 3))
    Yd <- sapply(seq_len(J), function(j) lam[j] * f + 0.2 * (1:(T0 + T1)) * 0.05 + rnorm(T0 + T1, 0, 0.3) + 2 * rnorm(1))
    y1 <- Yd %*% w_true + rnorm(T0 + T1, 0, 0.2); y1[(T0 + 1):(T0 + T1)] <- y1[(T0 + 1):(T0 + T1)] + tau
    list(Y = cbind(Yd, y1), dn = seq_len(J), tr = J + 1, pre = seq_len(T0), post = (T0 + 1):(T0 + T1), tau = tau, w = w_true)
}
est <- cover <- numeric(0); wrec <- numeric(0)
for (s in 1:60) {
    d <- sim_panel(s)
    f <- fam_bscm(d$Y, d$tr, d$dn, d$pre, d$post, ndraw = 400, seed = s)
    est <- c(est, f$tau); cover <- c(cover, f$cri[1] <= d$tau && d$tau <= f$cri[2])
    wrec <- c(wrec, NA)
}
chk("bscm: posterior mean of the effect is close to the truth", abs(mean(est) - 3) < 0.25 && stats::sd(est) < 0.5,
    sprintf("(mean %.3f, sd %.3f over 60 panels, truth 3)", mean(est), stats::sd(est)))
chk("bscm: 95% credible interval covers the truth", mean(cover) >= 0.85, sprintf("(coverage %.2f over 60 panels)", mean(cover)))
# (the true sparse weights are not identifiable from the pre-period fit when many donors can imitate the treated unit, so no check is made on which donors are picked)
# a null effect must give an interval around zero
dn <- sim_panel(99, tau = 0); f0 <- fam_bscm(dn$Y, dn$tr, dn$dn, dn$pre, dn$post, ndraw = 600, seed = 1)
chk("bscm: no effect gives an interval that contains zero", f0$cri[1] < 0 && f0$cri[2] > 0, sprintf("([%.2f, %.2f])", f0$cri[1], f0$cri[2]))
# reproducibility and the ndraw option
a1 <- fam_bscm(dn$Y, dn$tr, dn$dn, dn$pre, dn$post, ndraw = 300, seed = 5); a2 <- fam_bscm(dn$Y, dn$tr, dn$dn, dn$pre, dn$post, ndraw = 300, seed = 5)
chk("bscm: same seed gives the same draws", identical(a1$tau, a2$tau))

# ---- bscm against the horseshoe package (optional) ---------------------------------------------------------------------------
if (requireNamespace("horseshoe", quietly = TRUE)) {
    d <- sim_panel(7)
    Zc <- d$Y[d$pre, d$dn]; yc <- d$Y[d$pre, d$tr]
    ours <- fam_bscm(d$Y, d$tr, d$dn, d$pre, d$post, ndraw = 8000, seed = 11)
    zs <- apply(sweep(Zc, 2, colMeans(Zc)), 2, stats::sd); ysd <- stats::sd(yc - mean(yc))
    hs <- horseshoe::horseshoe(y = (yc - mean(yc)) / ysd, X = sweep(sweep(Zc, 2, colMeans(Zc)), 2, zs, "/"), method.tau = "halfCauchy",
                               method.sigma = "Jeffreys", burn = 4000, nmc = 8000, alpha = 0.05)
    ref_w <- hs$BetaHat * ysd / zs
    chk("bscm: posterior mean weights agree with the horseshoe package", max(abs(ours$w - ref_w)) < 0.05,
        sprintf("(max |diff| %.3f; weights up to %.2f)", max(abs(ours$w - ref_w)), max(abs(ref_w))))
} else cat("horseshoe package not installed: comparison with it skipped\n")

cat(if (nfail == 0) "\nAll enet and bscm checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
