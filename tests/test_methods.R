# Run from the module root:  Rscript tests/test_methods.R
# Simulated panel: 21 units, 30 periods, treatment at period 21, true effect = 5.
source("R/common.R")
set.seed(42)
Tn <- 30; N <- 21; t0 <- 21
units <- paste0("u", 1:N)
f <- cumsum(rnorm(Tn, 0, 0.5))
lam <- rnorm(N)
alpha <- rnorm(N, 0, 2)
df <- expand.grid(time = 1:Tn, unit = units, stringsAsFactors = FALSE)
df$y <- alpha[match(df$unit, units)] + lam[match(df$unit, units)] * f[df$time] +
    rnorm(nrow(df), 0, 0.2)
df$y[df$unit == "u1" & df$time >= t0] <- df$y[df$unit == "u1" & df$time >= t0] + 5

pan <- fam_panel(df, "y", "unit", "time", "u1", t0)
stopifnot(length(pan$donors) == N - 1, length(pan$pre) == t0 - 1, length(pan$post) == Tn - t0 + 1)

simplex <- c("penalized", "demeaned")
for (m in c(simplex, "ife", "enet", "bscm", "mc")) {
    fit <- fam_fit(m, pan$Y, pan$treated, pan$donors, pan$pre, pan$post,
                   lambda = 0.01, nfac = 2, alpha = 0.5)
    stopifnot(is.finite(fit$att))
    if (m == "bscm") stopifnot(length(fit$cri) == 2, fit$cri[1] < fit$cri[2])
    if (m %in% simplex) stopifnot(abs(sum(fit$w) - 1) < 1e-6, all(fit$w >= -1e-10))
    cat(sprintf("%-10s ATT = %6.3f   pre RMSPE = %6.3f   ratio = %6.2f\n",
                m, fit$att, fit$rmspe_pre, fit$ratio))
    stopifnot(abs(fit$att - 5) < 2)
}

# Penalty lambda = 0 must give finite simplex weights.
f0 <- fam_fit("penalized", pan$Y, pan$treated, pan$donors, pan$pre, pan$post, lambda = 0)
stopifnot(all(is.finite(f0$w)))

# Unbalanced panels must be rejected.
bad <- df[-5, ]
stopifnot(inherits(try(fam_panel(bad, "y", "unit", "time", "u1", t0), silent = TRUE), "try-error"))
fp <- fam_fit("penalized", pan$Y, pan$treated, pan$donors, pan$pre, pan$post, lambda = 0.01)
pars <- list(lambda = 0.01, nfac = 1, alpha = 0.5, ndraw = 200)
cf <- fam_conformal("penalized", pan$Y, pan$treated, pan$donors, pan$pre, pan$post,
                    fp$att, fp$gap, pars, ngrid = 41)
stopifnot(cf$p0 >= 0, cf$p0 <= 1, all(cf$pv >= 0 & cf$pv <= 1))
cat(sprintf("conformal p(0) = %.3f   CI = [%.3f, %.3f]\n", cf$p0, cf$lower, cf$upper))
# ---- covariates -----------------------------------------------------------------
# y = alpha_i + lambda_i f_t + 2 * x_it + noise, effect 5. IFE must recover beta ~ 2.
set.seed(7)
x <- matrix(rnorm(Tn * N), Tn, N) + matrix(lam, Tn, N, byrow = TRUE)
df2 <- df
df2$x1 <- as.vector(x)[match(paste(df2$time, df2$unit), paste(rep(1:Tn, N), rep(units, each = Tn)))]
df2$y <- df2$y + 2 * df2$x1
df2$m1 <- as.numeric(factor(df2$unit)) / N       # time-invariant unit covariate
pc <- fam_panel(df2, "y", "unit", "time", "u1", t0, covars = c("x1", "m1"))
stopifnot(all(dim(pc$cov) == c(Tn, N, 2)))
for (m in c("penalized", "demeaned")) {
    f <- fam_fit(m, pc$Y, pc$treated, pc$donors, pc$pre, pc$post, cov = pc$cov, covw = 1)
    stopifnot(is.finite(f$att), abs(sum(f$w) - 1) < 1e-6, all(f$w >= -1e-10))
    f0 <- fam_fit(m, pc$Y, pc$treated, pc$donors, pc$pre, pc$post, cov = pc$cov, covw = 0)
    fn <- fam_fit(m, pc$Y, pc$treated, pc$donors, pc$pre, pc$post)
    stopifnot(max(abs(f0$w - fn$w)) < 1e-6)     # zero covariate weight = no covariates
    cat(sprintf("%-10s with covariates: ATT = %6.3f   (without: %6.3f)\n", m, f$att, fn$att))
}
pi <- fam_panel(df2, "y", "unit", "time", "u1", t0, covars = "x1")
fi <- fam_fit("ife", pi$Y, pi$treated, pi$donors, pi$pre, pi$post, nfac = 1, cov = pi$cov)
cat(sprintf("ife        beta(x1) = %5.3f (true 2)   ATT = %6.3f\n", fi$beta[["x1"]], fi$att))
stopifnot(abs(fi$beta[["x1"]] - 2) < 0.3, abs(fi$att - 5) < 1.5)
# a donor treated as "treated" (placebo call) and missing-value handling
fj <- fam_fit("penalized", pc$Y, "u2", setdiff(pc$donors, "u2"), pc$pre, pc$post, cov = pc$cov)
stopifnot(is.finite(fj$att))
df3 <- df2; df3$x1[5] <- NA
p3 <- fam_panel(df3, "y", "unit", "time", "u1", t0, covars = "x1")
stopifnot(is.finite(fam_fit("penalized", p3$Y, p3$treated, p3$donors, p3$pre, p3$post, cov = p3$cov)$att))
stopifnot(inherits(try(fam_fit("ife", p3$Y, p3$treated, p3$donors, p3$pre, p3$post, cov = p3$cov), silent = TRUE), "try-error"))
# ---- synthdid port --------------------------------------------------------------------
source("R/sdid_engine.R")
Ys <- t(pan$Y[, c(pan$donors, pan$treated)])
es <- sdid_estimate(Ys, length(pan$donors), length(pan$pre))
stopifnot(abs(as.numeric(es) - 5) < 1.5, abs(sum(attr(es, "weights")$omega) - 1) < 1e-8,
          abs(sum(attr(es, "weights")$lambda) - 1) < 1e-8)
stopifnot(abs(mean(sdid_effect_curve(es)) - as.numeric(es)) < 1e-10)
set.seed(1)
cat(sprintf("sdid       ATT = %6.3f   placebo SE = %5.3f\n", as.numeric(es), sdid_se(es, "placebo", 30)))
colnames(x) <- units
Xs <- array(t(x[, c(pan$donors, pan$treated)]), dim = c(N, Tn, 1))
esx <- sdid_estimate(Ys + 2 * Xs[, , 1], length(pan$donors), length(pan$pre), X = Xs)
cat(sprintf("sdid + cov beta = %5.3f (true 2 on x added to Y)\n", attr(esx, "weights")$beta))
# ---- augsynth port (needs osqp and Matrix) ----------------------------------------------
if (requireNamespace("osqp", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) {
    source("R/augsynth_engine.R")
    fa <- asyn_fit(pan$Y, pan$times, pan$treated, t0)
    stopifnot(abs(mean(asyn_predict(fa, att = TRUE)[(t0):Tn]) - 5) < 1.5, abs(sum(fa$weights) - 1) < 1e-6)
    fn <- asyn_fit(pan$Y, pan$times, pan$treated, t0, progfunc = "none")
    stopifnot(all(fn$weights >= -1e-6), abs(sum(fn$weights) - 1) < 1e-6)
    cat(sprintf("augsynth   ATT = %6.3f   lambda = %.3g\n", mean(asyn_predict(fa, att = TRUE)[t0:Tn]), fa$lambda))
}
# ---- multisynth port ---------------------------------------------------------------------
if (requireNamespace("osqp", quietly = TRUE) && requireNamespace("Matrix", quietly = TRUE)) {
    source("R/augsynth_engine.R"); source("R/multisynth_engine.R")
    set.seed(10)
    Ns <- 24; Ts <- 16
    us <- sprintf("u%02d", 1:Ns)
    adopt_s <- setNames(c(rep(9, 2), rep(12, 2), rep(NA, Ns - 4)), us)
    lam_s <- rnorm(Ns); f_s <- cumsum(rnorm(Ts, 0, .5))
    Ym <- sapply(1:Ns, function(i) rnorm(1, 0, 2) + lam_s[i] * f_s + 0.1 * (1:Ts) + rnorm(Ts, 0, .3) +
                     ifelse(!is.na(adopt_s[i]) & 1:Ts >= adopt_s[i], 3, 0))
    colnames(Ym) <- us
    fm <- ms_fit(Ym, 1:Ts, adopt_s)
    stopifnot(all(abs(colSums(fm$weights) - 1) < 1e-8), all(fm$weights >= 0))
    att_m <- ms_predict(fm, att = TRUE)
    cat(sprintf("multisynth ATT = %6.3f  nu = %.3f\n", att_m[nrow(att_m), 1], fm$nu))
    stopifnot(abs(att_m[nrow(att_m), 1] - 3) < 1.5)
    set.seed(1)
    sm <- ms_summary(fm, "bootstrap", n_boot = 30)
    stopifnot(all(sm$lower[!is.na(sm$lower)] <= sm$upper[!is.na(sm$upper)]))
}
# ---- classic synthetic control (copied from the synthcontrol module) --------------------
if (requireNamespace("quadprog", quietly = TRUE) && file.exists("data/prop99_long.csv")) {
    source("R/scm_engine.R")
    pp <- read.csv("data/prop99_long.csv")
    pn <- scm_build_panel(pp, "cigsale", "state", "year", "California", 1989)
    Xp <- scm_predictor_matrix(pn$df, "state", "year", character(0), pn$units, pn$times[pn$pre])
    fp0 <- scm_fit(pn$Y, treat = "California", pre = pn$pre, post = pn$post, X = Xp, starts = 2, seed = 1)
    cat(sprintf("scm        ATT = %6.3f   pre RMSPE = %5.3f
", scm_att(fp0), fp0$rmspe_pre))
    stopifnot(abs(scm_att(fp0) - (-19.51)) < 0.1, abs(sum(fp0$weights) - 1) < 1e-6)
}
# ---- agreement with reference implementations (California Proposition 99) ----------------
# Reference values were computed once with the packages named below (R 4.6, 2026-11):
# pensynth (penalized), augsynth with fixedeff = TRUE and progfunc = "none" (demeaned),
# gsynth (IFE, two-way fixed effects, CV = FALSE) and MCPanel::mcnnm_fit (matrix completion).
if (file.exists("data/prop99_long.csv")) {
    pp <- read.csv("data/prop99_long.csv")
    pq <- fam_panel(pp, "cigsale", "state", "year", "California", 1989)
    att_of <- function(m, ...) fam_fit(m, pq$Y, pq$treated, pq$donors, pq$pre, pq$post, ...)$att
    chk <- function(name, got, want, tol) {
        cat(sprintf("%-34s ours %9.5f   reference %9.5f\n", name, got, want))
        stopifnot(abs(got - want) < tol)
    }
    # pensynth's lambda is half of ours; both scale the outcomes by the SD across donors
    chk("penalized lambda 0.01", att_of("penalized", lambda = 0.01), -22.2606046, 5e-3)
    chk("penalized lambda 0.1",  att_of("penalized", lambda = 0.1),  -23.0896033, 5e-3)
    chk("penalized lambda 1",    att_of("penalized", lambda = 1),    -25.3583312, 5e-3)
    chk("demeaned (augsynth)",   att_of("demeaned"),                 -11.1090420, 1e-4)
    chk("ife r = 1 (gsynth)",    att_of("ife", nfac = 1),            -13.8912234, 1e-6)
    chk("ife r = 2 (gsynth)",    att_of("ife", nfac = 2),            -0.4042069, 1e-6)
    chk("ife r = 3 (gsynth)",    att_of("ife", nfac = 3),            -1.8727750, 1e-6)
    chk("matrix completion 0.1 (MCPanel)", att_of("mc", mcl = 0.1),  -20.0518881, 5e-3)
    chk("matrix completion 0.3 (MCPanel)", att_of("mc", mcl = 0.3),  -21.5666670, 5e-3)
}

# ---- elastic net: exact solution of the paper's objective ---------------------------------
if (requireNamespace("glmnet", quietly = TRUE) && file.exists("data/prop99_long.csv")) {
    en <- fam_enet(pq$Y, pq$treated, pq$donors, pq$pre, pq$post, alpha = 0.5, cv_alpha = FALSE)
    a <- en$extra$alpha; lam <- en$extra$lambda; X <- pq$Y[pq$pre, pq$donors]; y <- pq$Y[pq$pre, pq$treated]
    p <- ncol(X); xm <- colMeans(X); ym <- mean(y)
    Xa <- rbind(sweep(X, 2, xm), sqrt(lam * (1 - a) / 2) * diag(p))            # augmented-lasso identity
    fa <- suppressWarnings(glmnet::glmnet(Xa, c(y - ym, rep(0, p)), alpha = 1, standardize = FALSE, intercept = FALSE,
                                          lambda = lam * a / (2 * nrow(Xa)), thresh = 1e-14, maxit = 1e7))
    wx <- as.vector(as.matrix(coef(fa))[-1, 1])
    cat(sprintf("elastic net: max |weights - exact solution| = %.2e (weights up to %.3f)\n", max(abs(wx - en$w)), max(abs(wx))))
    stopifnot(max(abs(wx - en$w)) < 1e-4, abs(en$extra$intercept - (ym - sum(xm * wx))) < 1e-3)
}

# Compare Estimators: every method must give exactly the numbers of the underlying engine call.
for (f in c("R/scm_engine.R", "R/sdid_engine.R", "R/augsynth_engine.R", "R/compare_run.R")) source(f)
ctx <- list(Y = pan$Y, tr = pan$treated, dn = pan$donors, pre = pan$pre, post = pan$post,
            times = pan$times, t0 = t0)
co <- list(lambda = 0.01, nFactors = 2, mcLambda = 0.1, bayesDraws = 500, optStarts = 2, seed = 1)
for (k in names(compare_methods)) {
    r <- compare_one(k, ctx, co, FALSE, function() NULL)
    stopifnot(is.finite(r$att), length(r$path) == Tn, length(r$gap) == Tn, abs(r$att - 5) < 3)
    cat(sprintf("compare %-5s ATT = %6.3f\n", k, r$att))
}
r <- compare_one("pen", ctx, co, FALSE, function() NULL)
d <- fam_fit("penalized", pan$Y, pan$treated, pan$donors, pan$pre, pan$post, 0.01, 2, 0.5, 500, 0.1, NULL, 1, TRUE, NULL)
stopifnot(identical(r$att, d$att), identical(r$path, d$path))
r <- compare_one("scm", ctx, co, FALSE, function() NULL)
d <- scm_fit(pan$Y, pan$treated, pan$pre, pan$post, X = NULL, starts = 2, seed = 1)
stopifnot(identical(r$att, mean(d$gap[pan$post])))
r <- compare_one("sdid", ctx, co, FALSE, function() NULL)
Ys <- t(pan$Y[, c(pan$donors, pan$treated)])
stopifnot(abs(r$att - as.numeric(sdid_estimate(Ys, length(pan$donors), length(pan$pre)))) < 1e-12)

cat("All scmfamily engine checks passed.\n")
