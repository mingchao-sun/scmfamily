# Agreement with the reference packages (computed live; each block is skipped if its package is missing).
# Run from the module root:  Rscript tests/test_reference.R
#   synthdid   (GitHub synth-inference/synthdid)  -> sdid, sc, did estimates and weights
#   augsynth   (GitHub ebenmichael/augsynth)      -> ridge ASCM and multisynth
#   Synth      (CRAN)                             -> classic SCM
for (f in c("R/common.R", "R/scm_engine.R", "R/sdid_engine.R", "R/augsynth_engine.R", "R/multisynth_engine.R")) source(f)
nfail <- 0
chk <- function(name, got, want, tol) {
    d <- max(abs(got - want))
    ok <- is.finite(d) && d < tol
    cat(sprintf("%-52s %s   max |diff| = %.2e (tol %.0e)\n", name, if (ok) "ok  " else "FAIL", d, tol))
    if (!ok) nfail <<- nfail + 1
}
have <- function(p) requireNamespace(p, quietly = TRUE)
pp <- read.csv("data/prop99_long.csv")
pq <- fam_panel(pp, "cigsale", "state", "year", "California", 1989)

# ---- synthdid: California Proposition 99 ---------------------------------------------------
if (have("synthdid")) {
    suppressMessages(library(synthdid))
    data("california_prop99")
    st <- panel.matrices(california_prop99)
    Ys <- t(pq$Y[, c(pq$donors, pq$treated)])[rownames(st$Y), ]
    ref <- list(sdid = synthdid_estimate(st$Y, st$N0, st$T0), sc = sc_estimate(st$Y, st$N0, st$T0),
                did = did_estimate(st$Y, st$N0, st$T0))
    our <- list(sdid = sdid_estimate(Ys, st$N0, st$T0), sc = sdid_sc_estimate(Ys, st$N0, st$T0),
                did = sdid_did_estimate(Ys, st$N0, st$T0))
    for (k in names(ref)) chk(paste("synthdid estimate:", k), as.numeric(our[[k]]), as.numeric(ref[[k]]), 1e-6)
    chk("synthdid unit weights (omega)",  attr(our$sdid, "weights")$omega,  attr(ref$sdid, "weights")$omega,  1e-6)
    chk("synthdid time weights (lambda)", attr(our$sdid, "weights")$lambda, attr(ref$sdid, "weights")$lambda, 1e-6)
    chk("synthdid effect curve", sdid_effect_curve(our$sdid), as.numeric(synthdid_effect_curve(ref$sdid)), 1e-6)
} else cat("synthdid not installed: skipped\n")

# ---- augsynth: single treated unit, ridge ASCM ------------------------------------------------
if (have("augsynth") && have("dplyr")) {
    suppressMessages(library(augsynth))
    d <- pp; d$trt <- as.integer(d$state == "California" & d$year >= 1989)
    msg <- capture.output(ra <- augsynth(cigsale ~ trt, state, year, d, progfunc = "Ridge", scm = TRUE))
    fa <- asyn_fit(pq$Y, pq$times, pq$treated, 1989)
    chk("augsynth ridge: lambda", fa$lambda, ra$lambda, 1e-6)
    chk("augsynth ridge: weights", fa$weights[, 1], ra$weights[names(fa$weights[, 1]), 1], 1e-6)
    chk("augsynth ridge: ATT", mean(asyn_predict(fa, att = TRUE)[pq$post]),
        summary(ra, inf_type = "jackknife")$average_att$Estimate, 1e-6)
    # no augmentation (plain SCM weights through the same QP)
    d2 <- capture.output(rn <- augsynth(cigsale ~ trt, state, year, d, progfunc = "None", scm = TRUE))
    fn <- asyn_fit(pq$Y, pq$times, pq$treated, 1989, progfunc = "none")
    chk("augsynth none: weights", fn$weights[, 1], rn$weights[names(fn$weights[, 1]), 1], 1e-6)
}

# ---- multisynth: staggered adoption (minimum-wage panel) -----------------------------------------
if (have("augsynth") && have("dplyr") && file.exists("data/mpdta_minwage.csv")) {
    mp <- read.csv("data/mpdta_minwage.csv")
    mp$county <- as.character(mp$county)
    mp <- mp[is.na(mp$first_treat) | !mp$first_treat %in% c(0, 2004), ]            # the 2004 cohort has one pre-period; fixed effects need two
    msg <- capture.output(rm <- multisynth(lemp ~ treated, county, year, mp[, c("county", "year", "lemp", "treated")], n_leads = 3, nu = 0.5, fixedeff = TRUE))
    Ym <- tapply(mp$lemp, list(mp$year, mp$county), sum)
    ad <- tapply(mp$first_treat, mp$county, function(z) z[1]); ad[ad == 0] <- NA
    fm <- ms_fit(Ym, as.numeric(rownames(Ym)), ad, n_leads = 3, nu = 0.5)
    ra <- summary(rm, inf_type = "none")$att
    ours <- ms_predict(fm, att = TRUE)
    ra_all <- ra[ra$Level == "Average", ]                        # event times -4..1, then the overall average
    chk("multisynth: average effect by event time + overall", ours[seq_len(nrow(ra_all)), 1], ra_all$Estimate, 1e-4)
}

# ---- Synth: classic SCM with the outcome in every pre-period as predictor -------------------------
if (have("Synth")) {
    suppressMessages(library(Synth))
    pp$id <- as.numeric(factor(pp$state)); ct <- unique(pp$id[pp$state == "California"])
    dp <- dataprep(foo = pp, predictors = NULL, dependent = "cigsale", unit.variable = "id", time.variable = "year",
                   special.predictors = lapply(1970:1988, function(y) list("cigsale", y, "mean")),
                   treatment.identifier = ct, controls.identifier = setdiff(unique(pp$id), ct),
                   time.predictors.prior = 1970:1988, time.optimize.ssr = 1970:1988,
                   unit.names.variable = "state", time.plot = 1970:2000)
    so <- synth(dp, verbose = FALSE)
    gap <- dp$Y1plot - dp$Y0plot %*% so$solution.w
    ref_att <- mean(gap[as.character(1989:2000), ])
    Xp <- scm_predictor_matrix(scm_build_panel(pp, "cigsale", "state", "year", "California", 1989)$df, "state", "year",
                               character(0), pq$donors, pq$times[pq$pre])
    pn <- scm_build_panel(pp, "cigsale", "state", "year", "California", 1989)
    Xp <- scm_predictor_matrix(pn$df, "state", "year", character(0), pn$units, pn$times[pn$pre])
    fo <- scm_fit(pn$Y, treat = "California", pre = pn$pre, post = pn$post, X = Xp, starts = 5, seed = 1)
    cat(sprintf("Synth: ATT %.4f, pre-RMSPE %.4f   module: ATT %.4f, pre-RMSPE %.4f\n", ref_att,
                sqrt(mean(gap[as.character(1970:1988), ]^2)), scm_att(fo), fo$rmspe_pre))
    chk("Synth: ATT", scm_att(fo), ref_att, 0.5)
    chk("Synth: pre-treatment RMSPE (module must fit at least as well)", max(0, fo$rmspe_pre - sqrt(mean(gap[as.character(1970:1988), ]^2))), 0, 0.05)
}

cat(if (nfail == 0) "\nAll reference checks passed.\n" else sprintf("\n%d reference check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
