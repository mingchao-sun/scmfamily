# Headless smoke test of every jamovi analysis (the .b.R classes), not just the engines.
# Best run in jamovi's own R with only jamovi/bundled libraries, as the installed module sees them:
#   R_LIBS_USER=notthere R_LIBS="<modules>/scmfamily/R;<jamovi>/Resources/modules/jmv/R;<jamovi>/Resources/modules/base/R" \
#     "<jamovi>/Frameworks/R/bin/Rscript.exe" tests/smoke.R
suppressMessages({library(jmvcore); library(R6); library(ggplot2)})
src <- c("R/common.R", "R/guide.R", "R/method_notes.R", "R/scm_help_text.R", "R/scm_engine.R", "R/sdid_engine.R",
         "R/augsynth_engine.R", "R/multisynth_engine.R", "R/scm_plots.R", "R/ascm_run.R", "R/sdid_run.R", "R/stagascm_run.R",
         "R/compare_run.R")
for (f in src) if (file.exists(f)) source(f)
for (f in c(list.files("R", "[.]h[.]R$", full.names = TRUE), list.files("R", "[.]b[.]R$", full.names = TRUE))) source(f)
p99 <- read.csv("data/prop99_long.csv"); mp <- read.csv("data/mpdta_minwage.csv")
mp$county <- as.character(mp$county)
bad <- 0
run <- function(label, fn, data, ...) {
    a <- list(...)
    fm <- formals(fn)
    for (n in setdiff(names(fm)[vapply(fm, function(z) identical(z, quote(expr = )), logical(1))], c("data", names(a)))) a[n] <- list(NULL)
    r <- tryCatch({ do.call(fn, c(list(data = data), a)); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("%-26s %s\n", label, r)); if (r != "ok") bad <<- bad + 1
}
base <- list(outcome = "cigsale", unit = "state", time = "year", treatedUnit = "California", treatmentTime = 1989)
for (k in c("scm", "pensc", "demsc", "ascm", "enet", "bscm", "ife", "mc", "sdid", "compare")) {
    # cheap settings: no placebo / robustness loops, plots off where the option exists
    extra <- list()
    fn <- get(k)
    for (o in c("placebo", "loo", "timing", "preLength", "inTime", "conformal", "jackknife", "bootstrap")) if (o %in% names(formals(fn))) extra[[o]] <- FALSE
    do.call(run, c(list(k, fn, p99), base, extra))
    # again with the guide and method notes turned on
    do.call(run, c(list(paste(k, "+ notes"), fn, p99), base, extra, list(showGuide = TRUE, showHelp = TRUE)))
}
run("stagascm (guide settings)", stagascm, mp, outcome = "lemp", unit = "county", time = "year", adoption = "first_treat", fixedeff = FALSE, timeCohort = TRUE)
run("stagascm + notes", stagascm, mp, outcome = "lemp", unit = "county", time = "year", adoption = "first_treat", fixedeff = FALSE, timeCohort = TRUE, showGuide = TRUE, showHelp = TRUE)
cat(if (bad == 0) "\nAll analyses ran.\n" else sprintf("\n%d analysis run(s) failed.\n", bad))
if (bad > 0) quit(status = 1)
