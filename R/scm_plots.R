
# scm_plots.R ---------------------------------------------------------------
# Base-graphics figures for the classic synthetic control analysis, drawn in the
# same style as the other analyses of the module (see fam_draw_paths in common.R).
# Each function takes the slim state built in scm_engine.R and draws directly.

scm_plot_paths <- function(st, ylab = "Outcome") {
    yl <- range(c(st$treated, st$synth), na.rm = TRUE)
    op <- graphics::par(mar = c(6.5, 4.1, 4.1, 2.1), xpd = FALSE)
    on.exit(graphics::par(op))
    graphics::plot(st$years, st$treated, type = "l", lwd = 2, ylim = yl,
                   xlab = "Time", ylab = ylab, main = "Treated unit vs counterfactual")
    graphics::lines(st$years, st$synth, lty = 2, lwd = 2, col = "steelblue")
    graphics::abline(v = st$years[st$post1], lty = 3)
    graphics::par(xpd = NA)
    graphics::legend("bottom", inset = c(0, -0.33), horiz = TRUE,
                     c("Treated", "Counterfactual"), lty = c(1, 2), lwd = 2,
                     col = c("black", "steelblue"), bty = "n")
    invisible(TRUE)
}

# Gap plot. Placebo gaps are drawn only for placebos whose pre-treatment RMSPE
# is at most k times the treated unit's (Abadie et al., 2010 practice).
scm_plot_gaps <- function(st, ylab = "Gap (treated - synthetic)") {
    G <- NULL
    if (!is.null(st$gaps)) {
        keep <- st$pre_rmspe <= st$k * st$rmspe_pre
        G <- st$gaps[, keep, drop = FALSE]
    }
    rng <- range(c(st$gap, G), na.rm = TRUE)
    graphics::plot(st$years, st$gap, type = "n", ylim = rng, xlab = "Time",
                   ylab = ylab, main = "Gap: treated minus counterfactual")
    if (!is.null(G)) for (j in seq_len(ncol(G))) graphics::lines(st$years, G[, j], col = "grey70")
    graphics::lines(st$years, st$gap, lwd = 2)
    graphics::abline(h = 0, lty = 3)
    graphics::abline(v = st$years[st$post1], lty = 3)
    if (!is.null(G))
        graphics::mtext(sprintf("Grey lines: placebo gaps with pre-fit RMSPE <= %g x treated (%d shown).",
                                st$k, ncol(G)), side = 1, line = 4, cex = 0.7)
    invisible(TRUE)
}

# Donor weights as horizontal bars (weights = named numeric vector).
scm_plot_weights <- function(weights, top = 15) {
    w <- rev(utils::head(sort(weights, decreasing = TRUE), top))
    op <- graphics::par(mar = c(4.1, 8, 4.1, 2.1))
    on.exit(graphics::par(op))
    graphics::barplot(as.numeric(w), names.arg = names(w), horiz = TRUE, las = 1,
                      col = "steelblue", border = NA, xlab = "Weight", main = "Donor weights")
    invisible(TRUE)
}

# Leave-one-out dot plot: ATT when each donor is removed, with the full-sample ATT.
scm_plot_loo <- function(loo, att_main) {
    d <- loo[order(loo$att), ]
    op <- graphics::par(mar = c(6.5, 8, 4.1, 2.1))
    on.exit(graphics::par(op))
    graphics::plot(d$att, seq_len(nrow(d)), pch = 19, yaxt = "n", ylab = "",
                   xlim = range(c(d$att, att_main, 0), na.rm = TRUE),
                   xlab = "Average post-treatment gap", main = "Leave-one-out")
    graphics::axis(2, at = seq_len(nrow(d)), labels = as.character(d$removed), las = 1,
                   cex.axis = min(1, 22 / nrow(d)))
    graphics::abline(v = 0, lty = 3)
    graphics::abline(v = att_main, lty = 2, col = "steelblue", lwd = 2)
    graphics::mtext("Dashed line: full-sample ATT.", side = 1, line = 4, cex = 0.7)
    invisible(TRUE)
}

# Scenario plot for timing, pre-period length and in-time placebo checks.
scm_plot_scenario <- function(df, att_main, xlab) {
    graphics::plot(df$x, df$att, type = "b", pch = 19, xlab = xlab,
                   ylim = range(c(df$att, att_main, 0), na.rm = TRUE),
                   ylab = "Average post-treatment gap", main = "Robustness check")
    graphics::abline(h = 0, lty = 3)
    graphics::abline(h = att_main, lty = 2, col = "steelblue", lwd = 2)
    graphics::mtext("Dashed line: full-sample ATT.", side = 1, line = 4, cex = 0.7)
    invisible(TRUE)
}
