# jamovi analysis class; estimation is in augsynth_engine.R (a port of the single-period part
# of the 'augsynth' package) and ascm_run.R. Renderers are written out here (not referenced
# at load time) because R sources files alphabetically and the shared code loads after this file.
ascmClass <- R6::R6Class(
    "ascmClass",
    inherit = ascmBase,
    private = list(
        .run = function() ascm_run(self, function() private$.checkpoint(flush = FALSE)),
        .plotPaths = function(image, ...) fam_draw_paths(image$state),
        .plotGaps = function(image, ...) fam_draw_gaps(image$state)
    )
)
