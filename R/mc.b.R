# jamovi analysis class; the estimation engine lives in common.R.
# Renderers are written out here (not referenced at load time) because R
# sources files alphabetically and common.R loads after this file.
mcClass <- R6::R6Class(
    "mcClass",
    inherit = mcBase,
    private = list(
        .run = function() fam_run(self, "mc", function() private$.checkpoint(flush = FALSE)),
        .plotPaths = function(image, ...) fam_draw_paths(image$state),
        .plotGaps = function(image, ...) fam_draw_gaps(image$state)
    )
)
