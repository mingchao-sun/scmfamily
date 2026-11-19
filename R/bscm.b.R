# jamovi analysis class; the estimation engine lives in common.R.
bscmClass <- R6::R6Class(
    "bscmClass",
    inherit = bscmBase,
    private = list(
        .run = function() fam_run(self, "bscm", function() private$.checkpoint(flush = FALSE)),
        .plotPaths = function(image, ...) fam_draw_paths(image$state),
        .plotGaps = function(image, ...) fam_draw_gaps(image$state)
    )
)
