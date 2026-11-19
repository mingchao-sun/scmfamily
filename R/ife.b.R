# jamovi analysis class; the estimation engine lives in common.R.
ifeClass <- R6::R6Class(
    "ifeClass",
    inherit = ifeBase,
    private = list(
        .run = function() fam_run(self, "ife", function() private$.checkpoint(flush = FALSE)),
        .plotPaths = function(image, ...) fam_draw_paths(image$state),
        .plotGaps = function(image, ...) fam_draw_gaps(image$state)
    )
)
