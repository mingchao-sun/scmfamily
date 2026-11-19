# jamovi analysis class; the estimation engine lives in common.R.
enetClass <- R6::R6Class(
    "enetClass",
    inherit = enetBase,
    private = list(
        .run = function() fam_run(self, "enet", function() private$.checkpoint(flush = FALSE)),
        .plotPaths = function(image, ...) fam_draw_paths(image$state),
        .plotGaps = function(image, ...) fam_draw_gaps(image$state)
    )
)
