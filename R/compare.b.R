# jamovi analysis class; the engine lives in compare_run.R (loaded after this file, so the
# renderers are written out here).
compareClass <- R6::R6Class(
    "compareClass",
    inherit = compareBase,
    private = list(
        .run = function() compare_run(self, function() private$.checkpoint(flush = FALSE)),
        .plotPaths = function(image, ...) compare_draw_paths(image$state),
        .plotGaps = function(image, ...) compare_draw_gaps(image$state),
        .plotAtt = function(image, ...) compare_draw_att(image$state)
    )
)
