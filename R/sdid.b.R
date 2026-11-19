# jamovi analysis class; estimation is in sdid_engine.R (a port of the 'synthdid' package)
# and sdid_run.R. Renderers are written out here (not referenced at load time) because R
# sources files alphabetically and the shared code loads after this file.
sdidClass <- R6::R6Class(
    "sdidClass",
    inherit = sdidBase,
    private = list(
        .run = function() sdid_run(self, function() private$.checkpoint(flush = FALSE)),
        .plotPaths = function(image, ...) fam_draw_paths(image$state),
        .plotGaps = function(image, ...) fam_draw_gaps(image$state)
    )
)
