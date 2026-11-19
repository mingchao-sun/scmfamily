# jamovi analysis class; estimation is in multisynth_engine.R (a port of the multisynth part of
# the 'augsynth' package) and stagascm_run.R. The renderer is written out here (not referenced
# at load time) because R sources files alphabetically and the shared code loads after this file.
stagascmClass <- R6::R6Class(
    "stagascmClass",
    inherit = stagascmBase,
    private = list(
        .run = function() stagascm_run(self, function() private$.checkpoint(flush = FALSE)),
        .plotEvent = function(image, ...) stagascm_draw_event(image$state)
    )
)
