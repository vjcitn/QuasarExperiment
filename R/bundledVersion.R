#' Report the version of the bundled quasar binary
#'
#' Reads the `VERSION` metadata file shipped alongside the platform-specific
#' bundled binary in `inst/mac_arm_bin/`.  Returns `NA` when no bundled binary
#' is present for the current platform.
#'
#' @return A named character vector with fields `quasar_version`,
#'   `source_commit`, `source_repo`, `build_date`, `platform`, and
#'   `dynamic_libs`, or `NA` if no bundled binary metadata is found.
#'
#' @examples
#' bundledQuasarVersion()
#'
#' @export
bundledQuasarVersion <- function() {
    vfile <- system.file("mac_arm_bin", "VERSION",
                         package = "QuasarExperiment")
    if (!nzchar(vfile))
        return(NA_character_)
    lines <- readLines(vfile)
    parts <- strsplit(lines[nzchar(lines)], ": ", fixed = TRUE)
    vals  <- vapply(parts, `[[`, character(1), 2L)
    names(vals) <- vapply(parts, `[[`, character(1), 1L)
    vals
}
