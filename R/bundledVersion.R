#' Report the version of the bundled quasar binary
#'
#' Reads the `VERSION` metadata file shipped alongside the platform-specific
#' bundled binary (`inst/mac_arm_bin/` on macOS ARM64, `inst/windows_bin/` on
#' Windows).  Returns `NA` when no bundled binary is present for the current
#' platform.
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
    si <- Sys.info()
    vfile <- if (identical(si[["sysname"]], "Darwin") &&
                     identical(si[["machine"]], "arm64")) {
        system.file("mac_arm_bin", "VERSION", package = "QuasarExperiment")
    } else if (identical(si[["sysname"]], "Windows")) {
        system.file("windows_bin", "VERSION", package = "QuasarExperiment")
    } else {
        ""
    }
    if (!nzchar(vfile))
        return(NA_character_)
    lines <- readLines(vfile)
    parts <- strsplit(lines[nzchar(lines)], ": ", fixed = TRUE)
    vals  <- vapply(parts, `[[`, character(1), 2L)
    names(vals) <- vapply(parts, `[[`, character(1), 1L)
    vals
}
