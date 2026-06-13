#' Extract variant-level results from a quasar run
#'
#' @param x The list returned by [runQuasar()].
#' @return A [GenomicRanges::GRanges] with one entry per feature–variant
#'   pair, containing metadata columns `feature_id`, `snp_id`, `ref`, `alt`,
#'   `maf`, `beta`, `se`, and `pvalue`.
#' @seealso [runQuasar()], [regionHits()]
#' @export
variantHits <- function(x) {
    if (!is.list(x) || is.null(x[["variants"]]))
        stop("'x' must be the list returned by runQuasar()")
    x[["variants"]]
}

#' Extract region-level results from a quasar run
#'
#' @param x The list returned by [runQuasar()].
#' @return A [GenomicRanges::GRanges] with one entry per phenotype feature,
#'   containing metadata columns `feature_id` and `pvalue`.
#' @seealso [runQuasar()], [variantHits()]
#' @export
regionHits <- function(x) {
    if (!is.list(x) || is.null(x[["regions"]]))
        stop("'x' must be the list returned by runQuasar()")
    x[["regions"]]
}
