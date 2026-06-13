#' Extract variant-level results from a quasar run
#'
#' @param x The list returned by [runQuasar()].
#' @return A [GenomicRanges::GRanges] with one entry per feature–variant
#'   pair, containing metadata columns `feature_id`, `snp_id`, `ref`, `alt`,
#'   `maf`, `beta`, `se`, and `pvalue`.
#' @seealso [runQuasar()], [regionHits()]
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' if (nzchar(Sys.which("quasar"))) {
#'     res <- runQuasar(qe, model = "lm", mode = "cis")
#'     variantHits(res)
#' }
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
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' if (nzchar(Sys.which("quasar"))) {
#'     res <- runQuasar(qe, model = "lm", mode = "cis")
#'     regionHits(res)
#' }
#' @export
regionHits <- function(x) {
    if (!is.list(x) || is.null(x[["regions"]]))
        stop("'x' must be the list returned by runQuasar()")
    x[["regions"]]
}
