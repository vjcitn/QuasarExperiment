#' Subset a QuasarExperiment
#'
#' Subsetting propagates to all slots: the inherited
#' [SummarizedExperiment::RangedSummarizedExperiment] slots (assay, rowRanges,
#' colData), the lazy `geno` matrix (samples × variants), and the `grm`.
#' The `variantRanges` and PLINK file are not subsetted — variant-level
#' subsetting is left to quasar's cis-window logic.
#'
#' @param x A [QuasarExperiment].
#' @param i Row (feature) index — integer, logical, or character.
#' @param j Column (sample) index — integer, logical, or character.
#' @param ... Ignored.
#' @param drop Ignored (kept for S4 compatibility).
#'
#' @return A [QuasarExperiment] with the selected features and/or samples.
#'
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' qe[1:5, 1:10]
#'
#' @export
setMethod("[", "QuasarExperiment", function(x, i, j, ..., drop = FALSE) {
    # delegate SE subsetting (handles missing i or j gracefully)
    se_sub <- if (missing(i) && missing(j))
        callNextMethod()
    else if (missing(i))
        callNextMethod(x, , j, ..., drop = drop)
    else if (missing(j))
        callNextMethod(x, i, , ..., drop = drop)
    else
        callNextMethod(x, i, j, ..., drop = drop)

    # subset geno rows (samples) to match the new colData
    j_idx <- if (missing(j)) seq_len(ncol(x)) else j
    new_geno <- x@geno[j_idx, , drop = FALSE]

    # subset grm if present
    new_grm <- if (!is.null(x@grm)) x@grm[j_idx, j_idx, drop = FALSE] else NULL

    se_sub@geno <- new_geno
    se_sub@grm  <- new_grm
    se_sub
})
