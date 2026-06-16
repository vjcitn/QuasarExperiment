#' Add gene symbols to rowRanges
#'
#' Looks up HGNC gene symbols for Ensembl gene IDs in `rownames(x)` using an
#' [ensembldb::EnsDb] annotation object and adds them as a `gene_name` column
#' in `mcols(rowRanges(x))`.  Version suffixes (e.g. `.3` in
#' `ENSG00000100181.3`) are stripped before lookup.
#'
#' @param x A [QuasarExperiment] (or any
#'   [SummarizedExperiment::RangedSummarizedExperiment] with Ensembl IDs as
#'   `rownames`).
#' @param ensdb An [ensembldb::EnsDb] object.  Defaults to
#'   `EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79` if that package is installed.
#'
#' @return `x` with an additional `gene_name` metadata column in
#'   `rowRanges(x)`.
#'
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv"),
#'     genome      = "hg38"
#' )
#' if (requireNamespace("EnsDb.Hsapiens.v79", quietly = TRUE))
#'     qe <- addGeneSymbols(qe)
#'
#' @importFrom AnnotationDbi select
#' @importFrom GenomeInfoDb genome
#' @export
addGeneSymbols <- function(x, ensdb = NULL) {
    if (is.null(ensdb)) {
        if (!requireNamespace("EnsDb.Hsapiens.v79", quietly = TRUE))
            stop("Package 'EnsDb.Hsapiens.v79' is needed. Install with:\n",
                 "  BiocManager::install('EnsDb.Hsapiens.v79')")
        ensdb <- EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79
    }

    ids <- rownames(x)
    ids_clean <- sub("\\..*$", "", ids)   # strip version suffix

    hits <- AnnotationDbi::select(
        ensdb,
        keys    = ids_clean,
        keytype = "GENEID",
        columns = c("GENEID", "GENENAME")
    )
    # one row per gene; match back in rowRanges order
    idx <- match(ids_clean, hits[["GENEID"]])
    mcols(rowRanges(x))[["gene_name"]] <- hits[["GENENAME"]][idx]
    x
}
