# ---- geno ---------------------------------------------------------------

#' Access the lazy genotype matrix
#'
#' @param x A [QuasarExperiment].
#' @param value A [BEDMatrix::BEDMatrix] or matrix-like object
#'   (samples × variants).
#' @return `geno()` returns the `geno` slot (a `BEDMatrix`, samples × variants).
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' geno(qe)
#' geno(qe)[1:3, 1:5]
#' @export
setGeneric("geno", function(x) standardGeneric("geno"))

#' @rdname geno
#' @export
setGeneric("geno<-", function(x, value) standardGeneric("geno<-"))

#' @rdname geno
setMethod("geno", "QuasarExperiment", function(x) x@geno)

#' @rdname geno
setReplaceMethod("geno", "QuasarExperiment", function(x, value) {
    x@geno <- value
    validObject(x)
    x
})

# ---- variantRanges ------------------------------------------------------

#' Access variant genomic ranges
#'
#' @param x A [QuasarExperiment].
#' @param value A [GenomicRanges::GRanges] with length equal to `nrow(geno(x))`.
#' @return `variantRanges()` returns a `GRanges` with one range per variant.
#'   Metadata columns include `snp_id`, `ref`, and `alt`.
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' variantRanges(qe)
#' @export
setGeneric("variantRanges", function(x) standardGeneric("variantRanges"))

#' @rdname variantRanges
#' @export
setGeneric("variantRanges<-", function(x, value) standardGeneric("variantRanges<-"))

#' @rdname variantRanges
setMethod("variantRanges", "QuasarExperiment", function(x) x@variantRanges)

#' @rdname variantRanges
setReplaceMethod("variantRanges", "QuasarExperiment", function(x, value) {
    x@variantRanges <- value
    validObject(x)
    x
})

# ---- grm ----------------------------------------------------------------

#' Access the genetic relatedness matrix
#'
#' @param x A [QuasarExperiment].
#' @param value A square matrix or `Matrix` (samples × samples),
#'   or `NULL`.
#' @return `grm()` returns the GRM as a `dspMatrix` or base matrix, or `NULL` if none was provided.
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv"),
#'     grmFile     = file.path(exdir, "grm-n100.tsv")
#' )
#' dim(grm(qe))
#' qe_no_grm <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' is.null(grm(qe_no_grm))
#' @export
setGeneric("grm", function(x) standardGeneric("grm"))

#' @rdname grm
#' @export
setGeneric("grm<-", function(x, value) standardGeneric("grm<-"))

#' @rdname grm
setMethod("grm", "QuasarExperiment", function(x) x@grm)

#' @rdname grm
setReplaceMethod("grm", "QuasarExperiment", function(x, value) {
    x@grm <- value
    validObject(x)
    x
})

# ---- plinkPrefix --------------------------------------------------------

#' Access the PLINK file prefix
#'
#' @param x A [QuasarExperiment].
#' @param value A length-one character string.
#' @return `plinkPrefix()` returns the path prefix for the `.bed`/`.bim`/`.fam`
#'   file set.
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' plinkPrefix(qe)
#' @export
setGeneric("plinkPrefix", function(x) standardGeneric("plinkPrefix"))

#' @rdname plinkPrefix
#' @export
setGeneric("plinkPrefix<-", function(x, value) standardGeneric("plinkPrefix<-"))

#' @rdname plinkPrefix
setMethod("plinkPrefix", "QuasarExperiment", function(x) x@plinkPrefix)

#' @rdname plinkPrefix
setReplaceMethod("plinkPrefix", "QuasarExperiment", function(x, value) {
    x@plinkPrefix <- value
    validObject(x)
    x
})
