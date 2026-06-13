#' @import methods
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame SimpleList mcols
#' @importFrom Matrix Matrix forceSymmetric
#' @importFrom BEDMatrix BEDMatrix
#' @importFrom utils head read.table write.table
#' @importFrom stats setNames
NULL

#' QuasarExperiment: a SummarizedExperiment for quasar eQTL mapping
#'
#' `QuasarExperiment` extends [SummarizedExperiment::RangedSummarizedExperiment]
#' to hold all data needed for the quasar eQTL mapping tool:
#'
#' - **assays**: phenotype matrix (features × samples), e.g. mean or sum
#'   pseudobulk counts
#' - **rowRanges**: genomic coordinates for each phenotype feature (gene)
#' - **colData**: per-sample metadata including covariates
#' - **geno**: lazy [BEDMatrix::BEDMatrix] over the PLINK `.bed` file
#'   (variants × samples); data are not loaded until subscripted
#' - **variantRanges**: [GenomicRanges::GRanges] for each variant, with
#'   metadata columns `snp_id`, `ref`, `alt`
#' - **grm**: optional packed-symmetric genetic relatedness matrix
#'   (`dspMatrix` from the Matrix package), or `NULL`
#' - **plinkPrefix**: path prefix for the PLINK file set, passed directly to
#'   the quasar binary so the genotype data are never materialised in R
#'
#' @slot geno A [BEDMatrix::BEDMatrix] (samples × variants) or any
#'   matrix-like object.
#' @slot variantRanges A [GenomicRanges::GRanges] of length equal to
#'   `nrow(geno)`.
#' @slot grm A `dspMatrix` (Matrix package) of dimension samples × samples, or `NULL`.
#' @slot plinkPrefix A length-one character string giving the PLINK file
#'   prefix.
#'
#' @seealso [QuasarExperiment()] for the constructor, [runQuasar()] to invoke
#'   the quasar binary.
#'
#' @exportClass QuasarExperiment
setClass("QuasarExperiment",
    contains = "RangedSummarizedExperiment",
    representation(
        geno         = "ANY",
        variantRanges = "GRanges",
        grm          = "ANY",
        plinkPrefix  = "character"
    )
)

setValidity("QuasarExperiment", function(object) {
    msg <- character(0)

    # BEDMatrix is samples x variants: rows = samples, cols = variants
    ns_geno <- nrow(object@geno)
    ns_se   <- ncol(object)
    if (ns_geno != ns_se)
        msg <- c(msg, sprintf(
            "nrow(geno) (%d) must equal ncol(SummarizedExperiment) (%d)",
            ns_geno, ns_se))

    nv_geno <- ncol(object@geno)
    nv <- length(object@variantRanges)
    if (nv_geno != nv)
        msg <- c(msg, sprintf(
            "ncol(geno) (%d) must equal length(variantRanges) (%d)", nv_geno, nv))

    if (!is.null(object@grm)) {
        if (!is(object@grm, "Matrix") && !is.matrix(object@grm))
            msg <- c(msg, "'grm' must be a matrix or Matrix object, or NULL")
        else {
            nr <- nrow(object@grm)
            if (nr != ncol(object))
                msg <- c(msg, sprintf(
                    "nrow(grm) (%d) must equal ncol(experiment) (%d)",
                    nr, ncol(object)))
        }
    }

    if (length(object@plinkPrefix) != 1L || is.na(object@plinkPrefix))
        msg <- c(msg, "'plinkPrefix' must be a single non-NA character string")

    if (length(msg)) msg else TRUE
})
