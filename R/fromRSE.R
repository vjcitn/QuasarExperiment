#' Construct a QuasarExperiment from a RangedSummarizedExperiment
#'
#' Converts an existing [SummarizedExperiment::RangedSummarizedExperiment] into
#' a [QuasarExperiment] by attaching PLINK genotype data and, optionally, a
#' genetic relatedness matrix.  This is the recommended entry point for users
#' who already hold their pseudobulk or bulk RNA-seq phenotype data in a
#' Bioconductor SE object (e.g. from `tximeta`, `muscat`, or `DESeq2`).
#'
#' @param se A [SummarizedExperiment::RangedSummarizedExperiment] whose
#'   `rowRanges()` carry genomic coordinates for each feature and whose
#'   `colnames()` match the individual IDs (IID) in the PLINK `.fam` file.
#' @param plinkPrefix Path prefix for the PLINK file set (`.bed`, `.bim`,
#'   `.fam`).
#' @param covColumns Character vector naming columns of `colData(se)` to treat
#'   as covariates.  These columns are written to the covariate file passed to
#'   quasar.  Defaults to all columns of `colData(se)` except `fam_index` (if
#'   present).
#' @param grmFile Optional path to a GRM TSV file (same format as accepted by
#'   [QuasarExperiment()]).  Pass `NULL` (default) to omit the GRM.
#' @param assayName Name of the assay in `se` to use as the phenotype matrix.
#'   Defaults to the first assay.
#' @param featureIdColumn Name of a metadata column in `rowRanges(se)` to use
#'   as the `phenotype_id`.  If `NULL` (default), `rownames(se)` are used.
#'
#' @return A [QuasarExperiment] object.
#'
#' @details
#' **Sample matching**: `colnames(se)` are matched against the IID column of
#' the `.fam` file.  A warning is issued for any samples present in `se` but
#' absent from the `.fam`.
#'
#' **Coordinates**: `rowRanges(se)` must be a `GRanges` with valid `seqnames`,
#' `start`, and `end`.  These are used directly as the feature coordinates;
#' no coordinate conversion is performed (Bioconductor uses 1-based closed
#' ranges throughout).
#'
#' @seealso [QuasarExperiment()] for construction from flat files,
#'   [runQuasar()] to invoke the quasar binary.
#'
#' @examples
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#'
#' # Build a minimal RangedSummarizedExperiment from the example data
#' qe_ref <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv")
#' )
#' se <- as(qe_ref, "RangedSummarizedExperiment")
#'
#' # Reconstruct from the SE
#' qe2 <- QuasarExperimentFromRSE(
#'     se          = se,
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     covColumns  = c("int", "sex", "age", "expr_pc1", "expr_pc2",
#'                     "geno_pc1", "geno_pc2", "geno_pc3",
#'                     "geno_pc4", "geno_pc5", "geno_pc6")
#' )
#' qe2
#'
#' @export
QuasarExperimentFromRSE <- function(se, plinkPrefix,
                                    covColumns      = NULL,
                                    grmFile         = NULL,
                                    assayName       = NULL,
                                    featureIdColumn = NULL) {

    if (!is(se, "RangedSummarizedExperiment"))
        stop("'se' must be a RangedSummarizedExperiment")

    # ---- assay --------------------------------------------------------------
    if (is.null(assayName))
        assayName <- SummarizedExperiment::assayNames(se)[1]
    if (!assayName %in% SummarizedExperiment::assayNames(se))
        stop("assayName '", assayName, "' not found in assays(se)")

    pheno_mat <- as.matrix(SummarizedExperiment::assay(se, assayName))
    mode(pheno_mat) <- "double"

    sample_ids  <- colnames(se)
    feature_ids <- if (!is.null(featureIdColumn)) {
        mcols(SummarizedExperiment::rowRanges(se))[[featureIdColumn]]
    } else {
        rownames(se)
    }
    if (is.null(feature_ids) || anyNA(feature_ids))
        stop("Could not determine feature IDs: set rownames(se) or supply ",
             "'featureIdColumn'")
    rownames(pheno_mat) <- feature_ids

    # ---- rowRanges ----------------------------------------------------------
    rr <- SummarizedExperiment::rowRanges(se)
    if (!is(rr, "GRanges") || length(rr) == 0L)
        stop("rowRanges(se) must be a non-empty GRanges")
    names(rr) <- feature_ids
    if (!"phenotype_id" %in% names(S4Vectors::mcols(rr)))
        rr$phenotype_id <- feature_ids

    # ---- colData / covariates -----------------------------------------------
    cd <- SummarizedExperiment::colData(se)
    if (is.null(covColumns)) {
        covColumns <- setdiff(colnames(cd), "fam_index")
    } else {
        missing_cols <- setdiff(covColumns, colnames(cd))
        if (length(missing_cols))
            stop("covColumns not found in colData(se): ",
                 paste(missing_cols, collapse = ", "))
    }
    cov_df <- cd[, covColumns, drop = FALSE]

    # ---- PLINK / BEDMatrix --------------------------------------------------
    bim_file <- paste0(plinkPrefix, ".bim")
    fam_file <- paste0(plinkPrefix, ".fam")
    for (f in c(paste0(plinkPrefix, ".bed"), bim_file, fam_file))
        if (!file.exists(f))
            stop("PLINK file not found: ", f)

    geno_mat <- BEDMatrix::BEDMatrix(plinkPrefix)

    fam     <- read.table(fam_file, header = FALSE,
                          col.names = c("fid", "iid", "pat", "mat", "sex", "phen"))
    fam_idx <- match(sample_ids, fam[["iid"]])
    missing_geno <- sample_ids[is.na(fam_idx)]
    if (length(missing_geno))
        warning("Samples in se not found in .fam file: ",
                paste(missing_geno, collapse = ", "))
    cov_df[["fam_index"]] <- fam_idx

    # ---- variant ranges (from .bim) -----------------------------------------
    bim <- read.table(bim_file, header = FALSE,
                      col.names = c("chrom", "snp_id", "cm", "pos", "alt", "ref"),
                      colClasses = c("character", "character", "numeric",
                                     "integer",   "character", "character"))
    var_gr <- GRanges(
        seqnames = bim[["chrom"]],
        ranges   = IRanges::IRanges(start = bim[["pos"]], width = 1L),
        snp_id   = bim[["snp_id"]],
        ref      = bim[["ref"]],
        alt      = bim[["alt"]]
    )
    names(var_gr) <- bim[["snp_id"]]

    # ---- GRM (optional) -----------------------------------------------------
    grm_mat <- NULL
    if (!is.null(grmFile)) {
        grm_raw <- read.table(grmFile, header = TRUE, sep = "\t",
                              check.names = FALSE)
        rownames(grm_raw) <- grm_raw[["sample_id"]]
        grm_raw[["sample_id"]] <- NULL
        grm_raw  <- grm_raw[sample_ids, sample_ids, drop = FALSE]
        grm_mat  <- Matrix::forceSymmetric(Matrix::Matrix(as.matrix(grm_raw)))
    }

    # ---- assemble -----------------------------------------------------------
    se_new <- SummarizedExperiment::SummarizedExperiment(
        assays    = S4Vectors::SimpleList(setNames(list(pheno_mat), assayName)),
        rowRanges = rr,
        colData   = cov_df
    )

    new("QuasarExperiment",
        se_new,
        geno          = geno_mat,
        variantRanges = var_gr,
        grm           = grm_mat,
        plinkPrefix   = as.character(plinkPrefix)
    )
}
