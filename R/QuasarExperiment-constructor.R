#' Construct a QuasarExperiment from quasar input files
#'
#' Reads phenotype, covariate, variant metadata, and (optionally) GRM files
#' into a [QuasarExperiment] object. Genotype data is represented lazily via
#' [BEDMatrix::BEDMatrix]; the binary PLINK `.bed` file is not loaded into
#' memory.
#'
#' @param plinkPrefix Path prefix for the PLINK file set (`.bed`, `.bim`,
#'   `.fam`). May be an absolute or relative path.
#' @param phenoFile Path to the phenotype BED file. The file is tab-separated
#'   with columns `#chr`, `start`, `end`, `phenotype_id`, then one column per
#'   sample. Coordinates are 0-based half-open (BED convention).
#' @param covFile Path to the covariate TSV file. The file is tab-separated
#'   with a header row; the first column (`sample_id`) identifies samples and
#'   subsequent columns are covariates.
#' @param grmFile Optional path to the GRM TSV file. The file is
#'   tab-separated: the first column is `sample_id` and remaining columns are
#'   the same sample IDs in the header. Pass `NULL` (default) to omit the GRM.
#' @param assayName Name to assign the phenotype matrix in `assays()`.
#'   Defaults to `"pheno"`.
#'
#' @return A [QuasarExperiment] object.
#'
#' @details
#' **Sample ordering**: samples are ordered by their appearance in the
#' phenotype file. Covariates and the GRM are reordered to match. The
#' `BEDMatrix` column order follows the `.fam` file; a mapping is stored in
#' `colData(qe)$fam_index` so that the correct genotype columns can be
#' extracted if needed.
#'
#' **BED coordinates**: the phenotype BED file uses 0-based half-open
#' coordinates. These are converted to 1-based closed ranges for
#' `rowRanges()`, consistent with Bioconductor conventions.
#'
#' @examples
#' \dontrun{
#' exdir <- system.file("extdata", package = "QuasarExperiment")
#' qe <- QuasarExperiment(
#'     plinkPrefix = file.path(exdir, "chr22-n100"),
#'     phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
#'     covFile     = file.path(exdir, "cov-n100.tsv"),
#'     grmFile     = file.path(exdir, "grm-n100.tsv")
#' )
#' qe
#' }
#'
#' @export
QuasarExperiment <- function(plinkPrefix, phenoFile, covFile,
                             grmFile = NULL, assayName = "pheno") {
    # ---- phenotype BED --------------------------------------------------
    pheno_raw <- read.table(phenoFile, header = TRUE, sep = "\t",
                            check.names = FALSE, comment.char = "")
    colnames(pheno_raw)[1] <- "chr"

    sample_ids <- colnames(pheno_raw)[-(1:4)]
    n_samples  <- length(sample_ids)
    feature_ids <- pheno_raw[["phenotype_id"]]

    pheno_mat <- as.matrix(pheno_raw[, -(1:4), drop = FALSE])
    rownames(pheno_mat) <- feature_ids
    mode(pheno_mat) <- "double"

    # rowRanges: BED is 0-based half-open → 1-based closed
    row_gr <- GRanges(
        seqnames = pheno_raw[["chr"]],
        ranges   = IRanges::IRanges(
            start = pheno_raw[["start"]] + 1L,
            end   = pheno_raw[["end"]]
        ),
        phenotype_id = feature_ids
    )
    names(row_gr) <- feature_ids

    # ---- covariates -----------------------------------------------------
    cov_raw <- read.table(covFile, header = TRUE, sep = "\t",
                          check.names = FALSE)
    rownames(cov_raw) <- cov_raw[["sample_id"]]
    cov_raw[["sample_id"]] <- NULL

    # reorder to match phenotype sample order
    missing_cov <- setdiff(sample_ids, rownames(cov_raw))
    if (length(missing_cov))
        stop("Samples in phenotype file missing from covariate file: ",
             paste(missing_cov, collapse = ", "))
    cov_df <- S4Vectors::DataFrame(cov_raw[sample_ids, , drop = FALSE],
                                   row.names = sample_ids)

    # ---- PLINK / BEDMatrix ----------------------------------------------
    bed_file <- paste0(plinkPrefix, ".bed")
    bim_file <- paste0(plinkPrefix, ".bim")
    fam_file <- paste0(plinkPrefix, ".fam")
    for (f in c(bed_file, bim_file, fam_file))
        if (!file.exists(f))
            stop("PLINK file not found: ", f)

    geno_mat <- BEDMatrix::BEDMatrix(plinkPrefix)

    # store index mapping phenotype samples → fam rows in colData
    fam <- read.table(fam_file, header = FALSE,
                      col.names = c("fid", "iid", "pat", "mat", "sex", "phen"))
    fam_ids <- fam[["iid"]]
    fam_idx <- match(sample_ids, fam_ids)
    missing_geno <- sample_ids[is.na(fam_idx)]
    if (length(missing_geno))
        warning("Samples in phenotype file not found in .fam file: ",
                paste(missing_geno, collapse = ", "))
    cov_df[["fam_index"]] <- fam_idx

    # ---- variant ranges (from .bim) -------------------------------------
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

    # ---- GRM (optional) -------------------------------------------------
    grm_mat <- NULL
    if (!is.null(grmFile)) {
        grm_raw <- read.table(grmFile, header = TRUE, sep = "\t",
                              check.names = FALSE)
        rownames(grm_raw) <- grm_raw[["sample_id"]]
        grm_raw[["sample_id"]] <- NULL
        grm_raw <- grm_raw[sample_ids, sample_ids, drop = FALSE]
        grm_mat <- Matrix::Matrix(as.matrix(grm_raw), forceSymmetric = TRUE)
    }

    # ---- assemble -------------------------------------------------------
    se <- SummarizedExperiment::SummarizedExperiment(
        assays   = S4Vectors::SimpleList(setNames(list(pheno_mat), assayName)),
        rowRanges = row_gr,
        colData  = cov_df
    )

    new("QuasarExperiment",
        se,
        geno          = geno_mat,
        variantRanges = var_gr,
        grm           = grm_mat,
        plinkPrefix   = as.character(plinkPrefix)
    )
}
