#' Run the quasar eQTL mapping binary
#'
#' Writes temporary input files from a [QuasarExperiment], invokes the quasar
#' binary, reads the results, and returns them as a named list of
#' [GenomicRanges::GRanges] objects.
#'
#' @param x A [QuasarExperiment].
#' @param binary Path to the quasar executable. Defaults to `"quasar"` (i.e.
#'   resolved via `PATH`). Can be set globally with
#'   `options(QuasarExperiment.binary = "/path/to/quasar")`.
#' @param model One of `"lm"`, `"nb_glm"`, or `"lmm"`. Passed to
#'   `--model`.
#' @param mode One of `"cis"`, `"trans"`, or `"residualise"`. Passed to
#'   `--mode`.
#' @param assayName Name of the assay in `x` to use as the phenotype matrix.
#'   Defaults to the first assay.
#' @param useApl Logical; pass `--use-apl` flag (Cox-Reid adjusted profile
#'   likelihood for NB dispersion). Recommended for `model = "nb_glm"`.
#' @param outPrefix Character prefix for quasar output files written inside
#'   the temp directory. Defaults to `"quasar_run"`.
#' @param debug Logical; if `TRUE`, the temporary directory containing the
#'   intermediate input files (phenotype BED, covariate TSV, optional GRM) is
#'   not deleted after the run and its path is printed.  Use this to inspect
#'   exactly what is passed to the quasar binary.
#' @param ... Additional arguments passed verbatim to the quasar binary as
#'   `--key value` pairs (values coerced to character). Logical `TRUE` values
#'   emit the flag with no value.
#'
#' @return A list with two [GenomicRanges::GRanges] elements:
#'   \describe{
#'     \item{`variants`}{One range per feature–variant pair tested, with
#'       metadata columns `feature_id`, `snp_id`, `ref`, `alt`, `maf`,
#'       `beta`, `se`, `pvalue`.}
#'     \item{`regions`}{One range per feature (the lead signal), with
#'       metadata columns `feature_id`, `pvalue`.}
#'   }
#'
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
#'     regionHits(res)
#' }
#'
#' @export
setGeneric("runQuasar",
    function(x, binary = NULL, model = c("lm", "nb_glm", "lmm"),
             mode = c("cis", "trans", "residualise"),
             assayName = NULL, useApl = FALSE, outPrefix = "quasar_run",
             debug = FALSE, ...)
        standardGeneric("runQuasar")
)

#' @rdname runQuasar
setMethod("runQuasar", "QuasarExperiment",
    function(x, binary = NULL, model = c("lm", "nb_glm", "lmm"),
             mode = c("cis", "trans", "residualise"),
             assayName = NULL, useApl = FALSE, outPrefix = "quasar_run",
             debug = FALSE, ...) {

        model <- match.arg(model)
        mode  <- match.arg(mode)

        binary <- .resolve_binary(binary)
        if (!is.null(assayName)) {
            if (!assayName %in% SummarizedExperiment::assayNames(x))
                stop("assayName '", assayName, "' not found in assays(x)")
        } else {
            assayName <- SummarizedExperiment::assayNames(x)[1]
        }

        tmpdir <- tempfile("quasar_")
        dir.create(tmpdir)
        if (debug) {
            message("debug=TRUE: intermediate files kept in ", tmpdir)
        } else {
            on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
        }

        pheno_file <- file.path(tmpdir, "pheno.bed")
        cov_file   <- file.path(tmpdir, "cov.tsv")
        grm_file   <- if (!is.null(x@grm)) file.path(tmpdir, "grm.tsv") else NULL

        .write_pheno_bed(x, assayName, pheno_file)
        .write_cov_tsv(x, cov_file)
        if (!is.null(grm_file))
            .write_grm_tsv(x, grm_file)

        if (debug) {
            message("--- pheno.bed (first 2 lines) ---")
            message(paste(readLines(pheno_file, n = 2), collapse = "\n"))
            message("--- cov.tsv (first 2 lines) ---")
            message(paste(readLines(cov_file, n = 2), collapse = "\n"))
        }

        extra <- list(...)
        args <- c(
            "-p", x@plinkPrefix,
            "-b", pheno_file,
            "-c", cov_file,
            "--model", model,
            "--mode",  mode,
            "-o", file.path(tmpdir, outPrefix)
        )
        if (!is.null(grm_file))
            args <- c(args, "-g", grm_file)
        if (useApl)
            args <- c(args, "--use-apl")
        for (nm in names(extra)) {
            val <- extra[[nm]]
            if (isTRUE(val))
                args <- c(args, paste0("--", nm))
            else
                args <- c(args, paste0("--", nm), as.character(val))
        }

        # Temporarily unset KMP/OMP thread-limit variables that R or its
        # packages may have set, so quasar's OpenMP can use all cores.
        omp_limit_vars <- c("OMP_THREAD_LIMIT", "KMP_DEVICE_THREAD_LIMIT",
                            "KMP_TEAMS_THREAD_LIMIT", "KMP_ALL_THREADS")
        saved_env <- Sys.getenv(omp_limit_vars, names = TRUE, unset = NA_character_)
        vars_to_unset <- names(saved_env)[!is.na(saved_env)]
        if (length(vars_to_unset)) {
            Sys.unsetenv(vars_to_unset)
            on.exit(
                do.call(Sys.setenv, as.list(saved_env[vars_to_unset])),
                add = TRUE
            )
        }

        status <- system2(binary, args = args)
        if (status != 0L)
            stop("quasar exited with status ", status)

        .read_results(tmpdir, outPrefix, mode, x)
    }
)

# ---- internal helpers ---------------------------------------------------

.resolve_binary <- function(binary) {
    if (!is.null(binary))
        return(binary)

    # honour explicit user preference
    opt <- getOption("QuasarExperiment.binary")
    if (!is.null(opt))
        return(opt)

    # use bundled binary for known platforms
    si <- Sys.info()
    bundled <- if (identical(si[["sysname"]], "Darwin") &&
                       identical(si[["machine"]], "arm64")) {
        system.file("mac_arm_bin", "quasar", package = "QuasarExperiment")
    } else if (identical(si[["sysname"]], "Windows")) {
        system.file("windows_bin", "quasar.exe", package = "QuasarExperiment")
    } else {
        ""
    }
    if (nzchar(bundled) && file.exists(bundled)) {
        if (!identical(si[["sysname"]], "Windows"))
            Sys.chmod(bundled, "0755")
        return(bundled)
    }

    # fall back to PATH
    path_hit <- Sys.which("quasar")
    if (nzchar(path_hit))
        return(path_hit)

    stop(
        "quasar binary not found. ",
        "On macOS ARM64 a bundled binary is used automatically. ",
        "On other platforms, install quasar and ensure it is on PATH, ",
        "or set options(QuasarExperiment.binary = '/path/to/quasar')."
    )
}

.write_pheno_bed <- function(x, assayName, path) {
    rr  <- SummarizedExperiment::rowRanges(x)
    mat <- SummarizedExperiment::assay(x, assayName)
    # quasar expects integer chromosome numbers using PLINK coding:
    # autosomes 1-22, X=23, Y=24, XY=25, MT/M=26
    chr_raw <- sub("^chr", "", as.character(GenomicRanges::seqnames(rr)))
    chr_raw[chr_raw == "X"]              <- "23"
    chr_raw[chr_raw == "Y"]              <- "24"
    chr_raw[chr_raw == "XY"]             <- "25"
    chr_raw[chr_raw %in% c("M", "MT")]   <- "26"
    chr <- suppressWarnings(as.integer(chr_raw))
    if (anyNA(chr)) {
        bad <- unique(as.character(GenomicRanges::seqnames(rr))[is.na(chr)])
        warning("Dropping ", sum(is.na(chr)), " features with seqnames that ",
                "cannot be mapped to PLINK integer codes: ",
                paste(bad, collapse = ", "))
        keep <- !is.na(chr)
        rr   <- rr[keep]
        mat  <- mat[keep, , drop = FALSE]
        chr  <- chr[keep]
    }
    # convert back to 0-based half-open BED coordinates
    df <- data.frame(
        `#chr`        = chr,
        start         = GenomicRanges::start(rr) - 1L,
        end           = GenomicRanges::end(rr),
        phenotype_id  = names(rr),
        check.names   = FALSE,
        stringsAsFactors = FALSE
    )
    out <- cbind(df, mat[, colnames(x), drop = FALSE])
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

.write_cov_tsv <- function(x, path) {
    cd <- as.data.frame(SummarizedExperiment::colData(x))
    cd[["fam_index"]] <- NULL
    non_numeric <- names(cd)[!vapply(cd, is.numeric, logical(1))]
    if (length(non_numeric))
        stop("Non-numeric covariate columns: ",
             paste(non_numeric, collapse = ", "),
             ". Use QuasarExperimentFromRSE() with a covariateFormula ",
             "to dummy-code categorical variables via model.matrix().")
    out <- cbind(data.frame(sample_id = rownames(cd),
                             stringsAsFactors = FALSE), cd)
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

.write_grm_tsv <- function(x, path) {
    g <- as.matrix(x@grm)
    ids <- colnames(x)
    rownames(g) <- ids
    colnames(g) <- ids
    out <- cbind(data.frame(sample_id = ids, stringsAsFactors = FALSE),
                 as.data.frame(g))
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

.read_results <- function(tmpdir, outPrefix, mode, x) {
    stem <- file.path(tmpdir, outPrefix)
    var_file    <- paste0(stem, "-quasar-", mode, "-variant.txt")
    region_file <- paste0(stem, "-quasar-", mode, "-region.txt")

    # variant-level results
    vdf <- read.table(var_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
    var_gr <- GRanges(
        seqnames = vdf[["chrom"]],
        ranges   = IRanges::IRanges(start = vdf[["pos"]], width = 1L),
        feature_id = vdf[["feature_id"]],
        snp_id     = vdf[["snp_id"]],
        ref        = vdf[["ref"]],
        alt        = vdf[["alt"]],
        maf        = vdf[["maf"]],
        beta       = vdf[["beta"]],
        se         = vdf[["se"]],
        pvalue     = vdf[["pvalue"]]
    )

    # region-level (per-feature) results
    rdf <- read.table(region_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
    # use rowRanges of the input SE for feature coordinates
    feat_gr <- SummarizedExperiment::rowRanges(x)
    idx <- match(rdf[["feature_id"]], names(feat_gr))
    reg_gr <- feat_gr[idx]
    S4Vectors::mcols(reg_gr) <- S4Vectors::DataFrame(
        feature_id = rdf[["feature_id"]],
        pvalue     = rdf[["pvalue"]]
    )

    list(variants = var_gr, regions = reg_gr)
}
