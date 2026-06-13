#' Show a QuasarExperiment
#'
#' Displays a compact summary of the object dimensions and key slots.
#'
#' @param object A [QuasarExperiment].
#' @importFrom methods show
setMethod("show", "QuasarExperiment", function(object) {
    cat("class: QuasarExperiment\n")
    cat("features:", nrow(object), " samples:", ncol(object), "\n")
    cat("assays(", length(SummarizedExperiment::assayNames(object)), "):",
        paste(SummarizedExperiment::assayNames(object), collapse = ", "), "\n")
    cat("rowRanges: GRanges with", nrow(object), "features\n")
    cov_cols <- setdiff(colnames(SummarizedExperiment::colData(object)),
                        "fam_index")
    cat("colData(", length(cov_cols), ") covariates:",
        paste(head(cov_cols, 4), collapse = ", "),
        if (length(cov_cols) > 4) "..." else "", "\n")
    cat("geno:", nrow(object@geno), "variants x", ncol(object@geno),
        "samples [BEDMatrix - lazy]\n")
    cat("plinkPrefix:", object@plinkPrefix, "\n")
    cat("grm:", if (is.null(object@grm)) "none" else
        paste0(nrow(object@grm), " x ", ncol(object@grm), " matrix"), "\n")
    invisible(object)
})
