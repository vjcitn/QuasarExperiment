exdir <- system.file("extdata", package = "QuasarExperiment")

skip_if_no_extdata <- function() {
    skip_if(
        !file.exists(file.path(exdir, "mean-pheno-n100.bed")),
        "Example data not installed (run inst/scripts/install_extdata.R)"
    )
}

test_that("QuasarExperiment constructs without GRM", {
    skip_if_no_extdata()
    qe <- QuasarExperiment(
        plinkPrefix = file.path(exdir, "chr22-n100"),
        phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
        covFile     = file.path(exdir, "cov-n100.tsv")
    )
    expect_s4_class(qe, "QuasarExperiment")
    expect_equal(nrow(qe), 20L)
    expect_equal(ncol(qe), 100L)
    expect_null(grm(qe))
})

test_that("QuasarExperiment constructs with GRM", {
    skip_if_no_extdata()
    qe <- QuasarExperiment(
        plinkPrefix = file.path(exdir, "chr22-n100"),
        phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
        covFile     = file.path(exdir, "cov-n100.tsv"),
        grmFile     = file.path(exdir, "grm-n100.tsv")
    )
    expect_false(is.null(grm(qe)))
    expect_equal(dim(grm(qe)), c(100L, 100L))
})

test_that("accessors round-trip", {
    skip_if_no_extdata()
    qe <- QuasarExperiment(
        plinkPrefix = file.path(exdir, "chr22-n100"),
        phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
        covFile     = file.path(exdir, "cov-n100.tsv")
    )
    expect_equal(nrow(geno(qe)), length(variantRanges(qe)))
    expect_equal(ncol(geno(qe)), ncol(qe))
    expect_equal(plinkPrefix(qe), file.path(exdir, "chr22-n100"))
})

test_that("runQuasar lm returns GRanges", {
    skip_if_no_extdata()
    bin <- Sys.which("quasar")
    skip_if(bin == "", "quasar binary not on PATH")
    qe <- QuasarExperiment(
        plinkPrefix = file.path(exdir, "chr22-n100"),
        phenoFile   = file.path(exdir, "mean-pheno-n100.bed"),
        covFile     = file.path(exdir, "cov-n100.tsv")
    )
    res <- runQuasar(qe, binary = bin, model = "lm", mode = "cis")
    expect_s4_class(variantHits(res), "GRanges")
    expect_s4_class(regionHits(res), "GRanges")
    expect_equal(length(regionHits(res)), nrow(qe))
    expect_true("pvalue" %in% names(S4Vectors::mcols(variantHits(res))))
})
