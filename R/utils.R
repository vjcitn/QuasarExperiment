# Match phenotype sample IDs against .fam IIDs.
#
# Tries direct match first. If that fails completely, tries stripping a
# leading run of digits + underscore from the fam IIDs (handles pipelines
# that write the FID into the IID column, e.g. "0_HG00096" -> "HG00096").
# Stops with an informative message if neither attempt yields a complete
# match.
.match_fam_ids <- function(sample_ids, fam_iids) {
    idx <- match(sample_ids, fam_iids)
    if (!anyNA(idx)) return(idx)

    # try stripping leading numeric prefix from fam IIDs
    stripped <- sub("^[0-9]+_", "", fam_iids)
    if (!identical(stripped, fam_iids)) {
        idx2 <- match(sample_ids, stripped)
        if (!anyNA(idx2)) {
            message("Matched samples after stripping leading numeric prefix ",
                    "from .fam IIDs (e.g. '0_ID' -> 'ID')")
            return(idx2)
        }
    }

    # neither worked — report which IDs failed
    unmatched <- sample_ids[is.na(idx)]
    stop(
        length(unmatched), " sample(s) in phenotype file not found in .fam file.\n",
        "First unmatched: ", paste(head(unmatched, 5), collapse = ", "),
        if (length(unmatched) > 5) " ..." else "", "\n",
        "First .fam IIDs: ", paste(head(fam_iids, 5), collapse = ", "),
        if (length(fam_iids) > 5) " ..." else "", "\n",
        "Check that sample IDs in the phenotype file match the IID column ",
        "(column 2) of the .fam file."
    )
}
