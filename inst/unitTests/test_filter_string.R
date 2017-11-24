test_filter_string <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)

    ## mzXML
    fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                      package = "msdata")
    mzxml <- openMSfile(fl, backend = "pwiz")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(is.na(hdr$filterString)))
    checkTrue("filterString" %in% colnames(hdr))
    fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                      package = "msdata")
    mzxml <- openMSfile(fl, backend = "Ramp")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(is.na(hdr$filterString)))
    checkTrue("filterString" %in% colnames(hdr))

    ## mzML - with filter string present.
    fl <- system.file("proteomics",
                      "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz",
                      package = "msdata")
    mzml <- openMSfile(fl, backend = "pwiz")
    hdr <- header(mzml)
    mzR::close(mzml)
    checkTrue("filterString" %in% colnames(hdr))
    checkTrue(all(!is.na(hdr$filterString)))
    checkTrue(all(startsWith(hdr$filterString, "FTMS")))
}
