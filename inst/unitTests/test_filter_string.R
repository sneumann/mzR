test_filter_string <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)

    ## mzML - with filter string present.
    fl <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()
    mzml <- openMSfile(fl, backend = "pwiz")
    hdr <- header(mzml)
    mzR::close(mzml)
    checkTrue("filterString" %in% colnames(hdr))
    checkTrue(all(!is.na(hdr$filterString)))
    checkTrue(all(startsWith(hdr$filterString, "FTMS")))
}
