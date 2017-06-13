test_injection_time <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)
    ## mzXML
    fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                      package = "msdata")
    mzxml <- openMSfile(fl, backend = "pwiz")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(hdr$injectionTime == 0))
    checkTrue(any(colnames(hdr) == "injectionTime"))
    fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                      package = "msdata")
    mzxml <- openMSfile(fl, backend = "Ramp")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(hdr$injectionTime == 0))
    checkTrue(any(colnames(hdr) == "injectionTime"))

    ## CDF
    fl <- system.file("cdf", "ko15.CDF",
                      package = "msdata")
    mzxml <- openMSfile(fl, backend = "netCDF")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(hdr$injectionTime == -1))
    checkTrue(any(colnames(hdr) == "injectionTime"))

    ## mzML - with injection time present.
    fl <- system.file("proteomics",
                      "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz",
                      package = "msdata")
    mzxml <- openMSfile(fl, backend = "pwiz")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(hdr$injectionTime != 0))
    checkTrue(any(colnames(hdr) == "injectionTime"))
    
}
