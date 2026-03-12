test_injection_time <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)
    fl <- MsDataHub::PestMix1_DDA.mzML()
    mzxml <- openMSfile(fl, backend = "pwiz")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(hdr$injectionTime == 0))
    checkTrue(any(colnames(hdr) == "injectionTime"))

    ## CDF
    fl <- MsDataHub::ko15.CDF()
    mzxml <- openMSfile(fl, backend = "netCDF")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(hdr$injectionTime == -1))
    checkTrue(any(colnames(hdr) == "injectionTime"))

    ## mzML - with injection time present.
    fl <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()
    mzxml <- openMSfile(fl, backend = "pwiz")
    hdr <- header(mzxml)
    mzR::close(mzxml)
    checkTrue(all(hdr$injectionTime != 0))
    checkTrue(any(colnames(hdr) == "injectionTime"))    
}
