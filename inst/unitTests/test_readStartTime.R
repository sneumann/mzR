test_runStartTimeStamp_pwiz <- function() {
    ## Works on (some) mzML, what with mzXML?
    fh <- mzR::openMSfile(MsDataHub::PestMix1_DDA.mzML())
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(!is.na(run_info$startTimeStamp))
    mzR::close(fh)

    fh <- mzR::openMSfile(MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz())
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(!is.na(run_info$startTimeStamp))
    checkTrue(is.character(run_info$startTimeStamp))
    mzR::close(fh)
}

test_runStartTimeStamp_cdf <- function() {
    fh <- mzR::openMSfile(MsDataHub::ko15.CDF())
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(is.character(run_info$startTimeStamp))
    mzR::close(fh)
}


