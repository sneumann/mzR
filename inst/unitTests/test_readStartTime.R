test_runStartTimeStamp_pwiz <- function() {
    ## Works on (some) mzML, what with mzXML?
    library(msdata)
    library(mzR)
    f <- system.file("lockmass/LockMass_test.mzXML", package = "msdata")
    fh <- mzR::openMSfile(f)
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(is.na(run_info$startTimeStamp))
    mzR::close(fh)
    
    f <- system.file("microtofq/MM14.mzML", package = "msdata")
    fh <- mzR::openMSfile(f)
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(is.na(run_info$startTimeStamp))
    mzR::close(fh)
    
    f <- system.file("microtofq/MM8.mzML", package = "msdata")
    fh <- mzR::openMSfile(f)
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(!is.na(run_info$startTimeStamp))
    checkTrue(is.character(run_info$startTimeStamp))
    mzR::close(fh)
    
    f <- system.file("proteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz", package = "msdata")
    fh <- mzR::openMSfile(f)
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(!is.na(run_info$startTimeStamp))
    checkTrue(is.character(run_info$startTimeStamp))
    mzR::close(fh)
}

test_runStartTimeStamp_cdf <- function() {
    ## Can not extract from CDF
    f <- system.file("cdf/ko15.CDF", package = "msdata")
    fh <- mzR::openMSfile(f, backend = "netCDF")
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(is.character(run_info$startTimeStamp))
    mzR::close(fh)
}

test_runStartTimeStamp_ramp <- function() {
    ## Can not extract from ramp
    f <- system.file("iontrap/extracted.mzData", package = "msdata")
    fh <- mzR::openMSfile(f, backend = "Ramp")
    run_info <- runInfo(fh)
    checkTrue(any(names(run_info) == "startTimeStamp"))
    checkTrue(is.na(run_info$startTimeStamp))
    mzR::close(fh)
}


