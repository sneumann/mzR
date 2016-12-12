test_oldramp.mzXML <- function() {
    cdfpath <- system.file("threonine", package = "msdata")
    filename <- list.files(cdfpath, pattern="threonine_i2_e35_pH_tree.mzXML",
                       full.names=TRUE, recursive = TRUE)

    rampid <- mzR:::rampOpen(filename)
    if (rampid < 0)
       stop("Could not open mzXML/mzData file")

    on.exit(mzR:::rampClose(rampid))
    rawdata <- mzR:::rampRawData(rampid)
    mzR:::rampClose(rampid)
}

test_oldramp.mzML <- function() {
    cdfpath <- system.file("microtofq", package = "msdata")
    filename <- list.files(cdfpath, pattern="MM14.mzML",
                       full.names=TRUE, recursive = TRUE)

    rampid <- mzR:::rampOpen(filename)
    if (rampid < 0)
       stop("Could not open mzXML/mzData file")

    on.exit(mzR:::rampClose(rampid))
    rawdata <- mzR:::rampRawData(rampid)
    mzR:::rampClose(rampid)
}

test_oldramp.mzData <- function() {
    cdfpath <- system.file("microtofq", package = "msdata")
    filename <- list.files(cdfpath, pattern="MM14.mzdata$",
                       full.names=TRUE, recursive = TRUE)

    rampid <- mzR:::rampOpen(filename)
    if (rampid < 0)
       stop("Could not open mzXML/mzData file")

    on.exit(mzR:::rampClose(rampid))
    rawdata <- mzR:::rampRawData(rampid)
    mzR:::rampClose(rampid)
}

test_oldramp.mzData.gz <- function() {
    cdfpath <- system.file("microtofq", package = "msdata")
    filename <- list.files(cdfpath, pattern = "MM14.mzdata.gz",
                       full.names = TRUE, recursive = TRUE)

    rampid <- mzR:::rampOpen(filename)
    if (rampid < 0)
       stop("Could not open mzXML/mzData file")

    on.exit(mzR:::rampClose(rampid))
    rawdata <- mzR:::rampRawData(rampid)
    mzR:::rampClose(rampid)
}
