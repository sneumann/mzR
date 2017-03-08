test_mzXML <- function() {
    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML", package = "msdata")
    mzxml <- openMSfile(file, backend="pwiz")
    checkTrue(class(mzxml)=="mzRpwiz")

    show(mzxml)
    length(mzxml)
    runInfo(mzxml)
    instrumentInfo(mzxml)
    peaks(mzxml)
    peaks(mzxml,1)
    peaks(mzxml, 2:3)
    peaksCount(mzxml)
    header(mzxml)
    header(mzxml,1)
    header(mzxml, 2:3)
    fileName(mzxml)
    close(mzxml) 
}

test_mzML <- function() {
    file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    mzml <- openMSfile(file, backend="pwiz")
    checkTrue(class(mzml)=="mzRpwiz")
    show(mzml)
    length(mzml)
    runInfo(mzml)
    instrumentInfo(mzml)
    peaks(mzml)
    peaks(mzml,1)
    peaks(mzml,2:3)
    peaksCount(mzml)
    header(mzml)
    header(mzml,1)
    header(mzml,2:3)

    checkTrue(ncol(header(mzml))>4)
    checkTrue(length(header(mzml,1))>4)
    checkTrue(ncol(header(mzml,2:3))>4)

    ## Check polarity reporting
    checkTrue(all(header(mzml)$polarity==1))

    fileName(mzml)
    close(mzml)
}
