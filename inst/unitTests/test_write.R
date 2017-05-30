## Testing to write stuff.

dontrun_write <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)
    test_folder = "/Users/jo/Desktop/"
    ## file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    ## mzml <- openMSfile(file, backend="pwiz")
    ## mzR:::writeMSfile(mzml, filename = "/Users/jo/Desktop/test.mzML",
    ##                   outformat = "mzml")
    ## mzR::close(mzml)
    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                        package = "msdata")
    mzxml <- openMSfile(file, backend = "pwiz")
    pks <- peaks(mzxml)
    hdr <- header(mzxml)
    mzR::close(mzxml)

    fnew <- paste0(test_folder, "test_copyWrite.mzML")
    mzR:::copyWriteMSfile(filename = fnew, originalFile = file, header = hdr,
                          data = pks)
    ##mzR:::replaceSpectrumList(mzml, hdr, pks)
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    mzR::close(mzml_new)

    ## Check that content is the same
    checkEquals(pks_new, pks)
    checkEquals(hdr_new, hdr)  ## polarity is OK here
        
    ## Now, let's pick selected spectra instead.
    hdr_sub <- hdr[c(1, 3, 5), ]
    pks_sub <- pks[c(1, 3, 5)]
    ## Using MSnbase...
}

dontrun_test_writeSpectrumList <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)
    test_folder = "/Users/jo/Desktop/"
    test_file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                             package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)

    ## Test writing the data.
    out_file <- paste0(test_folder, "test_write.mzML")
    mzR:::writeSpectrumList(out_file, header = hdr, data = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)  ## polarity does not match!
    checkEquals(pks, pks_2)

    ## Test subsetting.
    hdr_sub <- hdr[c(1, 3, 5), ]
    hdr_sub$seqNum <- 1:nrow(hdr_sub)
    pks_sub <- pks[c(1, 3, 5)]
    mzR:::writeSpectrumList(out_file, header = hdr_sub, data = pks_sub)
    in_file <- openMSfile(out_file)
    hdr_sub_2 <- header(in_file)
    pks_sub_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
    
    ## Test writing as other file type.
}
