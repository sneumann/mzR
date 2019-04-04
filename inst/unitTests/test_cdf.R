test_open <- function() {
    file <- system.file('cdf/ko15.CDF', package = "msdata")
    cdf <- openMSfile(file, backend="netCDF")    
    checkTrue(class(cdf)=="mzRnetCDF")
    close(cdf)
  }

test_length <- function() {
  file <- system.file('cdf/ko15.CDF', package = "msdata")
  cdf <- openMSfile(file, backend="netCDF")     
  
  checkEquals(length(cdf),1278)
  
  close(cdf)
}

test_peaks <- function() {
  file <- system.file('cdf/ko15.CDF', package = "msdata")
  cdf <- openMSfile(file, backend="netCDF")       

  checkTrue(class(cdf)=="mzRnetCDF")

  p <- peaks(cdf)
  checkTrue(length(p)==1278)
  checkTrue(all(colnames(p)==c("mz", "intensity")))
  
  p <- peaks(cdf,1)
  checkTrue(ncol(p)==2)
  checkTrue(nrow(p)==429)
  
  p <- peaks(cdf,2:3)
  checkTrue(length(p)==2)

  p <- peaks(cdf,1278)
  checkTrue(ncol(p)==2)
  checkTrue(nrow(p)==40)

  ## Can we check that this indeed throws an error:
  ## Error in seq.default(rawdata$scanindex[scans] + 1, min(rawdata$scanindex[scans +  : 'from' must be a finite number
  ## p <- peaks(cdf,1279)


  ri1 <- list(scanCount = 1278L, lowMz = 200, highMz = 600,
              dStartTime = 2501.378, dEndTime = 4499.824,
              msLevels = NA, startTimeStamp = "2004,06,01,10:28:03+0800")
  ri2 <- runInfo(cdf)
  checkTrue(all.equal(ri1, ri2))
  
  ii1 <- list(model = "                               ",
              manufacturer = "Agilent Technologies",
              ionisation = "Electrospray Ionization",
              detector = "Conversion Dynode Electron Multiplier",
              analyzer = NA)
  ii2 <- instrumentInfo(cdf)

  checkTrue(all.equal(ii1, ii2))
  
  close(cdf)
}


test_header <- function() { 
  file <- system.file('cdf/ko15.CDF', package = "msdata")
  cdf <- openMSfile(file, backend="netCDF")        

  h <- header(cdf)
  checkEquals(ncol(h), 26)
  checkEquals(nrow(h), 1278)
  checkTrue(any(colnames(h) == "centroided"))
  checkTrue(all(is.na(h$centroided)))
  checkTrue(any(colnames(h) == "ionMobilityDriftTime"))
  checkTrue(all(h$ionMobilityDriftTime == -1))
  
  checkTrue(any(colnames(h) == "spectrumId"))
  checkEquals(h$spectrumId, paste0("scan=", h$acquisitionNum))
  
  h <- header(cdf, 1)
  checkEquals(ncol(h), 26)
  checkEquals(nrow(h), 1)

  h <- header(cdf, 2:3)
  checkEquals(ncol(h), 26)
  checkEquals(nrow(h), 2)

  close(cdf)
}

test_chromatogram <- function() {
    file <- system.file('cdf/ko15.CDF', package = "msdata")
    x <- openMSfile(file, backend="netCDF")        
    suppressWarnings(
        chr <- chromatogram(x)
    )
    checkTrue(length(chr) == 0)
    suppressWarnings(
        chr <- chromatograms(x)
    )
    checkTrue(length(chr) == 0)
    close(x)
}

test_chromatogramHeader <- function() {
    file <- system.file('cdf/ko15.CDF', package = "msdata")
    x <- openMSfile(file, backend="netCDF")        
    suppressWarnings(
        ch <- chromatogramHeader(x)
    )
    checkTrue(nrow(ch) == 0)
    close(x)
}


if (FALSE) {
    suite <- RUnit::defineTestSuite(name = paste("mzR", "RUnit Tests"), 
                                    dirs = "/vol/R/BioC/devel/mzR/inst/unitTests",
                                    testFileRegexp = "^test_.*\\.R$", rngKind = "default", 
                                    rngNormalKind = "default")

    result <- RUnit::runTestSuite(suite)
    RUnit::printTextProtocol(result, showDetails = FALSE)

}
