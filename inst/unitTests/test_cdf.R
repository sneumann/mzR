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
  
  close(cdf)
}


test_header <- function() { 
  file <- system.file('cdf/ko15.CDF', package = "msdata")
  cdf <- openMSfile(file, backend="netCDF")        

  h <- header(cdf)
  checkEquals(ncol(h), 21)
  checkEquals(nrow(h), 1278)

  checkTrue(any(colnames(h) == "spectrumId"))
  checkEquals(h$spectrumId, paste0("scan=", h$acquisitionNum))
  
  h <- header(cdf, 1)
  checkEquals(ncol(h), 21)
  checkEquals(nrow(h), 1)

  h <- header(cdf, 2:3)
  checkEquals(ncol(h), 21)
  checkEquals(nrow(h), 2)

  close(cdf)
}
