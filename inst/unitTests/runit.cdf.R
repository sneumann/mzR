test.open <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    cdf <- openMSfile(file)    
    checkTrue(class(cdf)=="mzRnetCDF")
    close(cdf)
  }

test.length <- function() {
  file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
  cdf <- openMSfile(file)    
  
  checkEquals(length(cdf),1278)
  
  close(cdf)
}

test.peaks <- function() {
  file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
  cdf <- openMSfile(file)    

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


test.header <- function() {
  file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
  cdf <- openMSfile(file)    

  h <- header(cdf)
  checkEquals(ncol(h), 19)
  checkEquals(nrow(h), 1278)

  h <- header(cdf, 1)
  checkEquals(ncol(h), 19)
  checkEquals(nrow(h), 1)

  h <- header(cdf, 2:3)
  checkEquals(ncol(h), 19)
  checkEquals(nrow(h), 2)

  close(cdf)
}
