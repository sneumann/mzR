\name{openMSfile}
\alias{openMSfile}
\alias{initializeRamp}
\alias{isInitialized}
\alias{fileName}

\title{
  Create and check mzR objects from netCDF, mzXML, mzData or mzML files.
}
\usage{
 openMSfile(filename, backend=c("Ramp", "netCDF"), verbose = FALSE)

 initializeRamp(object)

 isInitialized(object)

 fileName(object)

}
\arguments{
  \item{filename}{ Path name of the netCDF, mzData, mzXML or mzML file to
    read. }
  \item{backend}{ A \code{character} specifiying with backend API to
    use. Currently only 'Ramp' and 'netCDF' are available. 'pwiz' will be
    added in future version. }
  \item{object}{ An instantiated mzR object. }
  \item{verbose}{ Enable verbose output. }
}

\description{
  The \code{openMSfile} constructor will create a new format-specifc
  \code{mzR} object, open 'filename' file and all header information is
  loaded as a Rcpp module and made accessible through the \code{ramp}
  slot of the resulting object. 

  % If the object was created manually with e.g. new("mzR"), then
  % initializeRamp() will open the given file. An mzR object can be
  % queried, whether it is assigned to a raw file with isInitialized(),
  % and finally fileName(object) will return the path to the file "behind"
  % the mzR object.
}

\author{
  Steffen Neumann and Laurent Gatto
}

\examples{
 library(msdata)
 filepath <- system.file("microtofq", package = "msdata")
 file <- list.files(filepath, pattern="MM14.mzML",
                    full.names=TRUE, recursive = TRUE)
  mz <- openMSfile(file)
  fileName(mz)
  isInitialized(mz)
  runInfo(mz)
  close(mz)
}
