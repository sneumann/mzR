\name{peaks}
\alias{peaks}
\alias{scans}
\alias{header}
\alias{peaksCount}
\alias{get3Dmap}

\title{
  Access the raw data from an \code{mzR} object.
}

\usage{
 header(object, scans, ...)

 peaksCount(object, scans, ...)

 peaks(object, ...) 

 scans(object, ...) ## same as peaks

 get3Dmap(object, scans, lowMz, highMz, resMz, ...)
}

\arguments{
  \item{object}{An instantiated \code{mzR} object.}
  \item{scans}{A \code{numeric} specifying which scans to
    return. Optional for the \code{header}, \code{peaks}, \code{scans} and
    \code{peaksCount} methods. If ommited, the requested data for all
    peaks is returned. } 
  \item{lowMz, highMz}{\code{Numeric}s defining the m/z range to be returned.}
  \item{resMz}{a \code{numeric} defining the m/z resolution.}
  \item{...}{Other arguments. A \code{scan} parameter can be passed to
    \code{peaks}.}
}

\description{

  Access the MS raw data. The \code{peaks}, \code{scans} (can be used
  interchangeably) and \code{peaksCount} functions return the (m/z,
  intensity) pairs and the number peaks in the
  spectrum/spectra. \code{peaks} and \code{scans} return a single matrix
  if \code{scans} is a \code{numeric} of length 1 and a list of matrices
  if several scans are asked for or no \code{scans} argument is provided
  (i.e all spectra in the oject are retured). \code{peaksCount} will
  return a numeric of length \code{n}.
  
  The \code{\link{header}} function returns a list containing
  \code{seqNum}, \code{acquisitionNum} , \code{msLevel},
  \code{peaksCount}, \code{totIonCurrent}, \code{retentionTime},
  \code{basePeakMZ}, \code{basePeakIntensity}, \code{collisionEnergy},
  \code{ionisationEnergy}, \code{lowM}, \code{highMZ},
  \code{precursorScanNum}, \code{precursorMZ}, \code{precursorCharge},
  \code{precursorIntensity}, \code{mergedScan},
  \code{mergedResultScanNum}, \code{mergedResultStartScanNum} and
  \code{mergedResultEndScanNum}, when available in the original file. If
  multiple scans are queried, a \code{data.frame} is returned with the
  scans reported along the rows.

 The \code{get3Dmap} function performs a simple resampling between
 \code{lowMz} and \code{highMz} with \code{reMz} resolution. A matrix of
 dimensions \code{length(scans)} times \code{seq(lowMz,highMz,resMz)} is
 returned. 
}

\seealso{
  \code{\link{instrumentInfo}} for metadata access and the
  \code{"\linkS4class{mzR}"} class.
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
 runInfo(mz)
 colnames(header(mz))
 close(mz)
}
