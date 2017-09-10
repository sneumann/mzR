\name{peaks}
\alias{peaks}
\alias{spectra}
\alias{header}
\alias{peaksCount}
\alias{get3Dmap}

\alias{header,mzRramp,missing-method}
\alias{header,mzRramp,numeric-method}
\alias{header,mzRpwiz,missing-method}
\alias{header,mzRpwiz,numeric-method}
\alias{header,mzRnetCDF,missing-method}
\alias{header,mzRnetCDF,numeric-method}

\alias{peaksCount,mzRramp,missing-method}
\alias{peaksCount,mzRramp,numeric-method}
\alias{peaksCount,mzRpwiz,missing-method}
\alias{peaksCount,mzRpwiz,numeric-method}
\alias{peaks,mzRramp-method}
\alias{spectra,mzRramp-method}
\alias{peaks,mzRpwiz-method}
\alias{spectra,mzRpwiz-method}
\alias{peaks,mzRnetCDF-method}
\alias{spectra,mzRnetCDF-method}

\alias{get3Dmap,mzRramp-method}
\alias{get3Dmap,mzRpwiz-method}

\alias{chromatogram,mzRpwiz-method}
\alias{chromatograms,mzRpwiz-method}
\alias{tic,mzRpwiz-method}
\alias{nChrom}
\alias{chromatogram}
\alias{chromatograms}
\alias{tic}


\title{
  Access the raw data from an \code{mzR} object.
}

\usage{

 header(object, scans, ...)

 peaksCount(object, scans, ...)

 peaks(object, ...) 

 spectra(object, ...) ## same as peaks

 get3Dmap(object, scans, lowMz, highMz, resMz, ...)

 chromatogram(object, ...) 

 chromatograms(object, ...) ## same as chromatogram

 tic(object, ...)

 nChrom(object)

}

\arguments{

  \item{object}{An instantiated \code{mzR} object.}
  
  \item{scans}{A \code{numeric} specifying which scans to
    return. Optional for the \code{header}, \code{peaks}, \code{scans}
    and \code{peaksCount} methods. If ommited, the requested data for
    all peaks is returned. }

  \item{lowMz, highMz}{\code{Numeric}s defining the m/z range to be
    returned.}
  
  \item{resMz}{a \code{numeric} defining the m/z resolution.}
  
  \item{...}{Other arguments. A \code{scan} parameter can be passed to
    \code{peaks}.}  }

\description{

  Access the MS raw data. The \code{peaks}, \code{spectra} (can be used
  interchangeably) and \code{peaksCount} functions return the (m/z,
  intensity) pairs and the number peaks in the
  spectrum/spectra. \code{peaks} and \code{spectra} return a single
  matrix if \code{scans} is a \code{numeric} of length 1 and a list of
  matrices if several scans are asked for or no \code{scans} argument is
  provided (i.e all spectra in the oject are retured). \code{peaksCount}
  will return a numeric of length \code{n}.
  
  The \code{\link{header}} function returns a list containing
  \code{seqNum}, \code{acquisitionNum}, \code{msLevel},
  \code{peaksCount}, \code{totIonCurrent}, \code{retentionTime} (in
  seconds), \code{basePeakMZ}, \code{basePeakIntensity},
  \code{collisionEnergy}, \code{ionisationEnergy}, \code{lowM},
  \code{highMZ}, \code{precursorScanNum}, \code{precursorMZ},
  \code{precursorCharge}, \code{precursorIntensity},
  \code{mergedScan}, \code{mergedResultScanNum},
  \code{mergedResultStartScanNum}, \code{mergedResultEndScanNum},
  \code{spectrumId} and
  \code{injectionTime} (ion injection time, in seconds) when available
  in the original file. If multiple scans are queried, a
  \code{data.frame} is returned with the scans reported along the
  rows.

  The \code{get3Dmap} function performs a simple resampling between
  \code{lowMz} and \code{highMz} with \code{reMz} resolution. A matrix
  of dimensions \code{length(scans)} times
  \code{seq(lowMz,highMz,resMz)} is returned.

  The \code{chromatogram} (\code{chromatograms}) accessors return
  chromatograms for the MS file handle. If a single index is provided,
  as \code{data.frame} containing the retention time (1st columns) and
  intensities (2nd column) is returned. The name of the former is always
  \code{time}, while the latter will depend on the run parameters.

  If more than 1 or no chromatogram indices are provided, then a list of
  chromatograms is returned; either those passed as argument or all of
  them. By default, the first (and possibly only) chromatogram is the
  total ion count, which can also be accessed with the \code{tic}
  method.

  The \code{nChrom} function returns the number of chromatograms,
  including the total ion chromatogram.

  Note that access to chromatograms is only supported in the \code{pwiz}
  backend.
  
}

\details{
  The column \code{acquisitionNum} in the \code{data.frame} returned by
  the \code{header} method contains the index during the scan in which
  the signal from the spectrum was measured. The \code{pwiz} backend
  extracts this number from the spectrum's ID provided in the mzML
  file. In contrast, column \code{seqNum} contains the index of each
  spectrum within the file and is thus consecutively numbered. Spectra
  from files with multiple MS levels are linked to each other \emph{via}
  their \code{acquisitionNum}: column \code{precursorScanNum} of an e.g. MS
  level 2 spectrum contains the \code{acquisitionNum} of the related MS
  level 1 spectrum.
}

\note{
  Spectrum identifiers are only specified in \emph{mzML} files, thus,
  for all other file types the column \code{"spectrumId"} of the result
  \code{data.frame} returned by \code{header} contains \code{"scan="}
  followed by the acquisition number of the spectrum. Also, only the
  \code{pwiz} backend supports extraction of the spectras' IDs from
  \emph{mzML} files. Thus, only \emph{mzML} files read with
  \code{backend = "pwiz"} provide the spectrum IDs defined in the file.
}

\seealso{
  \code{\link{instrumentInfo}} for metadata access and the
  \code{"\linkS4class{mzR}"} class.

  \code{\link{writeMSData}} and \code{\link{copyWriteMSData}} for
  functions to write MS data in \emph{mzML} or \code{mzXML} format.
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

 ## A shortgun LCMSMS experiment 
 f <- proteomics(full.names = TRUE,
                 pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
 x <- openMSfile(f, backend = "pwiz")
 x
 nChrom(x)
 head(tic(x))
 head(chromatogram(x, 1L)) ## same as tic(x)
 str(chromatogram(x)) ## as a list

 ## An MRM experiment 
 f <- proteomics(full.names = TRUE, pattern = "MRM")
 x <- openMSfile(f, backend = "pwiz")
 x
 nChrom(x)
 head(tic(x))
 head(chromatogram(x, 1L)) ## same as tic(x)
 str(chromatogram(x, 10:12))  
}
