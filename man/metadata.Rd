\name{metadata}
\alias{analyzer}
\alias{detector}
\alias{instrumentInfo}
\alias{ionisation}
\alias{manufacturer}
\alias{model}
\alias{runInfo}
\alias{chromatogramsInfo}
\alias{softwareInfo}
\alias{sampleInfo}
\alias{sourceInfo}
\alias{mzidInfo}
\alias{modifications}
\alias{psms}
\alias{substitutions}
\alias{database}
\alias{enzymes}
\alias{tolerance}
\alias{score}
\alias{para}
\alias{specParams}

\title{
  Access the metadata from an \code{mzR} object.
}

\usage{
runInfo(object)
chromatogramsInfo(object)
analyzer(object)
detector(object)
instrumentInfo(object)
ionisation(object)
softwareInfo(object)
sampleInfo(object)
sourceInfo(object)
model(object)
mzidInfo(object)
modifications(object, ...)
psms(object, ...)
substitutions(object)
database(object, ...)
enzymes(object)
tolerance(object, ...)
score(x, ...)
para(object)
specParams(object)
}
\arguments{
  \item{object}{An instantiated \code{mzR} object.}
  \item{x}{An instantiated \code{mzR} object.}
  \item{...}{ Additional arguments, currently ignored. }
}

\description{
  Accessors to the analytical setup metadata of a run.
  \code{runInfo} will show a summary of the experiment as a named list,
  including \code{scanCount}, \code{lowMZ}, \code{highMZ},
  \code{dStartTime}, \code{dEndTime} and \code{startTimeStamp}. Note
  that \code{startTimeStamp} can only be extracted from \emph{mzML}
  files using the \emph{pwiz} backend or from \emph{CDF} files. A
  \code{NA} is reported if its value is not available.
  The \code{instrumentInfo} method returns a named \code{list} including
  instrument manufacturer, model, ionisation technique, analyzer and
  detector. \code{mzRpwiz} will give more additional information including
  information on sample, software using and original source file. 
  These individual pieces of information can also be directly
  accessed by the specific methods.
  \code{mzidInfo} is used for the \code{mzR} object created from a mzid file.
  It returns basic information on this mzid file including file provider,
  creation date, software, database, enzymes and spectra data format.
  The \code{mzidInfo} will return the scoring results in identification.
  It will return different results for different searching software used.
}


\seealso{
  See for example \code{\link{peaks}} to access the data for the spectra
  in a \code{"\linkS4class{mzR}"} class.
}

\author{
  Steffen Neumann, Laurent Gatto and Qiang Kou
}

\examples{
 library(msdata)

 file <- system.file("microtofq/MM8.mzML", package = "msdata")
 mz <- openMSfile(file)
 fileName(mz)
 instrumentInfo(mz)
 runInfo(mz)
 close(mz)

 file <- system.file("cdf/ko15.CDF", package = "msdata")
 mz <- openMSfile(file, backend = "netCDF")
 fileName(mz)
 instrumentInfo(mz)
 runInfo(mz)
 close(mz)
 
 file <- system.file("mzid", "Tandem.mzid.gz", package="msdata")
 mzid <- openIDfile(file)
 softwareInfo(mzid)
 enzymes(mzid)
}
