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
\alias{peptides}
\alias{substitutions}
\alias{database}
\alias{enzymes}
\alias{tolerance}
\alias{score}

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
modifications(object)
peptides(object)
substitutions(object)
database(object)
enzymes(object)
tolerance(object)
score(object)
}
\arguments{
  \item{object}{An instantiated \code{mzR} object.}
}

\description{
  Accessors to the analytical setup metadata of a run.
  \code{runInfo} will show a summary of the experiment as a named list,
  including \code{scanCount}, \code{lowMZ}, \code{highMZ},
  \code{dStartTime} and \code{dEndTime}.
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
 filepath <- system.file("microtofq", package = "msdata")
 file <- list.files(filepath, pattern="MM14.mzML",
                    full.names=TRUE, recursive = TRUE)
  mz <- openMSfile(file, backend = "Ramp")
  fileName(mz)
  instrumentInfo(mz)
}
