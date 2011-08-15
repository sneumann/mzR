\name{metadata}
\alias{analyzer}
\alias{detector}
\alias{instrumentInfo}
\alias{ionisation}
\alias{manufacturer}
\alias{model}
\alias{runInfo}

\title{
  Access the metadata from an \code{mzR} object.
}

\usage{
runInfo(object)
analyzer(object)
detector(object)
instrumentInfo(object)
ionisation(object)
manufacturer(object)
model(object)
}
\arguments{
  \item{object}{An instantiated \code{mzR} object.}
}

\description{
  Accessors to the analytical setup metadata of a run.
  \code{runInfo} will show a summary of the experiment as a named list,
  including \code{scanCount}, \code{lowMZ}, \code{highMZ}, \code{startMZ},
  \code{endMZ}, \code{dStartTime} and \code{dEndTime}.
  The \code{instrumentInfo} method returns a named \code{list} including
  instrument manufacturer, model, ionisation technique, analyzer and
  detector. These individual pieces of information can also be directly
  accessed by the specific methods.
}


\seealso{
  See for example \code{\link{peaks}} to access the data for the spectra
  in a \code{"\linkS4class{mzR}"} class.
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
  runInfo(mz)
  close(mz)
}
