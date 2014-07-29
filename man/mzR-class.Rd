\name{mzR-class}
\Rdversion{1.1}
\docType{class}
\alias{mzR-class}
\alias{mzRramp-class}
\alias{mzRnetCDF-class}
\alias{mzRpwiz-class}
\alias{mzRident-class}
\alias{class:mzR}
\alias{class:mzRramp}
\alias{class:mzRnetCDF}
\alias{class:mzRpwiz}
\alias{class:mzRident}
\alias{close,mzRramp-method}
\alias{close,mzRnetCDF-method}
\alias{close}
\alias{analyzer,mzRramp-method}
\alias{analyzer,mzRpwiz-method}
\alias{detector,mzRramp-method}
\alias{detector,mzRpwiz-method}
\alias{fileName,mzR-method}
\alias{get3Dmap,mzRramp-method}
\alias{get3Dmap,mzRpwiz-method}
\alias{header,mzRramp,missing-method}
\alias{header,mzRramp,numeric-method}
\alias{header,mzRpwiz,missing-method}
\alias{header,mzRpwiz,numeric-method}
\alias{header,mzRnetCDF,missing-method}
\alias{header,mzRnetCDF,numeric-method}
\alias{initializeRamp,mzRramp-method}
\alias{instrumentInfo,mzRramp-method}
\alias{instrumentInfo,mzRpwiz-method}
\alias{ionisation,mzRramp-method}
\alias{ionisation,mzRpwiz-method}
\alias{isInitialized,mzRramp-method}
\alias{isInitialized,mzRnetCDF-method}
\alias{length,mzRramp-method}
\alias{length,mzRpwiz-method}
\alias{length,mzRnetCDF-method}
\alias{length}
\alias{manufacturer,mzRramp-method}
\alias{manufacturer,mzRpwiz-method}
\alias{model,mzRramp-method}
\alias{model,mzRpwiz-method}
\alias{peaksCount,mzRramp,missing-method}
\alias{peaksCount,mzRramp,numeric-method}
\alias{peaksCount,mzRpwiz,missing-method}
\alias{peaksCount,mzRpwiz,numeric-method}
\alias{peaks,mzRramp,missing-method}
\alias{peaks,mzRramp,numeric-method}
\alias{peaks,mzRpwiz,missing-method}
\alias{peaks,mzRpwiz,numeric-method}
\alias{peaks,mzRnetCDF,missing-method}
\alias{peaks,mzRnetCDF,numeric-method}
\alias{runInfo,mzRramp-method}
\alias{runInfo,mzRpwiz-method}
\alias{chromatogramsInfo,mzRpwiz-method}
\alias{writeMSfile,mzRpwiz-method}
\alias{softwareInfo,mzRpwiz-method}
\alias{softwareInfo,mzRident-method}
\alias{sampleInfo,mzRpwiz-method}
\alias{sourceInfo,mzRpwiz-method}
\alias{mzidInfo,mzRident-method}
\alias{peptides,mzRident-method}
\alias{modifications,mzRident-method}
\alias{substitutions,mzRident-method}
\alias{database,mzRident-method}
\alias{enzymes,mzRident-method}
\alias{sourceInfo,mzRident-method}
\alias{tolerance,mzRident-method}
\alias{score,mzRident-method}


\title{Class \code{mzR} and sub-classes}
\description{
  The class \code{mzR} is the main class for the common mass spectrometry
  formats. It is a virtual class and thus not supposed to be instanciated
  directly.

  The sub-classes implement specific APIs to access the underlying data
  and metadata in the files. Currently, \code{mzRramp} and \code{mzRpwiz}
  are available implementations. \code{mzRramp} uses the ISB 'RAMP' random access C/C++
  API, and \code{mzRpwiz} uses Proteowizard to access the relevant
  information in \code{mzData}, \code{mzXML}
  and \code{mzML} files. You can also open \code{mz5} file by using \code{mzRpwiz}.

  Additional sub-classes using the proteowizard API and netCDF are
  planned.
}

\section{Objects from the Class}{
  \code{mzR} is a virtual class, so instances cannot be created.

  Objects can be created by calls of the form \code{new("mzRramp",
    ...)}, but more often they will be created with
  \code{\link{openMSfile}}.  
  
  After creating a \code{mzR}, you can write it into a file.
  mzXML, mzML, mgf formats are supported. 
}

\section{Slots}{
  \describe{
    \item{\code{fileName}:}{Object of class \code{character} storing the
      original filename used when the instance was created. }
    \item{\code{backend}:}{ One of the implemented backens or
      \code{NULL}. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versioned"},
      from Biobase.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  Methods currently implemented for \code{mzR}
  \describe{
    \item{fileName}{\code{signature(object = "mzR")}: ... }
  }
  
  Methods currently implemented for \code{mzRramp}
  \describe{
    \item{analyzer}{\code{signature(object = "mzRramp")}: ... }
    \item{close}{\code{signature(con = "mzRramp")}: ... }
    \item{detector}{\code{signature(object = "mzRramp")}: ... }
    \item{fileName}{\code{signature(object = "mzRramp")}: ... }
    \item{get3Dmap}{\code{signature(object = "mzRramp")}: ... }
    \item{header}{\code{signature(object = "mzRramp", scans = "missing")}: ... }
    \item{header}{\code{signature(object = "mzRramp", scans = "numeric")}: ... }
    \item{header}{\code{signature(object = "mzRnetCDF", scans = "missing")}: ... }
    \item{header}{\code{signature(object = "mzRnetCDF", scans = "numeric")}: ... }
    \item{initializeRamp}{\code{signature(object = "mzRramp")}: ... }
    \item{instrumentInfo}{\code{signature(object = "mzRramp")}: ... }
    \item{ionisation}{\code{signature(object = "mzRramp")}: ... }
    \item{isInitialized}{\code{signature(object = "mzRramp")}: ... }
    \item{length}{\code{signature(x = "mzRramp")}: ... }
    \item{manufacturer}{\code{signature(object = "mzRramp")}: ... }
    \item{model}{\code{signature(object = "mzRramp")}: ... }
    \item{peaksCount}{\code{signature(object = "mzRramp", scans = "missing")}: ... }
    \item{peaksCount}{\code{signature(object = "mzRramp", scans = "numeric")}: ... }
    \item{peaks}{\code{signature(object = "mzRramp", scans = "missing")}: ... }
    \item{peaks}{\code{signature(object = "mzRramp", scans = "numeric")}: ... }
    \item{peaks}{\code{signature(object = "mzRnetCDF", scans = "missing")}: ... }
    \item{peaks}{\code{signature(object = "mzRnetCDF", scans = "numeric")}: ... }
    \item{runInfo}{\code{signature(object = "mzRramp")}: ... }

  }
  Methods currently implemented for \code{mzRpwiz}
  \describe{
    \item{analyzer}{\code{signature(object = "mzRpwiz")}: ... }
    \item{detector}{\code{signature(object = "mzRpwiz")}: ... }
    \item{instrumentInfo}{\code{signature(object = "mzRpwiz")}: ... }
    \item{ionisation}{\code{signature(object = "mzRpwiz")}: ... }
    \item{length}{\code{signature(x = "mzRpwiz")}: ... }
    \item{manufacturer}{\code{signature(object = "mzRpwiz")}: ... }
    \item{model}{\code{signature(object = "mzRpwiz")}: ... }
    \item{runInfo}{\code{signature(object = "mzRpwiz")}: ... }
    \item{chromatogramsInfo}{\code{signature(object = "mzRpwiz")}: ... }
    \item{writeMSfile}{\code{signature(object = "mzRpwiz")}: ... }
  }
  Methods currently implemented for \code{mzRident}
  \describe{
    \item{mzidInfo}{\code{signature(object = "mzRident")}: ... }
    \item{peptides}{\code{signature(object = "mzRident")}: ... }
    \item{modifications}{\code{signature(object = "mzRident")}: ... }
    \item{substitutions}{\code{signature(object = "mzRident")}: ... }
    \item{database}{\code{signature(x = "mzRident")}: ... }
    \item{enzymes}{\code{signature(object = "mzRident")}: ... }
    \item{sourceInfo}{\code{signature(object = "mzRident")}: ... }
    \item{tolerance}{\code{signature(object = "mzRident")}: ... }
    \item{score}{\code{signature(object = "mzRident")}: ... }
  }
}

\references{
  RAMP: http://tools.proteomecenter.org/wiki/index.php?title=Software:RAMP
  Proteowizard: http://proteowizard.sourceforge.net/
}

\author{
  Steffen Neumann, Laurent Gatto, Qiang Kou
}

\examples{
 library(msdata)
 filepath <- system.file("microtofq", package = "msdata")
 file <- list.files(filepath, pattern="MM14.mzML",
                    full.names=TRUE, recursive = TRUE)
  mzml <- openMSfile(file)
  close(mzml)
}

\keyword{classes}
