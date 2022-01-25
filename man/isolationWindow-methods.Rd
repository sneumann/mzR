\name{isolationWindow-methods}
\docType{methods}
\alias{isolationWindow}
\alias{isolationWindow-methods}
\alias{isolationWindow,character-method}
\alias{isolationWindow,mzRpwiz-method}

\title{ Returns the ion selection isolation window }

\description{

  The methods return matrices of lower (column \code{low}) and upper
  (column \code{high}) isolation window offsets. Matrices are returned
  as a list of length equal to the number of input files (provided as
  file names of raw mass spectrometry data objects, see below). By
  default (i.e when \code{unique. = TRUE} ), only unique offsets are
  returned, as they are expected to identical for all spectra per
  acquisition. If this is not the case, a message is displayed.
  
}


\section{Methods}{
  \describe{

    \item{\code{signature(object = "character", unique. = "logical",
	simplify = "logical")}}{Returns the isolation window for the
      file \code{object}. By default, only unique isolation windows
      are returned per file (\code{unique = TRUE}); if set to
      \code{FALSE}, a \code{matrix} with as many rows as there are MS2
      spectra. If only one file passed an input and \code{simplify} is
      set to \code{TRUE} (default), the resulting \code{list} of
      length 1 is simplified to a \code{matrix}.  }

    \item{\code{signature(object = "mzRpwiz", unique. = "logical",
	simplify = "logical")}}{As above for \code{mzRpwiz} objects. }

  }
}

\author{

  Laurent Gatto <lg390@cam.ac.uk> based on the functionality from the
  \code{msPurity:::get_isolation_offsets} function.

}

\examples{
library("msdata")
f <- msdata::proteomics(full.names = TRUE,
                        pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
isolationWindow(f)

rw <- openMSfile(f)
isolationWindow(rw)
str(isolationWindow(rw, unique = FALSE))
}

\keyword{methods}