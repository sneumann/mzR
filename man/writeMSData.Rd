\name{writeMSData}
\alias{writeMSData}
\alias{writeMSData,list,character-method}

\title{
  Write MS spectrum data to an MS file
}
\usage{

\S4method{writeMSData}{list,character}(object, file, header,
    backend = "pwiz", outformat = "mzml", rtime_seconds = TRUE,
    software_processing)
}
\arguments{
  \item{object}{
    \code{list} containing for each spectrum one \code{matrix} with
    columns \code{mz} (first column) and \code{intensity} (second
    column). See also \code{\link{peaks}} for the method that reads such
    data from an MS file.
  }
  \item{file}{
    \code{character(1)} defining the name of the file.
  }
  \item{header}{
    \code{data.frame} with the header data for the spectra. Has to be in
    the format as the \code{data.frame} returned by the
    \code{\link{header}} method.
  }
  \item{backend}{
    \code{character(1)} defining the backend that should be used for
    writing. Currently only \code{"pwiz"} backend is supported.
  }
  \item{outformat}{
    \code{character(1)} the format of the output file. One of
    \code{"mzml"} or \code{"mzxml"}.
  }
  \item{rtime_seconds}{
    \code{logical(1)} whether the retention time is provided in seconds
    or minutes (defaults to \code{TRUE}).
  }
  \item{software_processing}{
    \code{list} of \code{character} vectors (or single \code{character}
    vector). Each \code{character} vector providing information about
    the software that was used to process the data with optional
    additional description of processing steps. The length of each
    \code{character} vector has to be >= 3: the first element being the
    name of the software, the second string its version and the third
    element the MS CV ID of the software (or \code{"MS:-1"} if not
    known). All additional elements are optional and represent the MS CV
    ID of each processing step performed with the software.
  }
}

\description{
  \code{writeMSData} exports the MS spectrum data provided with
  parameters \code{header} and \code{data} to an MS file in mzML or
  mzXML format.
}

\author{
  Johannes Rainer
}

\seealso{
  \code{\link{copyWriteMSData}} for a function to copy general
  information from a MS data file and writing eventually modified MS
  data from that originating file.
}

\examples{
## Open a MS file and read the spectrum and header information
library(msdata)
fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
    package = "msdata")
ms_fl <- openMSfile(fl, backend = "pwiz")

## Get the spectra
pks <- spectra(ms_fl)
## Get the header
hdr <- header(ms_fl)

## Modify the spectrum data adding 100 to each intensity.
pks <- lapply(pks, function(z) {
    z[, 2] <- z[, 2] + 100
    z
})

## Write the data to a mzML file.
out_file <- tempfile()
writeMSData(object = pks, file = out_file, header = hdr)
}
