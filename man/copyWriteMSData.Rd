\name{copyWriteMSData}
\alias{copyWriteMSData}

\title{
  Write MS spectrum data to a MS file copying metadata from the
  originating file
}

\usage{
  copyWriteMSData(object, file, original_file, header, backend =
    "pwiz", outformat = "mzml", rtime_seconds = TRUE, software_processing)
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

  \item{original_file}{
    \code{character(1)} with the name of the original file from which
    the spectrum data was first read.
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
  Copy general information from the originating MS file and write this,
  along with the provided spectra data, to a new file. The expected
  workflow is the following: data is first loaded from an MS file,
  e.g. using \code{\link{peaks}} and \code{\link{header}} methods,
  processed in R and then saved again to an MS file providing the
  (eventually) manipulated spectra and header data with arguments
  \code{header} and \code{data}.
}

\note{
  \code{copyWriteMSData} supports at present copying data from
  \code{mzXML} and \code{mzML} and exporting to \code{mzML}. Copying and
  exporting to \code{mzXML} can fail for some input files.

  The intention of this function is to copy data from an existing file
  and save it along with eventually modified data to a new file. To
  write new MS data files use the \code{\link{writeMSData}} function
  instead.
}

\author{
  Johannes Rainer
}

\seealso{
  \code{\link{writeMSData}} for a function to save MS data to a new mzML
  or mzXML file.
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

## Copy metadata and additional information from the originating file
## and save it, along with the modified data, to a new mzML file.
out_file <- tempfile()
copyWriteMSData(pks, file = out_file, original_file = fl,
    header = hdr)
}
