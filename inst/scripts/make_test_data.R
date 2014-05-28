library("mzR")
library("msdata")
filepath <- system.file("microtofq", package = "msdata")
file <- list.files(filepath, pattern="MM14.mzML",
                   full.names=TRUE, recursive = TRUE)
mz <- openMSfile(file)
fn <- paste0("microtof_MM14mzXML_",
             gsub(" ", "_", paste(gsub(":", "", date()),
                                  packageVersion("mzR"),
                                  sep = "-")),
             ".rda")

save(mz, file = file.path("../extdata", fn))
close(mz)
