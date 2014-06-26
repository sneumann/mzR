test.backends <- function() {
    library("mzR")
    library("msdata")
    f <- system.file("microtofq", "MSMSpos20_6.mzML", package = "msdata")
    mr <- openMSfile(f, backend = "ramp")
    checkTrue(validObject(mr))
    mp <- openMSfile(f, backend = "pwiz")
    checkTrue(validObject(mp))
    checkTrue(identical(peaks(mr), peaks(mp)))
}
