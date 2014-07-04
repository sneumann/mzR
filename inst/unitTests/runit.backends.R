test.backends <- function() {
    library("mzR")
    library("msdata")
    f <- system.file("lockmass", "LockMass_test.mzXML", package = "msdata")
    mr <- openMSfile(f, backend = "ramp")
    checkTrue(validObject(mr))
    mp <- openMSfile(f, backend = "pwiz")
    checkTrue(validObject(mp))
    checkTrue(identical(peaks(mr), peaks(mp)))
    checkTrue(identical(header(mr), header(mp)))
}
