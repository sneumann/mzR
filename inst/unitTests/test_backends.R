test_backends <- function() {
    f <- system.file("lockmass", "LockMass_test.mzXML", package = "msdata")
    mr <- openMSfile(f, backend = "Ramp")
    checkTrue(validObject(mr))
    mp <- openMSfile(f, backend = "pwiz")
    checkTrue(validObject(mp))

## Temporarily disabled in 1.99.5 because of SEGV on Windows
##    checkTrue(identical(peaks(mr), peaks(mp)))
##    checkTrue(identical(header(mr), header(mp)))
}

test_mz5 <- function() {
    f <- system.file("microtofq", "MM14.mz5", package = "msdata")
    m <- openMSfile(f)
    checkTrue(validObject(m))
    checkEquals(nrow(header(m)), 112)
}


    

