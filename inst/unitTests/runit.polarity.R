#library(mzR)
#require(msdata) || stop("Cannot find msdata package")
#require(RUnit)

test.polarity <- function()
{

    mzdatapath <- file.path(find.package("msdata"))
    ms1 <- paste(mzdatapath, "/microtofq/MM14.mzdata", sep="")

    id <- mzR:::rampOpen(ms1)
    raw <- mzR:::rampRawData(id)

    checkEquals(length(raw$polarity), 112)

}
    
