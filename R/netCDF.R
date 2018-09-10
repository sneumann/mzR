
netCDFIsFile <- function(filename) {
    result <- tryCatch({
        ncid <- nc_open(filename)
        result <- !is.null(ncid)
        netCDFClose(ncid)
        return(result)
    },
    error=function(cond) return(FALSE)
    )
}

netCDFOpen <- function(filename) {
    result <- nc_open(filename, write=FALSE) 
    return(result)
}

netCDFClose <- function(ncid) {
    result <- tryCatch({
##        closedncid <- nc_close(ncid)        
        return(TRUE)
    },
    error=function(cond) return(FALSE)
    )
}

netCDFVarID <- function(ncid, var) {
stop()
    ## result <- .C("NetCDFVarID",
    ##              as.integer(ncid),
    ##              as.character(var),
    ##              id = integer(1),
    ##              status = integer(1),
    ##              PACKAGE = "mzR")
    
    ## if (result$status)
    ##     return(structure(result$status, 
    ##                      errortext = netCDFStrError(result$status)))
    
    ## return(result$id)
}

netCDFVarLen <- function(ncid, var) {

    return(ncid$dim[[var]]$len)
    
}

netCDFVarDouble <- function(ncid, var) {
    as.vector(ncvar_get(ncid, varid=var))
}
netCDFVarInt <- function(ncid, var) {
    as.vector(ncvar_get(ncid, varid=var))
}
netCDFVarText <- function(ncid, var) {
    as.vector(ncvar_get(ncid, varid=var))
}

## netCDFVarDouble <- function(ncid, var) {

##     if (is.character(var))
##         var <- netCDFVarID(ncid, var)
    
##     if (!is.null(attr(var, "errortext")))
##         return(var)
    
##     len <- netCDFVarLen(ncid, var)
##     if (!is.null(attr(len, "errortext")))
##         return(len)
    
##     .C("NetCDFVarDouble",
##        as.integer(ncid),
##        as.integer(var),
##        data = double(len),
##        status = integer(1),
##        PACKAGE = "mzR")$data
## }

## netCDFVarInt <- function(ncid, var) {

##     if (is.character(var))
##         var <- netCDFVarID(ncid, var)
    
##     if (!is.null(attr(var, "errortext")))
##         return(var)
    
##     len <- netCDFVarLen(ncid, var)
##     if (!is.null(attr(len, "errortext")))
##         return(len)
    
##     .C("NetCDFVarInt",
##        as.integer(ncid),
##        as.integer(var),
##        data = integer(len),
##        status = integer(1),
##        PACKAGE = "mzR")$data
## }

## netCDFVarText <- function(ncid, var) {

##     if (is.character(var))
##         var <- netCDFVarID(ncid, var)
    
##     if (!is.null(attr(var, "errortext")))
##         return(var)
        
##     .C("NetCDFVarText",
##        as.integer(ncid),
##        as.integer(var),
##        data = character(1),
##        status = integer(1),
##        PACKAGE = "mzR")$data
## }

netCDFAttText <- function(ncid, att) {
    ncatt_get(ncid, varid=0, attname=att)$value
}


netCDFMSPoints <- function(ncid, scanIndex) {

    if (!is.integer(scanIndex)) scanIndex <- as.integer(scanIndex)

    return (cbind.data.frame(massValues=netCDFVarDouble(ncid, "mass_values"),
                             intensityValues=netCDFVarDouble(ncid, "intensity_values")))
    
    ## var <- netCDFVarID(ncid, "mass_values")
    ## if (!is.null(attr(var, "errortext")))
    ##     return(var)
    
    ## len <- netCDFVarLen(ncid, var)
    ## if (!is.null(attr(len, "errortext")))
    ##     return(len)
    
    ## .C("NetCDFMSPoints",
    ##    as.integer(ncid),
    ##    as.integer(length(scanIndex)),
    ##    scanIndex,
    ##    as.integer(len),
    ##    massValues = double(len),
    ##    intensityValues = double(len),
    ##    status = integer(1),
    ##    PACKAGE = "mzR")[c("massValues", "intensityValues")]
}

netCDFRawData <- function(ncid) {
    
    rt <- netCDFVarDouble(ncid, "scan_acquisition_time")
    tic <- netCDFVarDouble(ncid, "total_intensity")
    scanindex <- netCDFVarInt(ncid, "scan_index")   
    pointValues <- netCDFMSPoints(ncid, scanindex)

    startTimeStamp <- netCDFAttText(ncid, "netcdf_file_date_time_stamp")
    return(list(rt = rt, tic = tic, scanindex = scanindex, 
                mz = pointValues$massValues,
                intensity = pointValues$intensityValues,
                startTimeStamp = startTimeStamp))
}

netCDFRunInfo <- function(ncid) {

    ncraw <- netCDFRawData(ncid)

    return(list(scanCount = length(ncraw$scanindex),
                lowMz = min (ncraw$mz),
                highMz = max (ncraw$mz),
                dStartTime = min (ncraw$rt),
                dEndTime = max (ncraw$rt),
                msLevels = NA,
                startTimeStamp = ncraw$startTimeStamp))
}


netCDFInstrumentInfo <- function(ncid) {
    
    imodel <- netCDFVarText(ncid, "instrument_name")
    imanufacturer <- netCDFVarText(ncid, "instrument_mfr")
    iionisation <- netCDFAttText(ncid, "test_ionization_mode")
    idetector <- netCDFAttText(ncid, "test_detector_type")
    ianalyzer <- NA

    return(list(model = imodel, manufacturer=imanufacturer,
                ionisation = iionisation, detector = idetector,
                analyzer = ianalyzer))
}

