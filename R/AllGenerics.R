setGeneric("runInfo", function(object) standardGeneric("runInfo"))
setGeneric("instrumentInfo", function(object) standardGeneric("instrumentInfo"))
setGeneric("fileName", function(object) standardGeneric("fileName"))
setGeneric("manufacturer", function(object) standardGeneric("manufacturer"))
setGeneric("model", function(object) standardGeneric("model"))
setGeneric("ionisation", function(object) standardGeneric("ionisation"))
setGeneric("analyzer", function(object) standardGeneric("analyzer"))
setGeneric("detector", function(object) standardGeneric("detector"))
setGeneric("isInitialized", function(object) standardGeneric("isInitialized"))
setGeneric("initializeRamp",
           signature=c("object"),
           function(object) standardGeneric("initializeRamp"))
setGeneric("header", function(object,scans,...) standardGeneric("header"))
setGeneric("peaks", function(object,scans,...) standardGeneric("peaks"))
setGeneric("peaksCount", function(object,scans,...) standardGeneric("peaksCount"))
setGeneric("get3Dmap",
           signature=c("object"),
           function(object,scans,lowMz,highMz,resMz,...) standardGeneric("get3Dmap"))
