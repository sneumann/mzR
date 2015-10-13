setGeneric("runInfo", function(object) standardGeneric("runInfo"))
setGeneric("mzidInfo", function(object) standardGeneric("mzidInfo"))
setGeneric("para", function(object) standardGeneric("para"))
setGeneric("substitutions", function(object) standardGeneric("substitutions"))
setGeneric("enzymes", function(object) standardGeneric("enzymes"))
setGeneric("tolerance", function(object) standardGeneric("tolerance"))
setGeneric("instrumentInfo", function(object) standardGeneric("instrumentInfo"))
setGeneric("softwareInfo", function(object) standardGeneric("softwareInfo"))
setGeneric("sampleInfo", function(object) standardGeneric("sampleInfo"))
setGeneric("sourceInfo", function(object) standardGeneric("sourceInfo"))
setGeneric("writeMSfile", function(object, filename, outformat) standardGeneric("writeMSfile"))
setGeneric("chromatogramsInfo", function(object) standardGeneric("chromatogramsInfo"))
setGeneric("creationDate", function(object) standardGeneric("creationDate"))
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
setGeneric("peaksCount", function(object,scans,...) standardGeneric("peaksCount"))
setGeneric("get3Dmap",
           signature=c("object"),
           function(object,scans,lowMz,highMz,resMz,...) standardGeneric("get3Dmap"))


### BiocGenerics
## setGeneric("score", function(x, ...) standardGeneric("score"))
## setGeneric("fileName", function(object, ...) standardGeneric("fileName"))

### ProtGenerics
## setGeneric("psms", function(object, ...) standardGeneric("psms"))
## setGeneric("peaks", function(object, ...) standardGeneric("peaks"))
## setGeneric("database", function(object, ...) standardGeneric("database"))
## setGeneric("modifications", function(object, ...) standardGeneric("modifications"))
