##############################################################
## Defines supported backend APIs
##  - NULL: default
##  - C++Object: for Rcpp modules for ramp and pwiz backends
##  - ncdf4 for netCDF files
setOldClass("ncdf4")
setClassUnion("msAPI",
              c("C++Object","ncdf4", "NULL"))

##############################################################
## mzR main virtual class
## The individual backends are implemented in the different
## sub-classes

setClass("mzR",
         representation(fileName="character",
                        backend="msAPI",
                        "VIRTUAL"),        
         contains=c("Versioned"),
         prototype=prototype(
           fileName = "",
           new("Versioned", versions=c(mzR="0.2.0"))),
         validity=function(object) {
           msg <- validMsg(NULL,NULL)
           if (object@fileName == "")
             msg <- validMsg(msg,"Filename is missing.")
           if (is.null(msg)) TRUE
           else msg
         })


##############################################################
## mzRramp - ramp backend through RcppRamp module
setClass("mzRramp",
         representation(backend="C++Object"),
         contains=c("mzR"),
         prototype=prototype(
           new("Versioned", versions=c(mzR="0.1.0"))),
         validity=function(object) {
           msg <- validMsg(NULL,NULL)
           if (is.null(object@backend))
             msg <- validMsg(msg,"cRamp object not initialised.")
           if (!object@backend$OK())
             msg <- validMsg(msg,"cRamp object not OK.")
           if (object@fileName!=object@backend$getFilename())
             msg <- validMsg(msg,"R slot and ramp filenames do not match.") 
           if (is.null(msg)) TRUE
           else msg
         })


##############################################################
## mzRpwiz - pwiz backend through an Rcpp module 
setClass("mzRpwiz",
         representation(backend="C++Object"),
         contains=c("mzR"),
         prototype=prototype(
           new("Versioned", versions=c(mzR="0.0.1")))
         )


##############################################################
## mzRnetCDF - netCDF backend 
setClass("mzRnetCDF",
         representation(backend="ncdf4"),
         contains=c("mzR"),
         prototype=prototype(
             new("Versioned", versions=c(mzR="0.0.2"))),
         validity=function(object) {
             msg <- validMsg(NULL,NULL)
             if (is.null(object@backend))
                 msg <- validMsg(msg,"ncdf4 object not initialised.")
             if (is.null(object@backend$id))
                 msg <- validMsg(msg,"ncdf4 object is closed.")
             if (is.null(msg)) TRUE
             else msg }
         )

##############################################################
## mzRident - pwiz backend for mzid file
setClass("mzRident",
         representation(backend="C++Object"),
         contains=c("mzR"),
         prototype=prototype(
           new("Versioned", versions=c(mzR="0.0.1")))
         )
