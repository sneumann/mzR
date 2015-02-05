setMethod("mzidInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getIDInfo()))
          
setMethod("psms",
          signature=c("mzRident"),
          function(object) {
              psms <- object@backend$getPsmInfo()
              psms$acquisitionNum <-
                  as.numeric(sub("^.*=([[:digit:]]+)$", "\\1", psms$spectrumID))
              return(psms)
          })
          
setMethod("score",
          signature=c("mzRident"),
          function(object) return(object@backend$getScore()))
          
setMethod("para",
          signature=c("mzRident"),
          function(object) return(object@backend$getPara()))          
          

setMethod("modifications",
          signature=c("mzRident"),
          function(object) return(object@backend$getModInfo()))
          
setMethod("substitutions",
          signature=c("mzRident"),
          function(object) return(object@backend$getSubInfo()))
          
setMethod("softwareInfo",
          signature="mzRident",
          function(object) {
            info <- mzidInfo(object)           
            return(info$software)
          })
          
setMethod("database",
          signature=c("mzRident"),
          function(object) return(object@backend$getDB()))
        
setMethod("enzymes",
          signature="mzRident",
          function(object) {
            info <- mzidInfo(object)           
            return(as.data.frame(info$enzymes))
          })   

setMethod("sourceInfo",
          signature="mzRident",
          function(object) {
            info <- mzidInfo(object)           
            return(info$SpectraSource)
          })  
          
setMethod("tolerance",
          signature="mzRident",
          function(object) {
            info <- mzidInfo(object)   
            ll <- list()
            ll$'FragmentTolerance' <- info$FragmentTolerance
            ll$'ParentTolerance' <- info$ParentTolerance
            return(ll)
          }) 

setMethod("length",
          signature = "mzRident",
          function(x) return(nrow(psms(x))))


setMethod("show",
          signature = "mzRident",
          function(object) {
              filename <- fileName(object)
              cat("Identification file handle.\n")
              cat("Filename: ", basename(filename), "\n")
              cat("Number of psms: ", length(object), "\n")
          })
