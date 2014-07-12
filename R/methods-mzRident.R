setMethod("mzidInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getIDInfo()))
          
setMethod("pepInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getPepInfo()))

setMethod("modInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getModInfo()))
          
setMethod("subInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getSubInfo()))
          
setMethod("softwareInfo",
          signature="mzRident",
          function(object) {
            info <- mzidInfo(object)           
            return(info$software)
          })
