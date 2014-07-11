setMethod("mzidInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getIDInfo()))
          
setMethod("pepInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getPepInfo()))

setMethod("modInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getModInfo()))
