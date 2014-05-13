setMethod("length", 
          signature=c("mzRpwiz"),
          function(x) return(x@backend$getLastScan()))

setMethod("instrumentInfo",
          signature="mzRpwiz",
          function(object) 
          return(object@backend$getInstrumentInfo()))
