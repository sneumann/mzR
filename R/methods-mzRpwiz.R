setMethod("length", 
          signature=c("mzRpwiz"),
          function(x) return(x@backend$getLastScan()))

setMethod("instrumentInfo",
          signature="mzRpwiz",
          function(object) 
          return(object@backend$getInstrumentInfo()))

setMethod("manufacturer",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$manufacturer)
          })

setMethod("model",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$model)
          })

setMethod("ionisation",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$ionisation)
          })

setMethod("analyzer",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$analyzer)
          })

setMethod("detector",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$detector)
          })
