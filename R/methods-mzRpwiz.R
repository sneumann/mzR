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

setMethod("header",
          signature=c("mzRpwiz","missing"),
          function(object) return(object@backend$getAllScanHeaderInfo()))

setMethod("header",
          signature=c("mzRpwiz","numeric"),
          function(object, scans) {
            if (length(scans)==1) {
              return(object@backend$getScanHeaderInfo(scans))
            } else {
              return(data.frame(t(sapply(scans,function(x) unlist(object@backend$getScanHeaderInfo(x))))))
            }
          })
