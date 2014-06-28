setMethod("idInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getIDInfo()))
