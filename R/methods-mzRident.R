setMethod("dateInfo",
          signature=c("mzRident"),
          function(object) return(object@backend$getCreationDate()))
