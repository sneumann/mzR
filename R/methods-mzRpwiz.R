setMethod("length", 
          signature=c("mzRpwiz"),
          function(x) return(x@backend$getLastScan()))
