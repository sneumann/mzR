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
          
setMethod("database",
          signature="mzRident",
          function(object) {
            info <- mzidInfo(object)           
            return(info$database)
          })
        
setMethod("enzymes",
          signature="mzRident",
          function(object) {
            info <- mzidInfo(object)           
            return(info$enzymes)
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
