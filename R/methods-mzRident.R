setMethod("mzidInfo",
          signature = "mzRident",
          function(object) return(object@backend$getIDInfo()))

setMethod("psms",
          signature = "mzRident",
          function(object) {
              psms <- object@backend$getPsmInfo()
              specpars <- specParams(object)
              psms <- merge(psms, specpars, by="spectrumID", sort=FALSE)
              psms$acquisitionNum <-
                  as.numeric(sub("^.*=([[:digit:]]+)$", "\\1", psms$spectrumID))
              return(psms)
          })

setMethod("score",
          signature = "mzRident",
          function(x) return(x@backend$getScore()))

setMethod("specParams",
          signature = "mzRident",
          function(object) {
              pars <- object@backend$getSpecParams()
              if ("scan.number.s." %in%  names(pars))
                  if (!is.numeric(pars[, "scan.number.s."]))
                      pars[, "scan.number.s."]  <- as.numeric(as.character(pars[, "scan.number.s."]))
              return(pars)
          })

setMethod("para",
          signature = "mzRident",
          function(object) return(object@backend$getPara()))


setMethod("modifications",
          signature = "mzRident",
          function(object) return(object@backend$getModInfo()))

setMethod("substitutions",
          signature = "mzRident",
          function(object) return(object@backend$getSubInfo()))

setMethod("softwareInfo",
          signature = "mzRident",
          function(object) {
            info <- mzidInfo(object)
            return(info$software)
          })

setMethod("database",
          signature = "mzRident",
          function(object) return(object@backend$getDB()))

setMethod("enzymes",
          signature = "mzRident",
          function(object) {
            info <- mzidInfo(object)
            return(as.data.frame(info$enzymes))
          })

setMethod("sourceInfo",
          signature = "mzRident",
          function(object) {
            info <- mzidInfo(object)
            return(info$SpectraSource)
          })

setMethod("tolerance",
          signature = "mzRident",
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
