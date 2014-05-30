BUILT_RCPP_VERSION = package_version("0.11.1")

.onLoad <-
    function(libname, pkgname) {
      ## checkeing installed vs build-time Rcpp version
      if (Sys.info()['sysname'] %in% c("Darwin", "Windows")) {
        installedRcpp <- utils::packageVersion("Rcpp")
          if (installedRcpp != BUILT_RCPP_VERSION) { # use > instead of !=?            
            msg <- paste0("mzR has been built against a different Rcpp version (", BUILT_RCPP_VERSION, ")\n",
                          "than is installed on your system (", installedRcpp, "). This might lead to errors\n",
                          "when loading mzR. If you encounter such issues, please send a report,\n",
                          "including the output of sessionInfo() to the Bioc mailing list at \n",
                          "http://www.bioconductor.org/help/mailing-list. For details see also\n",
                          "https://github.com/sneumann/mzR/wiki/mzR-Rcpp-compiler-linker-issue.")
            warning(msg)            
          }
      }
      require2 <- require
      require2("methods", character.only = TRUE, quietly = TRUE)
      loadRcppModules()

    }
