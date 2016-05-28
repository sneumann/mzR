
#include <Rcpp.h>
#include "pwiz/data/identdata/Version.hpp"

using namespace Rcpp;

RcppExport SEXP mzR_pwiz_version() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(pwiz::identdata::Version::str());
    return __result;
END_RCPP
}
