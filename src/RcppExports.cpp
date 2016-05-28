
#include <Rcpp.h>
#include "pwiz/data/msdata/Version.hpp"

using namespace Rcpp;

RcppExport SEXP mzR_pwiz_version() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(pwiz::msdata::Version::str());
    return __result;
END_RCPP
}
