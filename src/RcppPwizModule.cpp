#include <Rcpp.h>
#include "RcppPwiz.h"


RCPP_MODULE(Pwiz){
	
  using namespace Rcpp;

  class_<RcppPwiz>( "Pwiz" )
    .constructor("Initialises a new Rccp pwiz object.")
    .method( "open", &RcppPwiz::open, "Opens a mass spec file (mzXML, mzData, etc.) and creates a pwiz object" )
    .method( "getFilename", &RcppPwiz::getFilename, "Returns the mass spec filename.")
    .method( "getInstrumentInfo", &RcppPwiz::getInstrumentInfo, "Reads the instrument information from the mzXML header." )
    .method( "getLastScan", &RcppPwiz::getLastScan, "Returns the last scan (not necessarily the number of scans because of missing scans)." )
    ;
}
