#include <Rcpp.h>
#include "RcppRamp.h"


RCPP_MODULE(Ramp){
  using namespace Rcpp;

  class_<RcppRamp>( "Ramp" )
    .constructor("Initialises a new Rccp Ramp object.")
    .method( "open", &RcppRamp::open, "Opens a mass spec file (mzXML, mzData, etc.) and creates a cRamp object" )
    .method( "close", &RcppRamp::close, "Closes the mzXML file. Releases the memory of the cRamp object." )
    .method( "getFilename", &RcppRamp::getFilename, "Returns the mass spec filename.")
    .method( "getRunInfo", &RcppRamp::getRunInfo, "Reads the run information from the mzXML header." )
    .method( "getInstrumentInfo", &RcppRamp::getInstrumentInfo, "Reads the instrument information from the mzXML header." )
    .method( "getScanHeaderInfo", &RcppRamp::getScanHeaderInfo, "Reads the header info for one mass spectrum." )
    .method( "getAllScanHeaderInfo", &RcppRamp::getAllScanHeaderInfo, "Reads the header info for all mass spectra." )
    .method( "getPeakList", &RcppRamp::getPeakList, 
	     "Performs a non-sequential parsing operation on an indexed mzXML file to obtain the peak list for a numbered scan." )
    .method( "get3DMap", &RcppRamp::get3DMap, "Reads al scans and returns them as a matrix." )
    .method( "getLastScan", &RcppRamp::getLastScan, "Returns the last scan (not necessarily the number of scans because of missing scans)." )
    .method( "OK", &RcppRamp::OK, "Checks the status of the object." )
    ;
}

