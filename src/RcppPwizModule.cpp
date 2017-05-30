#include <Rcpp.h>
#include "RcppPwiz.h"


RCPP_MODULE(Pwiz)
{

    using namespace Rcpp;

    class_<RcppPwiz>( "Pwiz" )
    .constructor("Initialises a new Rccp pwiz object.")
    .method( "open", &RcppPwiz::open, "Opens a mass spec file (mzXML, mzData, etc.) and creates a pwiz object" )
    .method( "close", &RcppPwiz::close, "Close the connection to a mass spec file (mzXML, mzData, etc.)." )
    // .method( "writeMSfile", &RcppPwiz::writeMSfile, "Write a pwiz object into files (mzXML, mzData, etc.)" )
    .method( "writeSpectrumList", &RcppPwiz::writeSpectrumList, "Write the spectrum list to an mzML file." )
    .method( "copyWriteMSfile", &RcppPwiz::copyWriteMSfile, "Copy general content from the originalting MS file, add the provided spectrum list and write the data to a new mzML file." )
    .method( "getFilename", &RcppPwiz::getFilename, "Returns the mass spec filename.")
    .method( "getInstrumentInfo", &RcppPwiz::getInstrumentInfo, "Reads the instrument information from a pwiz object" )
    .method( "getScanHeaderInfo", &RcppPwiz::getScanHeaderInfo, "Reads the header info for one mass spectrum." )
    .method( "getChromatogramsInfo", &RcppPwiz::getChromatogramsInfo, "Reads the chromatogram information.")
    .method( "getPeakList", &RcppPwiz::getPeakList,
             "Performs a non-sequential parsing operation on an indexed mzXML file to obtain the peak list for a numbered scan." )
    .method( "getAllScanHeaderInfo", &RcppPwiz::getAllScanHeaderInfo, "Reads the header info for all mass spectra." )
    .method( "get3DMap", &RcppPwiz::get3DMap, "Reads all scans and returns them as a matrix." )
    .method( "getLastScan", &RcppPwiz::getLastScan, "Returns the last scan (not necessarily the number of scans because of missing scans)." )
    .method( "getLastChrom", &RcppPwiz::getLastChrom, "Returns the index of the last chromatogram." )
    ;
}
