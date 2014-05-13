#include "RcppPwiz.h"

RcppPwiz::RcppPwiz() {
  msd = NULL;
  runInfo = Rcpp::List::create( );
  isInCacheRunInfo = FALSE;
  instrumentInfo = Rcpp::List::create( );
  isInCacheInstrumentInfo = FALSE;
  allScanHeaderInfo = Rcpp::List::create( );
  isInCacheAllScanHeaderInfo = FALSE;
  filename = Rcpp::StringVector::create( );
}

void RcppPwiz::open(const string& fileName) {

	filename = Rcpp::StringVector::create( fileName );
	msd = new MSDataFile(fileName);

}


Rcpp::StringVector RcppPwiz::getFilename (  ) {

  return filename;
}

int RcppPwiz::getLastScan() const {
  if (msd != NULL) {
	  SpectrumListPtr slp = msd->run.spectrumListPtr;
	  return slp->size();
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return -1;
}

Rcpp::List RcppPwiz::getInstrumentInfo ( ) {
  if (msd != NULL) {
    if (!isInCacheInstrumentInfo) {

      vector<InstrumentConfigurationPtr> icp = msd->instrumentConfigurationPtrs; // NULL for mzData

      if (icp.size() != 0) { 
      
      CVTranslator cvTranslator;
      LegacyAdapter_Instrument adapter(*icp[0], cvTranslator);
      
      instrumentInfo = Rcpp::List::create(
					  Rcpp::_["manufacturer"]  = std::string(adapter.manufacturer()),
					  Rcpp::_["model"]         = std::string(adapter.model()),
					  Rcpp::_["ionisation"]    = std::string(adapter.ionisation()),
					  Rcpp::_["analyzer"]      = std::string(adapter.analyzer()),
					  Rcpp::_["detector"]      = std::string(adapter.detector() )
					  ) ;

      } else {
      instrumentInfo = Rcpp::List::create(
					  Rcpp::_["manufacturer"]  = "",
					  Rcpp::_["model"]         = "",
					  Rcpp::_["ionisation"]    = "",
					  Rcpp::_["analyzer"]      = "",
					  Rcpp::_["detector"]      = ""
					  ) ;
      }
      isInCacheInstrumentInfo = TRUE;
    } 
    return(instrumentInfo);
  }
  Rprintf("Warning: pwiz not yet initialized.\n ");
  return instrumentInfo;
}
