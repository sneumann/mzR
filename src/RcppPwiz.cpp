#include "RcppPwiz.h"

RcppPwiz::RcppPwiz() {
  msd = NULL;
  runInfo = Rcpp::List::create( );
  isInCacheRunInfo = FALSE;
  instrumentInfo = Rcpp::List::create( );
  isInCacheInstrumentInfo = FALSE;
  allScanHeaderInfo = Rcpp::List::create( );
  isInCacheAllScanHeaderInfo = FALSE;
  filename = "";
}

void RcppPwiz::open(const string& fileName) {

	filename = fileName;
	msd = new MSDataFile(fileName);

}


string RcppPwiz::getFilename (  ) {

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

Rcpp::List RcppPwiz::getScanHeaderInfo ( int whichScan  ) {
	if (msd != NULL) {
		SpectrumListPtr slp = msd->run.spectrumListPtr;
		return slp->size();
		if ((whichScan <= 0) || (whichScan > slp->size())) {
			Rprintf("Index whichScan out of bounds [1 ... %d].\n", slp->size());
			return Rcpp::List::create( );
		}
		
		RAMPAdapter adapter(filename);
		ScanHeaderStruct header;
		adapter.getScanHeader(whichScan, header);

		return Rcpp::List::create(
				Rcpp::_["seqNum"]				= header.seqNum,
			    Rcpp::_["acquisitionNum"]		= header.acquisitionNum,
			    Rcpp::_["msLevel"]				= header.msLevel,
			    Rcpp::_["peaksCount"]     		= header.peaksCount,
			    Rcpp::_["totIonCurrent"]      	= header.totIonCurrent,
			    Rcpp::_["retentionTime"]      	= header.retentionTime,
			    Rcpp::_["basePeakMZ"]      		= header.basePeakMZ,
			    Rcpp::_["basePeakIntensity"]  	= header.basePeakIntensity,
			    Rcpp::_["collisionEnergy"]    	= header.collisionEnergy,
			    Rcpp::_["ionisationEnergy"]  	= header.ionisationEnergy,
			    Rcpp::_["lowMZ"]   				= header.lowMZ,
			    Rcpp::_["highMZ"]      			= header.highMZ,
			    Rcpp::_["precursorScanNum"]   	= header.precursorScanNum,
			    Rcpp::_["precursorMZ"]     		= header.precursorMZ,
			    Rcpp::_["precursorCharge"]      = header.precursorCharge,
			    Rcpp::_["precursorIntensity"] 	= header.precursorIntensity,
			    //Rcpp::_["scanType"]     		= header.scanType,
			    //Rcpp::_["activationMethod"]  	= header.activationMethod,
			    //Rcpp::_["possibleCharges"]   	= header.possibleCharges,
			    //Rcpp::_["numPossibleCharges"]	= header.numPossibleCharges,
			    //Rcpp::_["possibleChargesArray"]	= header.possibleChargesArray,
			    Rcpp::_["mergedScan"]      		= header.mergedScan,
			    Rcpp::_["mergedResultScanNum"]	= header.mergedResultScanNum,
			    Rcpp::_["mergedResultStartScanNum"] = header.mergedResultStartScanNum,
			    Rcpp::_["mergedResultEndScanNum"]   = header.mergedResultEndScanNum
			    //Rcpp::_["filePosition"]      	= header.filePosition
			    );
	}else{
		Rprintf("Warning: pwiz not yet initialized.\n ");
		return Rcpp::List::create( );
	}
}
