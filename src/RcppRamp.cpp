#include "RcppRamp.h"

#include "ListBuilder.h"


RcppRamp::RcppRamp() {
  ramp = NULL;
  runInfo = Rcpp::List::create( );
  isInCacheRunInfo = FALSE;
  instrumentInfo = Rcpp::List::create( );
  isInCacheInstrumentInfo = FALSE;
  allScanHeaderInfo = Rcpp::DataFrame::create( );
  isInCacheAllScanHeaderInfo = FALSE;
  filename = Rcpp::StringVector::create( );
}

RcppRamp::~RcppRamp() {
  RcppRamp::close();
}

void
RcppRamp::open( const char* fileName, bool declaredScansOnly ) {
  RcppRamp::close();
  ramp = new cRamp(fileName, declaredScansOnly);
  if (ramp->OK()) {
    filename = Rcpp::StringVector::create( fileName );
  } else {
    RcppRamp::close();
    Rprintf("Failed to open file.\n ");
  }
}

void
RcppRamp::close() {
  if (ramp != NULL) {
    delete ramp;
    ramp = NULL;
    runInfo = Rcpp::List::create( );
    isInCacheRunInfo = FALSE;
    instrumentInfo = Rcpp::List::create( );
    isInCacheInstrumentInfo = FALSE;
    allScanHeaderInfo = Rcpp::List::create( );
    isInCacheAllScanHeaderInfo = FALSE;
    filename = Rcpp::StringVector::create( );
  }
}


Rcpp::StringVector
RcppRamp::getFilename (  ) {
  if (ramp != NULL) {
    return filename;
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return filename;
}

Rcpp::List
RcppRamp::getRunInfo (  ) {
  if (ramp != NULL) {
    if (!isInCacheRunInfo) {
      // Rprintf("Read from disk.\n ");
      rampRunInfo *info = ramp->getRunInfo();
      RunHeaderStruct data = info->m_data;
      delete info;
      runInfo = Rcpp::List::create(
				   Rcpp::_["scanCount"]  = data.scanCount,
				   Rcpp::_["lowMZ"]      = data.lowMZ,
				   Rcpp::_["highMZ"]     = data.highMZ,
				   Rcpp::_["startMZ"]    = data.startMZ,
				   Rcpp::_["endMZ"]      = data.endMZ,
				   Rcpp::_["dStartTime"] = data.dStartTime,
				   Rcpp::_["dEndTime"]   = data.dEndTime
				   );
      isInCacheRunInfo = TRUE;
    } else {
      // Rprintf("Read from cache.\n ");
    }
    return runInfo;
  }
  Rprintf("Warning: Ramp not yet initialized.\n");
  return runInfo;
}

Rcpp::List
RcppRamp::getInstrumentInfo ( ) {
  if (ramp != NULL) {
    if (!isInCacheInstrumentInfo) {
      // Rprintf("Read from disk.\n ");
      rampInstrumentInfo *info = ramp->getInstrumentInfo(); // NULL for mzData

      if (info != NULL) { 
      InstrumentStruct * data = info->m_instrumentStructPtr;
      
      instrumentInfo = Rcpp::List::create(
					  Rcpp::_["manufacturer"]  = std::string(data->manufacturer),
					  Rcpp::_["model"]         = std::string(data->model),
					  Rcpp::_["ionisation"]    = std::string(data->ionisation),
					  Rcpp::_["analyzer"]      = std::string(data->analyzer),
					  Rcpp::_["detector"]      = std::string(data->detector)
					  ) ;
      delete info; 
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
    } else {
      // Rprintf("Read from cache.\n ");
    }
    return(instrumentInfo);
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return instrumentInfo;
}

Rcpp::List
RcppRamp::getScanHeaderInfo ( int whichScan  ) {
  if (ramp != NULL) {
    if ((whichScan <= 0) || (whichScan > ramp->getLastScan())) {
      Rprintf("Index whichScan out of bounds [1 ... %d].\n", ramp->getLastScan());
      return Rcpp::List::create( );
    }
    rampScanInfo *info = ramp->getScanHeaderInfo( whichScan );
    ScanHeaderStruct data = info->m_data;
    delete info;

    //    Rcpp::List header = Rcpp::List::create();
    ListBuilder header;

    header.add("seqNum",  Rcpp::wrap(data.seqNum));
    header.add("acquisitionNum", Rcpp::wrap(data.acquisitionNum));
    header.add("msLevel",       Rcpp::wrap(data.msLevel));
    header.add("polarity",       Rcpp::wrap(data.polarity));
    header.add("peaksCount",       Rcpp::wrap(data.peaksCount));
    header.add("totIonCurrent",       Rcpp::wrap(data.totIonCurrent));
    header.add("retentionTime",       Rcpp::wrap(data.retentionTime));
    header.add("basePeakMZ",       Rcpp::wrap(data.basePeakMZ));
    header.add("basePeakIntensity",       Rcpp::wrap(data.basePeakIntensity));
    header.add("collisionEnergy",       Rcpp::wrap(data.collisionEnergy));
    header.add("ionisationEnergy",       Rcpp::wrap(data.ionisationEnergy));
    header.add("lowMZ",       Rcpp::wrap(data.lowMZ));
    header.add("highMZ",       Rcpp::wrap(data.highMZ));
    header.add("precursorScanNum",       Rcpp::wrap(data.precursorScanNum));
    header.add("precursorMZ",       Rcpp::wrap(data.precursorMZ));
    header.add("precursorCharge",       Rcpp::wrap(data.precursorCharge));
    header.add("precursorIntensity",       Rcpp::wrap(data.precursorIntensity));
    //  header.add("scanType",       Rcpp::wrap(data.scanType));
    //  header.add("activationMethod",       Rcpp::wrap(data.activationMethod));
    //  header.add("possibleCharges",       Rcpp::wrap(data.possibleCharges));
    // header.add("numPossibleCharges",       Rcpp::wrap(data.numPossibleCharges));
    //  header.add("possibleChargesArray",       Rcpp::wrap(data.possibleChargesArray));
    header.add("mergedScan",       Rcpp::wrap(data.mergedScan));
    header.add("mergedResultScanNum",       Rcpp::wrap(data.mergedResultScanNum));
    header.add("mergedResultStartScanNum",       Rcpp::wrap(data.mergedResultStartScanNum));
    header.add("mergedResultEndScanNum",       Rcpp::wrap(data.mergedResultEndScanNum));
    //			     header.add("filePosition",       data.filePosition


  return  header;
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return Rcpp::List::create( );
}

Rcpp::DataFrame
RcppRamp::getAllScanHeaderInfo ( ) {
  if (ramp != NULL) {
    if (!isInCacheAllScanHeaderInfo) {
      // Rprintf("Read from disk.\n ");
      int N = ramp->getLastScan();
      rampScanInfo *info = ramp->getScanHeaderInfo( 1 );
      ScanHeaderStruct scanHeader;
      Rcpp::IntegerVector seqNum(N); // number in sequence observed file (1-based)
      Rcpp::IntegerVector acquisitionNum(N); // scan number as declared in File (may be gaps)
      Rcpp::IntegerVector  msLevel(N);
      Rcpp::IntegerVector  polarity(N);
      Rcpp::IntegerVector  peaksCount(N);
      Rcpp::NumericVector totIonCurrent(N);
      Rcpp::NumericVector retentionTime(N);        /* in seconds */
      Rcpp::NumericVector basePeakMZ(N);
      Rcpp::NumericVector basePeakIntensity(N);
      Rcpp::NumericVector collisionEnergy(N);
      Rcpp::NumericVector ionisationEnergy(N);
      Rcpp::NumericVector lowMZ(N);
      Rcpp::NumericVector highMZ(N);
      Rcpp::IntegerVector precursorScanNum(N); /* only if MS level > 1 */
      Rcpp::NumericVector precursorMZ(N);  /* only if MS level > 1 */
      Rcpp::IntegerVector precursorCharge(N);  /* only if MS level > 1 */
      Rcpp::NumericVector precursorIntensity(N);  /* only if MS level > 1 */
      // char scanType[SCANTYPE_LENGTH];
      // char activationMethod[SCANTYPE_LENGTH];
      // char possibleCharges[SCANTYPE_LENGTH];
      // int numPossibleCharges;
      // bool possibleChargesArray[CHARGEARRAY_LENGTH]; /* NOTE: does NOT include "precursorCharge" information; only from "possibleCharges" */
      Rcpp::IntegerVector mergedScan(N);  /* only if MS level > 1 */
      Rcpp::IntegerVector mergedResultScanNum(N); /* scan number of the resultant merged scan */
      Rcpp::IntegerVector mergedResultStartScanNum(N); /* smallest scan number of the scanOrigin for merged scan */
      Rcpp::IntegerVector mergedResultEndScanNum(N); /* largest scan number of the scanOrigin for merged scan */
      
      for (int whichScan=1; whichScan <= N; whichScan++) {
		readHeader(ramp->m_handle, ramp->m_scanOffsets[whichScan], &scanHeader);
		seqNum[whichScan-1] = scanHeader.seqNum;
		acquisitionNum[whichScan-1] = scanHeader.acquisitionNum;
		msLevel[whichScan-1] = scanHeader.msLevel;
		polarity[whichScan-1] = scanHeader.polarity;
		peaksCount[whichScan-1] = scanHeader.peaksCount;
		totIonCurrent[whichScan-1] = scanHeader.totIonCurrent;
		retentionTime[whichScan-1] = scanHeader.retentionTime;
		basePeakMZ[whichScan-1] = scanHeader.basePeakMZ;
		basePeakIntensity[whichScan-1] = scanHeader.basePeakIntensity;
		collisionEnergy[whichScan-1] = scanHeader.collisionEnergy;
		ionisationEnergy[whichScan-1] = scanHeader.ionisationEnergy;
		lowMZ[whichScan-1] = scanHeader.lowMZ;
		highMZ[whichScan-1] = scanHeader.highMZ;
		precursorScanNum[whichScan-1] = scanHeader.precursorScanNum;
		precursorMZ[whichScan-1] = scanHeader.precursorMZ;
		precursorCharge[whichScan-1] = scanHeader.precursorCharge;
		precursorIntensity[whichScan-1] = scanHeader.precursorIntensity;
		mergedScan[whichScan-1] = scanHeader.mergedScan;
		mergedResultScanNum[whichScan-1] = scanHeader.mergedResultScanNum;
		mergedResultStartScanNum[whichScan-1] = scanHeader.mergedResultStartScanNum;
		mergedResultEndScanNum[whichScan-1] = scanHeader.mergedResultEndScanNum;
      }

      ListBuilder header;
      header.add("seqNum", seqNum);
      header.add("acquisitionNum",           acquisitionNum);
      header.add("msLevel",                  msLevel);
      header.add("polarity",               polarity);
      header.add("peaksCount",               peaksCount);
      header.add("totIonCurrent",            totIonCurrent);
      header.add("retentionTime",            retentionTime);
      header.add("basePeakMZ",               basePeakMZ);
      header.add("basePeakIntensity",        basePeakIntensity);
      header.add("collisionEnergy",          collisionEnergy);
      header.add("ionisationEnergy",         ionisationEnergy);
      header.add("lowMZ",                    lowMZ);
      header.add("highMZ",                   highMZ);
      header.add("precursorScanNum",         precursorScanNum);
      header.add("precursorMZ",              precursorMZ);
      header.add("precursorCharge",          precursorCharge);
      header.add("precursorIntensity",       precursorIntensity);
      //  header.add("scanType",                 scanType);
      //  header.add("activationMethod",         activationMethod);
      //  header.add("possibleCharges",          possibleCharges);
      //  header.add("numPossibleCharges",       numPossibleCharges);
						  //  header.add("possibleChargesArray",     possibleChargesArray);
      header.add("mergedScan",               mergedScan);
      header.add("mergedResultScanNum",      mergedResultScanNum);
      header.add("mergedResultStartScanNum", mergedResultStartScanNum);
      header.add("mergedResultEndScanNum",   mergedResultEndScanNum);

      allScanHeaderInfo = header.get();
      isInCacheAllScanHeaderInfo = TRUE;
    } else {
      // Rprintf("Read from cache.\n ");
    }
    return(allScanHeaderInfo);
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return Rcpp::DataFrame::create( );
}

Rcpp::List
RcppRamp::getPeakList ( int whichScan ) {
  if (ramp != NULL) {
    if ((whichScan <= 0) || (whichScan > ramp->getLastScan())) {
      Rprintf("Index whichScan out of bounds [1 ... %d].\n", ramp->getLastScan());
      return Rcpp::List::create( );
    }
    rampPeakList *pl = ramp->getPeakList( whichScan );
    int peaksCount = 0;
    if (pl != NULL) {
      peaksCount = pl->getPeakCount();
    } 
    Rcpp::NumericMatrix peaks(peaksCount, 2);
    if (pl != NULL) {
      rampPeakInfoStruct *peak;
      peak = pl->getPeak(0);
      peaks(0,0) = peak->mz;
      peaks(0,1) = peak->intensity;
      for (int i=1; i < peaksCount; i++) {
	peak++;
	peaks(i,0) = peak->mz;
	peaks(i,1) = peak->intensity;
      }
      delete pl;
    }
    return Rcpp::List::create(
			      Rcpp::_["peaksCount"]  = peaksCount,
			      Rcpp::_["peaks"]  = peaks
			     ) ;
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return Rcpp::List::create( );
}

Rcpp::NumericMatrix
RcppRamp::get3DMap ( std::vector<int> scanNumbers, double whichMzLow, double whichMzHigh, double resMz ) {
  if (ramp != NULL) {
    double f = 1 / resMz;
    int low = round(whichMzLow * f);
    int high = round(whichMzHigh * f);
    int dmz = high - low + 1;
    int drt = scanNumbers.size();
    Rcpp::NumericMatrix map3d(drt, dmz);
    for (int i = 0; i < drt; i++) {
      for (int j = 0; j < dmz; j++) {
	map3d(i,j) = 0.0;
      }
    }
    // map3d = 0.0;
    int j=0;
    Rprintf("%d\n",1);
    for (int i = 0; i < scanNumbers.size(); i++) {
      rampPeakList *pl = ramp->getPeakList( scanNumbers[i] );
      int peaksCount = pl->getPeakCount();
      rampPeakInfoStruct *peak;
      peak = pl->getPeak(0);
      j = round(peak->mz * f) - low;
      if ((j >= 0) & (j < dmz)) {
	if (peak->intensity > map3d(i,j)) {
	  map3d(i,j) = peak->intensity;
	}
      }
      for (int k=1; k < peaksCount; k++) {
	peak++;
	j = round(peak->mz * f) - low;
	if ((j >= 0) & (j < dmz)) {
	  if (peak->intensity > map3d(i,j)) {
	    map3d(i,j) = peak->intensity;
	  }
	}
      }
      delete pl;
    }
    return(map3d);
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return Rcpp::NumericMatrix(0,0);
}

int
RcppRamp::getLastScan() const {
  if (ramp != NULL) {
    return ramp->getLastScan();
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return -1;
}

bool
RcppRamp::OK (  ) {
  if (ramp != NULL) {
    return ramp->OK();
  }
  // Rprintf("Warning: Ramp not yet initialized.\n ");
  return false;
}

