#include "RcppRamp.h"



RcppRamp::RcppRamp() {
  ramp = NULL;
  runInfo = Rcpp::List::create( );
  isInCacheRunInfo = FALSE;
  instrumentInfo = Rcpp::List::create( );
  isInCacheInstrumentInfo = FALSE;
  allScanHeaderInfo = Rcpp::List::create( );
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
    printf("Failed to open file.\n ");
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
  printf("Warning: Ramp not yet initialized.\n ");
  return filename;
}

Rcpp::List
RcppRamp::getRunInfo (  ) {
  if (ramp != NULL) {
    if (!isInCacheRunInfo) {
      // printf("Read from disk.\n ");
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
      // printf("Read from cache.\n ");
    }
    return runInfo;
  }
  printf("Warning: Ramp not yet initialized.\n");
  return runInfo;
}

Rcpp::List
RcppRamp::getInstrumentInfo ( ) {
  if (ramp != NULL) {
    if (!isInCacheInstrumentInfo) {
      // printf("Read from disk.\n ");
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
      // printf("Read from cache.\n ");
    }
    return(instrumentInfo);
  }
  printf("Warning: Ramp not yet initialized.\n ");
  return instrumentInfo;
}

Rcpp::List
RcppRamp::getScanHeaderInfo ( int whichScan  ) {
  if (ramp != NULL) {
    if ((whichScan <= 0) || (whichScan > ramp->getLastScan())) {
      printf("Index whichScan out of bounds [1 ... %d].\n", ramp->getLastScan());
      return Rcpp::List::create( );
    }
    rampScanInfo *info = ramp->getScanHeaderInfo( whichScan );
    ScanHeaderStruct data = info->m_data;
    delete info;
    return Rcpp::List::create(
			     Rcpp::_["seqNum"]              = data.seqNum,
			     Rcpp::_["acquisitionNum"]      = data.acquisitionNum,
			     Rcpp::_["msLevel"]      = data.msLevel,
			     Rcpp::_["peaksCount"]      = data.peaksCount,
			     Rcpp::_["totIonCurrent"]      = data.totIonCurrent,
			     Rcpp::_["retentionTime"]      = data.retentionTime,
			     Rcpp::_["basePeakMZ"]      = data.basePeakMZ,
			     Rcpp::_["basePeakIntensity"]      = data.basePeakIntensity,
			     Rcpp::_["collisionEnergy"]      = data.collisionEnergy,
			     Rcpp::_["ionisationEnergy"]      = data.ionisationEnergy,
			     Rcpp::_["lowMZ"]      = data.lowMZ,
			     Rcpp::_["highMZ"]      = data.highMZ,
			     Rcpp::_["precursorScanNum"]      = data.precursorScanNum,
			     Rcpp::_["precursorMZ"]      = data.precursorMZ,
			     Rcpp::_["precursorCharge"]      = data.precursorCharge,
			     Rcpp::_["precursorIntensity"]      = data.precursorIntensity,
			     //  Rcpp::_["scanType"]      = data.scanType,
			     //  Rcpp::_["activationMethod"]      = data.activationMethod,
			     //  Rcpp::_["possibleCharges"]      = data.possibleCharges,
			     // Rcpp::_["numPossibleCharges"]      = data.numPossibleCharges,
			     //  Rcpp::_["possibleChargesArray"]      = data.possibleChargesArray,
			     Rcpp::_["mergedScan"]      = data.mergedScan,
			     Rcpp::_["mergedResultScanNum"]      = data.mergedResultScanNum,
			     Rcpp::_["mergedResultStartScanNum"]      = data.mergedResultStartScanNum,
			     Rcpp::_["mergedResultEndScanNum"]      = data.mergedResultEndScanNum
			     //			     Rcpp::_["filePosition"]      = data.filePosition
			     ) ;
  }
  printf("Warning: Ramp not yet initialized.\n ");
  return Rcpp::List::create( );
}

Rcpp::DataFrame
RcppRamp::getAllScanHeaderInfo ( ) {
  if (ramp != NULL) {
    if (!isInCacheAllScanHeaderInfo) {
      // printf("Read from disk.\n ");
      int N = ramp->getLastScan();
      rampScanInfo *info = ramp->getScanHeaderInfo( 1 );
      ScanHeaderStruct scanHeader;
      Rcpp::IntegerVector seqNum(N); // number in sequence observed file (1-based)
      Rcpp::IntegerVector acquisitionNum(N); // scan number as declared in File (may be gaps)
      Rcpp::IntegerVector  msLevel(N);
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
      allScanHeaderInfo = Rcpp::DataFrame::create( 
						  Rcpp::_["seqNum"]                   = seqNum,
						  Rcpp::_["acquisitionNum"]           = acquisitionNum,
						  Rcpp::_["msLevel"]                  = msLevel,
						  Rcpp::_["peaksCount"]               = peaksCount,
						  Rcpp::_["totIonCurrent"]            = totIonCurrent,
						  Rcpp::_["retentionTime"]            = retentionTime,
						  Rcpp::_["basePeakMZ"]               = basePeakMZ,
						  Rcpp::_["basePeakIntensity"]        = basePeakIntensity,
						  Rcpp::_["collisionEnergy"]          = collisionEnergy,
						  Rcpp::_["ionisationEnergy"]         = ionisationEnergy,
						  Rcpp::_["lowMZ"]                    = lowMZ,
						  Rcpp::_["highMZ"]                   = highMZ,
						  Rcpp::_["precursorScanNum"]         = precursorScanNum,
						  Rcpp::_["precursorMZ"]              = precursorMZ,
						  Rcpp::_["precursorCharge"]          = precursorCharge,
						  Rcpp::_["precursorIntensity"]       = precursorIntensity,
						  //  Rcpp::_["scanType"]                 = scanType,
						  //  Rcpp::_["activationMethod"]         = activationMethod,
						  //  Rcpp::_["possibleCharges"]          = possibleCharges,
						  //  Rcpp::_["numPossibleCharges"]       = numPossibleCharges,
						  //  Rcpp::_["possibleChargesArray"]     = possibleChargesArray,
						  Rcpp::_["mergedScan"]               = mergedScan,
						  Rcpp::_["mergedResultScanNum"]      = mergedResultScanNum,
						  Rcpp::_["mergedResultStartScanNum"] = mergedResultStartScanNum,
						  Rcpp::_["mergedResultEndScanNum"]   = mergedResultEndScanNum
						   );
      isInCacheAllScanHeaderInfo = TRUE;
    } else {
      // printf("Read from cache.\n ");
    }
    return(allScanHeaderInfo);
  }
  printf("Warning: Ramp not yet initialized.\n ");
  return Rcpp::DataFrame::create( );
}

Rcpp::List
RcppRamp::getPeakList ( int whichScan ) {
  if (ramp != NULL) {
    if ((whichScan <= 0) || (whichScan > ramp->getLastScan())) {
      printf("Index whichScan out of bounds [1 ... %d].\n", ramp->getLastScan());
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
  printf("Warning: Ramp not yet initialized.\n ");
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
    printf("%d\n",1);
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
  printf("Warning: Ramp not yet initialized.\n ");
  return Rcpp::NumericMatrix(0,0);
}

int
RcppRamp::getLastScan() const {
  if (ramp != NULL) {
    return ramp->getLastScan();
  }
  printf("Warning: Ramp not yet initialized.\n ");
  return -1;
}

bool
RcppRamp::OK (  ) {
  if (ramp != NULL) {
    return ramp->OK();
  }
  // printf("Warning: Ramp not yet initialized.\n ");
  return false;
}

