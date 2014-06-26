#include "RcppPwiz.h"

RcppPwiz::RcppPwiz() {
  msd = NULL;
  runInfo = Rcpp::List::create();
  isInCacheRunInfo = FALSE;
  instrumentInfo = Rcpp::List::create();
  chromatogramsInfo = Rcpp::List::create();
  isInCacheInstrumentInfo = FALSE;
  allScanHeaderInfo = Rcpp::List::create();
  isInCacheAllScanHeaderInfo = FALSE;
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
  Rprintf("Warning: pwiz not yet initialized.\n ");
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
		if ((whichScan <= 0) || (whichScan > slp->size())) {
			Rprintf("Index out of bounds [1 ... %d].\n", slp->size());
			return Rcpp::List::create( );
		}
		SpectrumInfo info;
		info.update(*msd->run.spectrumListPtr->spectrum(whichScan - 1));
		SpectrumPtr s = slp->spectrum(whichScan - 1, true);
		vector<MZIntensityPair> pairs;
		s->getMZIntensityPairs(pairs);
		
		if(info.msLevel == 1){
			return Rcpp::List::create(
				Rcpp::_["id"]				= info.id,
				Rcpp::_["index"]			= info.index,
			    Rcpp::_["scanNumber"]		= info.scanNumber,
			    Rcpp::_["msLevel"]			= info.msLevel,
			    Rcpp::_["retentionTime"]    = info.retentionTime,
			    Rcpp::_["isZoomScan"]      	= info.isZoomScan,
			    Rcpp::_["filterString"]     = info.filterString,
			    Rcpp::_["peaksCount"]     	= pairs.size(),
			    Rcpp::_["lowMZ"]   			= info.mzLow,
			    Rcpp::_["highMZ"]      		= info.mzHigh,
			    Rcpp::_["massAnalyzerType"]	= info.massAnalyzerTypeAbbreviation(),
			    Rcpp::_["basePeakMZ"]		= info.basePeakMZ,
			    Rcpp::_["basePeakIntensity"]= info.basePeakIntensity,
			    Rcpp::_["totalIonCurrent"]	= info.totalIonCurrent,
			    Rcpp::_["thermoMonoisotopicMZ"]	= info.thermoMonoisotopicMZ,
			    Rcpp::_["ionInjectionTime"]	= info.ionInjectionTime);
		}else{
			return Rcpp::List::create(
				Rcpp::_["id"]				= info.id,
				Rcpp::_["index"]			= info.index,
			    Rcpp::_["scanNumber"]		= info.scanNumber,
			    Rcpp::_["msLevel"]			= info.msLevel,
			    Rcpp::_["retentionTime"]    = info.retentionTime,
			    Rcpp::_["isZoomScan"]      	= info.isZoomScan,
			    Rcpp::_["filterString"]     = info.filterString,
			    Rcpp::_["peaksCount"]     	= pairs.size(),
			    Rcpp::_["lowMZ"]   			= info.mzLow,
			    Rcpp::_["highMZ"]      		= info.mzHigh,
			    Rcpp::_["massAnalyzerType"]	= info.massAnalyzerTypeAbbreviation(),
			    Rcpp::_["basePeakMZ"]		= info.basePeakMZ,
			    Rcpp::_["basePeakIntensity"]= info.basePeakIntensity,
			    Rcpp::_["totalIonCurrent"]	= info.totalIonCurrent,
			    Rcpp::_["thermoMonoisotopicMZ"]	= info.thermoMonoisotopicMZ,
			    Rcpp::_["ionInjectionTime"]	= info.ionInjectionTime,
			    Rcpp::_["precursorIndex"]	= info.precursors[0].index,
			    Rcpp::_["precursorMZ"]		= info.precursors[0].mz,
			    Rcpp::_["precursorCharge"]	= info.precursors[0].charge,
			    Rcpp::_["precursorIntensity"]=info.precursors[0].intensity);
		}	
	}else{
		Rprintf("Warning: pwiz not yet initialized.\n ");
		return Rcpp::List::create( );
	}
}


Rcpp::DataFrame RcppPwiz::getAllScanHeaderInfo ( ) {
	if (msd != NULL) {
		if (!isInCacheAllScanHeaderInfo) {
			SpectrumListPtr slp = msd->run.spectrumListPtr;
			int N = slp->size();
			
			Rcpp::StringVector id(N);
			Rcpp::IntegerVector index(N);
			Rcpp::IntegerVector scanNumber(N);
			Rcpp::IntegerVector msLevel(N);
			Rcpp::NumericVector retentionTime(N);
			Rcpp::LogicalVector isZoomScan(N);
			Rcpp::NumericVector lowMZ(N);
			Rcpp::NumericVector highMZ(N);
			Rcpp::NumericVector peaksCount(N);
			Rcpp::StringVector filterString(N);
			Rcpp::StringVector massAnalyzerType(N);
			Rcpp::NumericVector basePeakMZ(N);
			Rcpp::NumericVector basePeakIntensity(N);
			Rcpp::NumericVector totalIonCurrent(N);
			Rcpp::NumericVector thermoMonoisotopicMZ(N);
			Rcpp::NumericVector ionInjectionTime(N);
			Rcpp::NumericVector precursorIndex(N);
			Rcpp::NumericVector precursorMZ(N);
			Rcpp::NumericVector precursorCharge(N);
			Rcpp::NumericVector precursorIntensity(N);
			
			SpectrumInfo info;
			SpectrumPtr s;
			vector<MZIntensityPair> pairs;
			
			for (int whichScan=1; whichScan <= N; whichScan++) {
				info.update(*msd->run.spectrumListPtr->spectrum(whichScan - 1));
				s = slp->spectrum(whichScan - 1, true);
				s->getMZIntensityPairs(pairs);
				id[whichScan-1] = info.id;
				index[whichScan-1] = info.index;
				scanNumber[whichScan-1] = info.scanNumber;
				msLevel[whichScan-1] = info.msLevel;
				retentionTime[whichScan-1] = info.retentionTime;
				isZoomScan[whichScan-1] = info.isZoomScan;
				lowMZ[whichScan-1] = info.mzLow;
				highMZ[whichScan-1] = info.mzHigh;
				peaksCount[whichScan-1] = pairs.size();
				filterString[whichScan-1] = info.filterString;
				massAnalyzerType[whichScan-1] = info.massAnalyzerTypeAbbreviation();
				basePeakMZ[whichScan-1] = info.basePeakMZ;
				basePeakIntensity[whichScan-1] = info.basePeakIntensity;
				totalIonCurrent[whichScan-1] = info.totalIonCurrent;
				thermoMonoisotopicMZ[whichScan-1] = info.thermoMonoisotopicMZ;
				ionInjectionTime[whichScan-1] = info.ionInjectionTime;

				if(info.msLevel > 1){
					precursorIndex[whichScan-1] = info.precursors[0].index;
					precursorMZ[whichScan-1] = info.precursors[0].mz;
					precursorCharge[whichScan-1] = info.precursors[0].charge;
					precursorIntensity[whichScan-1] = info.precursors[0].intensity;
				}else {
					precursorIndex[whichScan-1] = 0;
					precursorMZ[whichScan-1] = 0;
					precursorCharge[whichScan-1] = 0;
					precursorIntensity[whichScan-1] = 0;
				}

			}
			
			allScanHeaderInfo = Rcpp::DataFrame::create( 
								Rcpp::_["id"] = id,
								Rcpp::_["index"] = index,
								Rcpp::_["scanNumber"] = scanNumber,
								Rcpp::_["msLevel"] = msLevel,
								Rcpp::_["retentionTime"] = retentionTime,
								Rcpp::_["isZoomScan"] = isZoomScan,
								Rcpp::_["lowMZ"] = lowMZ,
								Rcpp::_["highMZ"] = highMZ,
								Rcpp::_["peaksCount"] = peaksCount,
								Rcpp::_["filterString"] = filterString,
								Rcpp::_["massAnalyzerType"] = massAnalyzerType,
								Rcpp::_["basePeakMZ"] = basePeakMZ,
								Rcpp::_["basePeakIntensity"] = basePeakIntensity,
								Rcpp::_["totalIonCurrent"] = totalIonCurrent,
								Rcpp::_["thermoMonoisotopicMZ"] = thermoMonoisotopicMZ,
								Rcpp::_["ionInjectionTime"] = ionInjectionTime,
								Rcpp::_["precursorIndex"] = precursorIndex,
								Rcpp::_["precursorMZ"] = precursorMZ,
								Rcpp::_["precursorCharge"] = precursorCharge,
								Rcpp::_["precursorIntensity"] = precursorIntensity);
			  isInCacheAllScanHeaderInfo = TRUE;
		}
		return(allScanHeaderInfo);
	}
	Rprintf("Warning: pwiz not yet initialized.\n ");
	return Rcpp::DataFrame::create( );
}

Rcpp::List RcppPwiz::getPeakList ( int whichScan ) {
	if (msd != NULL) {
		SpectrumListPtr slp = msd->run.spectrumListPtr;

		if ((whichScan <= 0) || (whichScan > slp->size())) {
			Rprintf("Index whichScan out of bounds [1 ... %d].\n", slp->size());
			return Rcpp::List::create( );
		}
		
		SpectrumPtr s = slp->spectrum(whichScan - 1, true);
		vector<MZIntensityPair> pairs;
		s->getMZIntensityPairs(pairs);

		Rcpp::NumericMatrix peaks(pairs.size(), 2);
		
		if(pairs.size()!=0){
			for (int i = 0; i < pairs.size(); i++) {
				MZIntensityPair p = pairs.at(i);
				peaks(i,0) = p.mz;
				peaks(i,1) = p.intensity;
			}
			
		}

    return Rcpp::List::create(
			      Rcpp::_["peaksCount"]  = pairs.size(),
			      Rcpp::_["peaks"]  = peaks
			     ) ;
  }
  Rprintf("Warning: pwiz not yet initialized.\n ");
  return Rcpp::List::create( );
}

Rcpp::List RcppPwiz::getChromatogramsInfo() {
  if (msd != NULL) {
    ChromatogramListPtr clp = msd->run.chromatogramListPtr;
	Rcpp::Rcout << clp->size() << endl;
    ChromatogramPtr c = clp->chromatogram(0, true);  
    vector<TimeIntensityPair> pairs;
    c->getTimeIntensityPairs(pairs);
      
    int N = pairs.size();
    Rcpp::NumericVector time(N);
    Rcpp::NumericVector intensity(N); 

    for(int i = 0; i < pairs.size(); i++)
    {
        TimeIntensityPair p = pairs.at(i);
        time[i] = p.time;
        intensity[i] = p.intensity;     
          
      }

    return Rcpp::List::create(
			    Rcpp::_["time"]	  = time,
			    Rcpp::_["intensity"]  = intensity);
    
  }
  Rprintf("Warning: pwiz not yet initialized.\n ");
  return Rcpp::List::create( );
}

Rcpp::NumericMatrix RcppPwiz::get3DMap ( std::vector<int> scanNumbers, double whichMzLow, double whichMzHigh, double resMz ) {
	if (msd != NULL) {
		
		SpectrumListPtr slp = msd->run.spectrumListPtr;
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

		int j=0;
		Rprintf("%d\n",1);
		for (int i = 0; i < scanNumbers.size(); i++) {
			SpectrumPtr s = slp->spectrum(scanNumbers[i] - 1, true);
			vector<MZIntensityPair> pairs;
			s->getMZIntensityPairs(pairs);
			

			for (int k=0; k < pairs.size(); k++) {
				MZIntensityPair p = pairs.at(k);
				j = round(p.mz * f) - low;
				if ((j >= 0) & (j < dmz)) {
					if (p.intensity > map3d(i,j)) {
						map3d(i,j) = p.intensity;
					}
				}
			}

		}
		return(map3d);
	}
	
	Rprintf("Warning: pwiz not yet initialized.\n ");
	return Rcpp::NumericMatrix(0,0);
}
