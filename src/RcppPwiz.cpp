#include "RcppPwiz.h"

RcppPwiz::RcppPwiz()
{
    msd = NULL;
    instrumentInfo = Rcpp::List::create();
    chromatogramsInfo = Rcpp::DataFrame::create();
    isInCacheInstrumentInfo = FALSE;
    allScanHeaderInfo = Rcpp::List::create();
    isInCacheAllScanHeaderInfo = FALSE;
}
/* Destructor*/
RcppPwiz::~RcppPwiz()
{
  RcppPwiz::close();
}

void RcppPwiz::open(const string& fileName)
{

    filename = fileName;
    msd = new MSDataFile(fileName);

}

/* Release all memory on close. */
void RcppPwiz::close()
{
  if (msd != NULL)
    {
      delete msd;
      msd = NULL;
      instrumentInfo = Rcpp::List::create();
      chromatogramsInfo = Rcpp::DataFrame::create();
      isInCacheInstrumentInfo = FALSE;
      allScanHeaderInfo = Rcpp::List::create();
      isInCacheAllScanHeaderInfo = FALSE;
    }
}

/*
void RcppPwiz::writeMSfile(const string& file, const string& format)
{
    if (msd != NULL)
    {
        if(format == "mgf")
        {
            std::ofstream* mgfOutFileP = new std::ofstream(file.c_str());
            Serializer_MGF serializerMGF;
            serializerMGF.write(*mgfOutFileP, *msd);
            mgfOutFileP->flush();
            mgfOutFileP->close();
        }
        else if(format == "mzxml")
        {
            std::ofstream mzXMLOutFileP(file.c_str());
            Serializer_mzXML::Config config;
            config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
            Serializer_mzXML serializerMzXML(config);
            serializerMzXML.write(mzXMLOutFileP, *msd);
        }
        else if(format == "mzml")
        {
            std::ofstream mzXMLOutFileP(file.c_str());
            Serializer_mzML::Config config;
            config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
            Serializer_mzML mzmlSerializer(config);
            mzmlSerializer.write(mzXMLOutFileP, *msd);
        }
        else
            Rcpp::Rcerr << format << " format not supported! Please try mgf, mzML, mzXML or mz5." << std::endl;
    }
    else
        Rcpp::Rcerr << "No pwiz object available! Please open a file first!" << std::endl;
}
*/

string RcppPwiz::getFilename() {
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

int RcppPwiz::getLastChrom() const {
    if (msd != NULL) {
      ChromatogramListPtr clp = msd->run.chromatogramListPtr;
      return clp->size();
    }
    Rprintf("Warning: pwiz not yet initialized.\n ");
    return -1;
}

Rcpp::List RcppPwiz::getInstrumentInfo ( )
{
    if (msd != NULL)
    {
        if (!isInCacheInstrumentInfo)
        {

            vector<InstrumentConfigurationPtr> icp = msd->instrumentConfigurationPtrs; // NULL for mzData

            if (icp.size() != 0)
            {

                CVTranslator cvTranslator;
                LegacyAdapter_Instrument adapter(*icp[0], cvTranslator);
                vector<SoftwarePtr> sp = msd->softwarePtrs;
                std::vector<SamplePtr> sample = msd->samplePtrs;
                std::vector<ScanSettingsPtr> scansetting = msd->scanSettingsPtrs;
                instrumentInfo = Rcpp::List::create(
                                     Rcpp::_["manufacturer"]  = std::string(adapter.manufacturer()),
                                     Rcpp::_["model"]         = std::string(adapter.model()),
                                     Rcpp::_["ionisation"]    = std::string(adapter.ionisation()),
                                     Rcpp::_["analyzer"]      = std::string(adapter.analyzer()),
                                     Rcpp::_["detector"]      = std::string(adapter.detector()),
                                     Rcpp::_["software"]      = sp[0]->id + " " + sp[0]->version,
                                     Rcpp::_["sample"]		  = (sample.size()>0?sample[0]->name+sample[0]->id:""),
                                     Rcpp::_["source"]        = (scansetting.size()>0?scansetting[0]->sourceFilePtrs[0]->location:"")
                                 ) ;

            }
            else
            {
                instrumentInfo = Rcpp::List::create(
                                     Rcpp::_["manufacturer"]  = "",
                                     Rcpp::_["model"]         = "",
                                     Rcpp::_["ionisation"]    = "",
                                     Rcpp::_["analyzer"]      = "",
                                     Rcpp::_["detector"]      = "",
                                     Rcpp::_["software"]      = "",
                                     Rcpp::_["sample"]		  = "",
                                     Rcpp::_["source"]		  = ""
                                 ) ;
            }

            isInCacheInstrumentInfo = TRUE;
        }
        return(instrumentInfo);
    }
    Rprintf("Warning: pwiz not yet initialized.\n ");
    return instrumentInfo;
}


Rcpp::DataFrame RcppPwiz::getScanHeaderInfo (Rcpp::IntegerVector whichScan)
{
    if (msd != NULL)
    {
      SpectrumListPtr slp = msd->run.spectrumListPtr;
      int N = slp->size();
      
      int N_scans = whichScan.size();
      
      ScanHeaderStruct scanHeader;
      RAMPAdapter * adapter = new  RAMPAdapter(filename);
      Rcpp::IntegerVector seqNum(N_scans); // number in sequence observed file (1-based)
      Rcpp::IntegerVector acquisitionNum(N_scans); // scan number as declared in File (may be gaps)
      Rcpp::IntegerVector msLevel(N_scans);
      Rcpp::IntegerVector polarity(N_scans);
      Rcpp::IntegerVector peaksCount(N_scans);
      Rcpp::NumericVector totIonCurrent(N_scans);
      Rcpp::NumericVector retentionTime(N_scans);        /* in seconds */
      Rcpp::NumericVector basePeakMZ(N_scans);
      Rcpp::NumericVector basePeakIntensity(N_scans);
      Rcpp::NumericVector collisionEnergy(N_scans);
      Rcpp::NumericVector ionisationEnergy(N_scans);
      Rcpp::NumericVector lowMZ(N_scans);
      Rcpp::NumericVector highMZ(N_scans);
      Rcpp::IntegerVector precursorScanNum(N_scans); /* only if MS level > 1 */
      Rcpp::NumericVector precursorMZ(N_scans);  /* only if MS level > 1 */
      Rcpp::IntegerVector precursorCharge(N_scans);  /* only if MS level > 1 */
      Rcpp::NumericVector precursorIntensity(N_scans);  /* only if MS level > 1 */
      //char scanType[SCANTYPE_LENGTH];
      //char activationMethod[SCANTYPE_LENGTH];
      //char possibleCharges[SCANTYPE_LENGTH];
      //int numPossibleCharges;
      //bool possibleChargesArray[CHARGEARRAY_LENGTH]; /* NOTE: does NOT include "precursorCharge" information; only from "possibleCharges" */
      Rcpp::IntegerVector mergedScan(N_scans);  /* only if MS level > 1 */
      Rcpp::IntegerVector mergedResultScanNum(N_scans); /* scan number of the resultant merged scan */
      Rcpp::IntegerVector mergedResultStartScanNum(N_scans); /* smallest scan number of the scanOrigin for merged scan */
      Rcpp::IntegerVector mergedResultEndScanNum(N_scans); /* largest scan number of the scanOrigin for merged scan */
      
      for (int i = 0; i < N_scans; i++)
	{
	  int current_scan = whichScan[i];
	  adapter->getScanHeader(current_scan - 1, scanHeader, false);
	  seqNum[i] = scanHeader.seqNum;
	  acquisitionNum[i] = scanHeader.acquisitionNum;
	  msLevel[i] = scanHeader.msLevel;
	  
	  SpectrumPtr sp = slp->spectrum(current_scan-1, false); // Is TRUE neccessary here ? 
	  CVParam param = sp->cvParamChild(MS_scan_polarity);
	  polarity[i] = (param.cvid==MS_negative_scan ? 0 : (param.cvid==MS_positive_scan ? +1 : -1 ) );
	  
	  peaksCount[i] = scanHeader.peaksCount;
	  totIonCurrent[i] = scanHeader.totIonCurrent;
	  retentionTime[i] = scanHeader.retentionTime;
	  basePeakMZ[i] = scanHeader.basePeakMZ;
	  basePeakIntensity[i] = scanHeader.basePeakIntensity;
	  collisionEnergy[i] = scanHeader.collisionEnergy;
	  ionisationEnergy[i] = scanHeader.ionisationEnergy;
	  lowMZ[i] = scanHeader.lowMZ;
	  highMZ[i] = scanHeader.highMZ;
	  precursorScanNum[i] = scanHeader.precursorScanNum;
	  precursorMZ[i] = scanHeader.precursorMZ;
	  precursorCharge[i] = scanHeader.precursorCharge;
	  precursorIntensity[i] = scanHeader.precursorIntensity;
	  mergedScan[i] = scanHeader.mergedScan;
	  mergedResultScanNum[i] = scanHeader.mergedResultScanNum;
	  mergedResultStartScanNum[i] = scanHeader.mergedResultStartScanNum;
	  mergedResultEndScanNum[i] = scanHeader.mergedResultEndScanNum;
	}
      // delete adapter issue #64
      delete adapter;
      adapter = NULL;
      
      Rcpp::List header(21);
      std::vector<std::string> names;
      int i = 0;
      names.push_back("seqNum");
      header[i++] = Rcpp::wrap(seqNum);
      names.push_back("acquisitionNum");
      header[i++] = Rcpp::wrap(acquisitionNum);
      names.push_back("msLevel");
      header[i++] = Rcpp::wrap(msLevel);
      names.push_back("polarity");
      header[i++] = Rcpp::wrap(polarity);
      names.push_back("peaksCount");
      header[i++] = Rcpp::wrap(peaksCount);
      names.push_back("totIonCurrent");
      header[i++] = Rcpp::wrap(totIonCurrent);
      names.push_back("retentionTime");
      header[i++] = Rcpp::wrap(retentionTime);
      names.push_back("basePeakMZ");
      header[i++] = Rcpp::wrap(basePeakMZ);
      names.push_back("basePeakIntensity");
      header[i++] = Rcpp::wrap(basePeakIntensity);
      names.push_back("collisionEnergy");
      header[i++] = Rcpp::wrap(collisionEnergy);
      names.push_back("ionisationEnergy");
      header[i++] = Rcpp::wrap(ionisationEnergy);
      names.push_back("lowMZ");
      header[i++] = Rcpp::wrap(lowMZ);
      names.push_back("highMZ");
      header[i++] = Rcpp::wrap(highMZ);
      names.push_back("precursorScanNum");
      header[i++] = Rcpp::wrap(precursorScanNum);
      names.push_back("precursorMZ");
      header[i++] = Rcpp::wrap(precursorMZ);
      names.push_back("precursorCharge");
      header[i++] = Rcpp::wrap(precursorCharge);
      names.push_back("precursorIntensity");
      header[i++] = Rcpp::wrap(precursorIntensity);
      names.push_back("mergedScan");
      header[i++] = Rcpp::wrap(mergedScan);
      names.push_back("mergedResultScanNum");
      header[i++] = Rcpp::wrap(mergedResultScanNum);
      names.push_back("mergedResultStartScanNum");
      header[i++] = Rcpp::wrap(mergedResultStartScanNum);
      names.push_back("mergedResultEndScanNum");
      header[i++] = Rcpp::wrap(mergedResultEndScanNum);
      
      header.attr("names") = names;
      
      return(header);
    }
    Rprintf("Warning: pwiz not yet initialized.\n ");
    return Rcpp::DataFrame::create( );
}

Rcpp::DataFrame RcppPwiz::getAllScanHeaderInfo ( )
{
    if (msd != NULL)
    {
        if (!isInCacheAllScanHeaderInfo)
        {
            SpectrumListPtr slp = msd->run.spectrumListPtr;
            int N = slp->size();

            allScanHeaderInfo = getScanHeaderInfo(Rcpp::seq(1, N));
            isInCacheAllScanHeaderInfo = TRUE;	    
        }
        return(allScanHeaderInfo);
    }
    Rprintf("Warning: pwiz not yet initialized.\n ");
    return Rcpp::DataFrame::create( );
}

Rcpp::List RcppPwiz::getPeakList ( int whichScan )
{
    if (msd != NULL)
    {
        SpectrumListPtr slp = msd->run.spectrumListPtr;

        if ((whichScan <= 0) || (whichScan > slp->size()))
        {
            Rprintf("Index whichScan out of bounds [1 ... %d].\n", slp->size());
            return Rcpp::List::create( );
        }

        SpectrumPtr s = slp->spectrum(whichScan - 1, true);
        vector<MZIntensityPair> pairs;
        s->getMZIntensityPairs(pairs);

        Rcpp::NumericMatrix peaks(pairs.size(), 2);

        if(pairs.size()!=0)
        {
            for (int i = 0; i < pairs.size(); i++)
            {
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

Rcpp::DataFrame RcppPwiz::getChromatogramsInfo( int whichChrom )
{
    if (msd != NULL) {
      ChromatogramListPtr clp = msd->run.chromatogramListPtr;
      if (clp.get() == 0) {
	Rcpp::Rcerr << "The direct support for chromatogram info is only available in mzML format." << std::endl;
	return Rcpp::DataFrame::create();
      } else if (clp->size() == 0) {
	Rcpp::Rcerr << "No available chromatogram info." << std::endl;
	return Rcpp::DataFrame::create();
      } else if ( (whichChrom < 0) || (whichChrom > clp->size()) ) {
	Rprintf("Index whichChrom out of bounds [0 ... %d].\n", (clp->size())-1);
	return Rcpp::DataFrame::create( );
      } else {
	std::vector<double> time;
	std::vector<double> intensity;
	ChromatogramPtr c = clp->chromatogram(whichChrom, true);
	vector<TimeIntensityPair> pairs;
	c->getTimeIntensityPairs (pairs);

	for (int i = 0; i < pairs.size(); i++) {
	  TimeIntensityPair p = pairs.at(i);
	  time.push_back(p.time);
	  intensity.push_back(p.intensity);
	}

	chromatogramsInfo = Rcpp::DataFrame::create(Rcpp::_["time"] = time,
						    Rcpp::_[c->id]  = intensity);

      }
      return(chromatogramsInfo);
    }
    Rprintf("Warning: pwiz not yet initialized.\n ");
    return Rcpp::DataFrame::create( );
}

Rcpp::NumericMatrix RcppPwiz::get3DMap ( std::vector<int> scanNumbers, double whichMzLow, double whichMzHigh, double resMz )
{
    if (msd != NULL)
    {

        SpectrumListPtr slp = msd->run.spectrumListPtr;
        double f = 1 / resMz;
        int low = round(whichMzLow * f);
        int high = round(whichMzHigh * f);
        int dmz = high - low + 1;
        int drt = scanNumbers.size();

        Rcpp::NumericMatrix map3d(drt, dmz);

        for (int i = 0; i < drt; i++)
        {
            for (int j = 0; j < dmz; j++)
            {
                map3d(i,j) = 0.0;
            }
        }

        int j=0;
        Rprintf("%d\n",1);
        for (int i = 0; i < scanNumbers.size(); i++)
        {
            SpectrumPtr s = slp->spectrum(scanNumbers[i] - 1, true);
            vector<MZIntensityPair> pairs;
            s->getMZIntensityPairs(pairs);

            for (int k=0; k < pairs.size(); k++)
            {
                MZIntensityPair p = pairs.at(k);
                j = round(p.mz * f) - low;
                if ((j >= 0) & (j < dmz))
                {
                    if (p.intensity > map3d(i,j))
                    {
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
