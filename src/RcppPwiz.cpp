#include "RcppPwiz.h"

RcppPwiz::RcppPwiz()
{
    msd = NULL;
    runInfo = Rcpp::List::create();
    isInCacheRunInfo = FALSE;
    instrumentInfo = Rcpp::List::create();
    chromatogramsInfo = Rcpp::List::create();
    isInCacheInstrumentInfo = FALSE;
    allScanHeaderInfo = Rcpp::List::create();
    isInCacheAllScanHeaderInfo = FALSE;
}

void RcppPwiz::open(const string& fileName)
{

    filename = fileName;
    msd = new MSDataFile(fileName);

}

void RcppPwiz::writeMSfile(const string& file, const string& format)
{
	if (msd != NULL){
		if(format == "mgf"){
			MSDataFile::WriteConfig writeConfig;
			Serializer_mz5 serializer(writeConfig);
			IterationListenerRegistry ilr;
			serializer.write(file, *msd, &ilr);
		}else if(format == "mz5"){

			Serializer_mz5 serializerMz5;
			serializerMz5.write(file, *msd);
		}else if(format == "mzxml"){
			std::ofstream mzXMLOutFileP(file.c_str()); 
			Serializer_mzXML::Config config;
			config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
			Serializer_mzXML serializerMzXML(config);
			serializerMzXML.write(mzXMLOutFileP, *msd);
		}else if(format == "mzml"){
			std::ofstream mzXMLOutFileP(file.c_str()); 
			Serializer_mzML::Config config;
			config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
			Serializer_mzML mzmlSerializer(config);
			mzmlSerializer.write(mzXMLOutFileP, *msd);
		}else
			Rcpp::Rcerr << format << " format not supported! Please try mgf, mzML, mzXML or mz5." << std::endl;
	}else
		Rcpp::Rcerr << "No pwiz object available! Please open a file first!" << std::endl;
}


string RcppPwiz::getFilename (  )
{

    return filename;
}

int RcppPwiz::getLastScan() const
{
    if (msd != NULL)
    {
        SpectrumListPtr slp = msd->run.spectrumListPtr;
        return slp->size();
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
                                     Rcpp::_["sample"]		  = sample[0]->name,
                                     Rcpp::_["source"]        = scansetting[0]->sourceFilePtrs[0]->location
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


Rcpp::List RcppPwiz::getScanHeaderInfo ( int whichScan  )
{
    if (msd != NULL)
    {
        SpectrumListPtr slp = msd->run.spectrumListPtr;
        if ((whichScan <= 0) || (whichScan > slp->size()))
        {
            Rprintf("Index out of bounds [1 ... %d].\n", slp->size());
            return Rcpp::List::create( );
        }

        RAMPAdapter * adapter = new  RAMPAdapter(filename);
        ScanHeaderStruct header;
        adapter->getScanHeader(whichScan - 1, header);

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
    }
    else
    {
        Rprintf("Warning: pwiz not yet initialized.\n ");
        return Rcpp::List::create( );
    }
}

Rcpp::DataFrame RcppPwiz::getAllScanHeaderInfo ( )
{
    if (msd != NULL)
    {
        if (!isInCacheAllScanHeaderInfo)
        {
            SpectrumListPtr slp = msd->run.spectrumListPtr;
            int N = slp->size();

            ScanHeaderStruct scanHeader;
            RAMPAdapter * adapter = new  RAMPAdapter(filename);
            Rcpp::IntegerVector seqNum(N); // number in sequence observed file (1-based)
            Rcpp::IntegerVector acquisitionNum(N); // scan number as declared in File (may be gaps)
            Rcpp::IntegerVector msLevel(N);
            Rcpp::IntegerVector peaksCount(N);
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
            //char scanType[SCANTYPE_LENGTH];
            //char activationMethod[SCANTYPE_LENGTH];
            //char possibleCharges[SCANTYPE_LENGTH];
            //int numPossibleCharges;
            //bool possibleChargesArray[CHARGEARRAY_LENGTH]; /* NOTE: does NOT include "precursorCharge" information; only from "possibleCharges" */
            Rcpp::IntegerVector mergedScan(N);  /* only if MS level > 1 */
            Rcpp::IntegerVector mergedResultScanNum(N); /* scan number of the resultant merged scan */
            Rcpp::IntegerVector mergedResultStartScanNum(N); /* smallest scan number of the scanOrigin for merged scan */
            Rcpp::IntegerVector mergedResultEndScanNum(N); /* largest scan number of the scanOrigin for merged scan */

            for (int whichScan=1; whichScan <= N; whichScan++)
            {
                adapter->getScanHeader(whichScan - 1, scanHeader);
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
                                    //Rcpp::_["scanType"]              	= scanType,
                                    //Rcpp::_["activationMethod"]       = activationMethod,
                                    //Rcpp::_["possibleCharges"]        = possibleCharges,
                                    //Rcpp::_["numPossibleCharges"]     = numPossibleCharges,
                                    //Rcpp::_["possibleChargesArray"]   = possibleChargesArray,
                                    Rcpp::_["mergedScan"]               = mergedScan,
                                    Rcpp::_["mergedResultScanNum"]      = mergedResultScanNum,
                                    Rcpp::_["mergedResultStartScanNum"] = mergedResultStartScanNum,
                                    Rcpp::_["mergedResultEndScanNum"]   = mergedResultEndScanNum
                                );
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

Rcpp::DataFrame RcppPwiz::getChromatogramsInfo()
{
    if (msd != NULL)
    {
        if (!isInCacheChromatogramsInfo)
        {
            SpectrumListPtr slp = msd->run.spectrumListPtr;
            int N = slp->size();
            ScanHeaderStruct scanHeader;
			RAMPAdapter * adapter = new  RAMPAdapter(filename);
			
            Rcpp::IntegerVector msLevel(N);
            Rcpp::NumericVector retentionTime(N);
            Rcpp::NumericVector basePeakIntensity(N);
            Rcpp::NumericVector totIonCurrent(N);

            SpectrumInfo info;
            SpectrumPtr s;
            vector<MZIntensityPair> pairs;

            for (int whichScan=1; whichScan <= N; whichScan++)
            {
				adapter->getScanHeader(whichScan - 1, scanHeader);
                msLevel[whichScan-1] = scanHeader.msLevel;
                retentionTime[whichScan-1] = scanHeader.retentionTime;
                basePeakIntensity[whichScan-1] = scanHeader.basePeakIntensity;
                totIonCurrent[whichScan-1] = scanHeader.totIonCurrent;
            }

            chromatogramsInfo = Rcpp::DataFrame::create(
                                    Rcpp::_["msLevel"] = msLevel,
                                    Rcpp::_["retentionTime"] = retentionTime,
                                    Rcpp::_["basePeakIntensity"] = basePeakIntensity,
                                    Rcpp::_["totIonCurrent"] = totIonCurrent);
            isInCacheChromatogramsInfo = TRUE;
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
