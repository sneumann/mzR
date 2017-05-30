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


// void RcppPwiz::writeMSfile(const string& file, const string& format)
// {
//     if (msd != NULL)
//     {
//         if(format == "mgf")
//         {
//             std::ofstream* mgfOutFileP = new std::ofstream(file.c_str());
//             Serializer_MGF serializerMGF;
//             serializerMGF.write(*mgfOutFileP, *msd);
//             mgfOutFileP->flush();
//             mgfOutFileP->close();
//         }
//         else if(format == "mzxml")
//         {
//             std::ofstream mzXMLOutFileP(file.c_str());
//             Serializer_mzXML::Config config;
//             config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
//             Serializer_mzXML serializerMzXML(config);
//             serializerMzXML.write(mzXMLOutFileP, *msd);
//         }
//         else if(format == "mzml")
//         {
//             std::ofstream mzXMLOutFileP(file.c_str());
//             Serializer_mzML::Config config;
//             config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
//             Serializer_mzML mzmlSerializer(config);
//             mzmlSerializer.write(mzXMLOutFileP, *msd);
//         }
//         else
//             Rcpp::Rcerr << format << " format not supported! Please try mgf, mzML, mzXML or mz5." << std::endl;
//     }
//     else
//         Rcpp::Rcerr << "No pwiz object available! Please open a file first!" << std::endl;
// }


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

        Rcpp::List res(21);
        std::vector<std::string> names;
        int i = 0;

        names.push_back("seqNum");
        res[i++] = Rcpp::wrap(header.seqNum);
        names.push_back("acquisitionNum");
        res[i++] = Rcpp::wrap(header.acquisitionNum);
        names.push_back("msLevel");
        res[i++] = Rcpp::wrap(header.msLevel);

        names.push_back("polarity");
	SpectrumPtr sp = slp->spectrum(whichScan - 1, true); // Is TRUE neccessary here ? 
	CVParam param = sp->cvParamChild(MS_scan_polarity);
	res[i++] = Rcpp::wrap( (param.cvid==MS_negative_scan ? 0 : (param.cvid==MS_positive_scan ? +1 : -1 ) ) );
	
	names.push_back("peaksCount");
        res[i++] = Rcpp::wrap(header.peaksCount);
        names.push_back("totIonCurrent");
        res[i++] = Rcpp::wrap(header.totIonCurrent);
        names.push_back("retentionTime");
        res[i++] = Rcpp::wrap(header.retentionTime);
        names.push_back("basePeakMZ");
        res[i++] = Rcpp::wrap(header.basePeakMZ);
        names.push_back("basePeakIntensity");
        res[i++] = Rcpp::wrap(header.basePeakIntensity);
        names.push_back("collisionEnergy");
        res[i++] = Rcpp::wrap(header.collisionEnergy);
        names.push_back("ionisationEnergy");
        res[i++] = Rcpp::wrap(header.ionisationEnergy);
        names.push_back("lowMZ");
        res[i++] = Rcpp::wrap(header.lowMZ);
        names.push_back("highMZ");
        res[i++] = Rcpp::wrap(header.highMZ);
        names.push_back("precursorScanNum");
        res[i++] = Rcpp::wrap(header.precursorScanNum);
        names.push_back("precursorMZ");
        res[i++] = Rcpp::wrap(header.precursorMZ);
        names.push_back("precursorCharge");
        res[i++] = Rcpp::wrap(header.precursorCharge);
        names.push_back("precursorIntensity");
        res[i++] = Rcpp::wrap(header.precursorIntensity);
        names.push_back("mergedScan");
        res[i++] = Rcpp::wrap(header.mergedScan);
        names.push_back("mergedResultScanNum");
        res[i++] = Rcpp::wrap(header.mergedResultScanNum);
        names.push_back("mergedResultStartScanNum");
        res[i++] = Rcpp::wrap(header.mergedResultStartScanNum);
        names.push_back("mergedResultEndScanNum");
        res[i++] = Rcpp::wrap(header.mergedResultEndScanNum);

	// delete adapter issue #64
	delete adapter;
	adapter = NULL;
        res.attr("names") = names;
        return res;
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
            Rcpp::IntegerVector polarity(N);
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

		SpectrumPtr sp = slp->spectrum(whichScan-1, true); // Is TRUE neccessary here ? 
		CVParam param = sp->cvParamChild(MS_scan_polarity);
		polarity[whichScan-1] = (param.cvid==MS_negative_scan ? 0 : (param.cvid==MS_positive_scan ? +1 : -1 ) );
			
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

            allScanHeaderInfo = header;
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

// copyWriteMSFile copies (general) content from the originating MS file and
// replaces the Spectrum list with the new data provided with arguments
// spctr_header and spctr_data.
// TODO: add strings that describe what processings have been done in R.
// We're copying:
// o fileDescription (adding also the originating MS file)
// o softwareList (adding also additional info)
// o instrumentConfigurationList
// o dataProcessingList (adding also additional info)
// o run: all "general" info from the run.
// Potential additional parameters:
// - centroided: whether spectra data is centroided: NA, TRUE, FALSE.
// - processings: string/vector with all processing steps.
// - software(s): name (and version?) of software.
void RcppPwiz::copyWriteMSfile(const string& file, const string& format,
			       const string& originalFile,
			       Rcpp::DataFrame spctr_header,
			       Rcpp::List spctr_data,
			       bool rtime_seconds) {
  MSDataFile *msd;
  msd = new MSDataFile(originalFile);
  MSData newmsd;
  newmsd.cvs = defaultCVList();

  // Break the header down into its elements/columns:
  Rcpp::IntegerVector seqNum = spctr_header["seqNum"];
  Rcpp::IntegerVector acquisitionNum = spctr_header["acquisitionNum"];
  Rcpp::IntegerVector msLevel = spctr_header["msLevel"];
  Rcpp::IntegerVector polarity = spctr_header["polarity"];
  Rcpp::IntegerVector peaksCount = spctr_header["peaksCount"];
  Rcpp::IntegerVector totIonCurrent = spctr_header["totIonCurrent"];
  Rcpp::NumericVector retentionTime = spctr_header["retentionTime"];
  Rcpp::NumericVector basePeakMZ = spctr_header["basePeakMZ"];
  Rcpp::NumericVector basePeakIntensity = spctr_header["basePeakIntensity"];
  Rcpp::NumericVector collisionEnergy = spctr_header["collisionEnergy"];
  Rcpp::NumericVector ionisationEnergy = spctr_header["ionisationEnergy"];
  Rcpp::NumericVector lowMZ = spctr_header["lowMZ"];
  Rcpp::NumericVector highMZ = spctr_header["highMZ"];
  Rcpp::IntegerVector precursorScanNum = spctr_header["precursorScanNum"];
  Rcpp::NumericVector precursorMZ = spctr_header["precursorMZ"];
  Rcpp::IntegerVector precursorCharge = spctr_header["precursorCharge"];
  Rcpp::NumericVector precursorIntensity = spctr_header["precursorIntensity"];
  Rcpp::NumericVector mergedScan = spctr_header["mergedScan"];
  // Skipping mergedResultScanNum, mergedResultStartScanNum and mergedResultEndScanNum

  // Copy data from the original file.
  // o fileDescription with: fileContent, sourceFileList
  //   TODO: if we did filtering on MS or centroiding we might have to adapt the
  //   fileContent.
  newmsd.fileDescription = msd->fileDescription;
  bool is_ms1 = false;
  bool is_msn = false;
  for (int i = 0; i < msLevel.size(); i++) {
    if (msLevel[i] == 1)
      is_ms1 = true;
    if (msLevel[i] > 1)
      is_msn = true;
  }
  if (is_ms1)
    newmsd.fileDescription.fileContent.set(MS_MS1_spectrum);
  if (is_msn)
    newmsd.fileDescription.fileContent.set(MS_MSn_spectrum);
  // The serializer adds also the original file here AND the newly written file.
  // TODO: check if we need to add processing steps here too.
  // newmsd.fileDescription.sourceFilePtrs.push_back();
  // o softwareList
  newmsd.softwarePtrs = msd->softwarePtrs;
  SoftwarePtr softRPack(new Software);
  softRPack->id = "mzR";
  // softRPack->version = "";
  newmsd.softwarePtrs.push_back(softRPack);
  // o instrumentConfigurationList
  vector<InstrumentConfigurationPtr> icp = msd->instrumentConfigurationPtrs;
  newmsd.instrumentConfigurationPtrs = icp;

  // sampleList
  newmsd.samplePtrs = msd->samplePtrs;

  // newmsd.instrumentConfigurationPtrs.push_back(icp);
  // o dataProcessingList
  newmsd.dataProcessingPtrs = msd->dataProcessingPtrs;
  // Add processing.
  DataProcessingPtr dpR(new DataProcessing);
  dpR->id = "R processing";

  ProcessingMethod procR;
  procR.order = newmsd.dataProcessingPtrs.size() + 1;
  procR.softwarePtr = softRPack;
  // TODO what has been done in R...
  dpR->processingMethods.push_back(procR);

  newmsd.dataProcessingPtrs.push_back(dpR);
  
  // From MSnbase::Spectrum        Column in the header
  // msLevel integer               $msLevel
  // peaksCount integer
  // rt numeric
  // acquisitionNum integer        $acquisitionNum
  // scanIndex integer             $seqNum
  // tic numeric                   $totIonCurrent
  // mz numeric                    peaks()[, 1]
  // intensity numeric             peaks()[, 2]
  // fromFile integer
  // centroided logical
  // smoothed logical
  // polarity integer              $polarity: 0 negative, 1 positive, -1 unknown
  // Spectrum2
  // merged numeric                $mergedScan
  // precScanNum integer           $precursorScanNum
  // precursorMz numeric           $precursorMz
  // precursorIntensity numeric    $precursorIntensity
  // precursorCharge integer       $precursorCharge
  // collisionEnergy numeric       $collisionEnergy
  
  // Initialize the run and fill with data from the original file.
  Run &original_run = msd->run;
  newmsd.run.id = original_run.id;
  newmsd.run.defaultInstrumentConfigurationPtr =
    original_run.defaultInstrumentConfigurationPtr;
  newmsd.run.samplePtr = original_run.samplePtr;
  newmsd.run.startTimeStamp = original_run.startTimeStamp;
  newmsd.run.defaultSourceFilePtr = original_run.defaultSourceFilePtr;

  // Now filling with new data
  addSpectrumList(newmsd, spctr_header, spctr_data, rtime_seconds);
  
  // Now save that one.
  // if(format == "mzml") {
  //   std::ofstream mzXMLOutFileP(file.c_str());
  //   Serializer_mzML::Config config;
  //   config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
  //   Serializer_mzML mzmlSerializer(config);
  //   mzmlSerializer.write(mzXMLOutFileP, newmsd);
  // }
  // else
  //   Rcpp::Rcerr << format << " format not supported! Please try mgf, mzML, mzXML or mz5." << std::endl;
  
  if(format == "mgf")
    {
      std::ofstream* mgfOutFileP = new std::ofstream(file.c_str());
      Serializer_MGF serializerMGF;
      serializerMGF.write(*mgfOutFileP, newmsd);
      mgfOutFileP->flush();
      mgfOutFileP->close();
    }
  else if(format == "mzxml")
    {
      std::ofstream mzXMLOutFileP(file.c_str());
      Serializer_mzXML::Config config;
      config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
      Serializer_mzXML serializerMzXML(config);
      serializerMzXML.write(mzXMLOutFileP, newmsd);
    }
  else if(format == "mzml")
    {
      std::ofstream mzXMLOutFileP(file.c_str());
      Serializer_mzML::Config config;
      config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
      Serializer_mzML mzmlSerializer(config);
      mzmlSerializer.write(mzXMLOutFileP, newmsd);
    }
  else
    Rcpp::Rcerr << format << " format not supported! Please try mgf, mzML, mzXML or mz5." << std::endl;

  // Cleanup.
  delete msd;
}

// writeSpectrumList: writes the provided spectrum data to a file.
void RcppPwiz::writeSpectrumList(const string& file, const string& format,
				 Rcpp::DataFrame spctr_header,
				 Rcpp::List spctr_data,
				 bool rtime_seconds) {
  MSData newmsd;
  newmsd.cvs = defaultCVList();

  // Break the header down into its elements/columns:
  Rcpp::IntegerVector seqNum = spctr_header["seqNum"];
  Rcpp::IntegerVector acquisitionNum = spctr_header["acquisitionNum"];
  Rcpp::IntegerVector msLevel = spctr_header["msLevel"];
  Rcpp::IntegerVector polarity = spctr_header["polarity"];
  Rcpp::IntegerVector peaksCount = spctr_header["peaksCount"];
  Rcpp::IntegerVector totIonCurrent = spctr_header["totIonCurrent"];
  Rcpp::NumericVector retentionTime = spctr_header["retentionTime"];
  Rcpp::NumericVector basePeakMZ = spctr_header["basePeakMZ"];
  Rcpp::NumericVector basePeakIntensity = spctr_header["basePeakIntensity"];
  Rcpp::NumericVector collisionEnergy = spctr_header["collisionEnergy"];
  Rcpp::NumericVector ionisationEnergy = spctr_header["ionisationEnergy"];
  Rcpp::NumericVector lowMZ = spctr_header["lowMZ"];
  Rcpp::NumericVector highMZ = spctr_header["highMZ"];
  Rcpp::IntegerVector precursorScanNum = spctr_header["precursorScanNum"];
  Rcpp::NumericVector precursorMZ = spctr_header["precursorMZ"];
  Rcpp::IntegerVector precursorCharge = spctr_header["precursorCharge"];
  Rcpp::NumericVector precursorIntensity = spctr_header["precursorIntensity"];
  Rcpp::NumericVector mergedScan = spctr_header["mergedScan"];
  // Skipping mergedResultScanNum, mergedResultStartScanNum and mergedResultEndScanNum

  bool is_ms1 = false;
  bool is_msn = false;
  for (int i = 0; i < msLevel.size(); i++) {
    if (msLevel[i] == 1)
      is_ms1 = true;
    if (msLevel[i] > 1)
      is_msn = true;
  }
  if (is_ms1)
    newmsd.fileDescription.fileContent.set(MS_MS1_spectrum);
  if (is_msn)
    newmsd.fileDescription.fileContent.set(MS_MSn_spectrum);
  // The serializer adds also the original file here AND the newly written file.
  // TODO: check if we need to add processing steps here too.
  // newmsd.fileDescription.sourceFilePtrs.push_back();
  // o softwareList
  SoftwarePtr softRPack(new Software);
  softRPack->id = "mzR";
  // softRPack->version = "";
  newmsd.softwarePtrs.push_back(softRPack);
  // Add processing.
  DataProcessingPtr dpR(new DataProcessing);
  dpR->id = "R processing";
  
  ProcessingMethod procR;
  procR.order = newmsd.dataProcessingPtrs.size() + 1;
  procR.softwarePtr = softRPack;
  // TODO what has been done in R...
  dpR->processingMethods.push_back(procR);

  newmsd.dataProcessingPtrs.push_back(dpR);
  
  // From MSnbase::Spectrum        Column in the header
  // msLevel integer               $msLevel
  // peaksCount integer
  // rt numeric
  // acquisitionNum integer        $acquisitionNum
  // scanIndex integer             $seqNum
  // tic numeric                   $totIonCurrent
  // mz numeric                    peaks()[, 1]
  // intensity numeric             peaks()[, 2]
  // fromFile integer
  // centroided logical
  // smoothed logical
  // polarity integer              $polarity: 0 negative, 1 positive, -1 unknown
  // Spectrum2
  // merged numeric                $mergedScan
  // precScanNum integer           $precursorScanNum
  // precursorMz numeric           $precursorMz
  // precursorIntensity numeric    $precursorIntensity
  // precursorCharge integer       $precursorCharge
  // collisionEnergy numeric       $collisionEnergy
  
  // Initialize the run and fill with data from the original file.
  newmsd.run.id = "Experiment 1";
  // newmsd.run.defaultInstrumentConfigurationPtr =
  //   original_run.defaultInstrumentConfigurationPtr;
  // newmsd.run.samplePtr = samplePtr;
  // newmsd.run.startTimeStamp = original_run.startTimeStamp;
  // newmsd.run.defaultSourceFilePtr = original_run.defaultSourceFilePtr;

  // Now filling with new data
  addSpectrumList(newmsd, spctr_header, spctr_data, rtime_seconds);

  if(format == "mgf")
    {
      std::ofstream* mgfOutFileP = new std::ofstream(file.c_str());
      Serializer_MGF serializerMGF;
      serializerMGF.write(*mgfOutFileP, newmsd);
      mgfOutFileP->flush();
      mgfOutFileP->close();
    }
  else if(format == "mzxml")
    {
      std::ofstream mzXMLOutFileP(file.c_str());
      Serializer_mzXML::Config config;
      config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
      Serializer_mzXML serializerMzXML(config);
      serializerMzXML.write(mzXMLOutFileP, newmsd);
    }
  else if(format == "mzml")
    {
      std::ofstream mzXMLOutFileP(file.c_str());
      Serializer_mzML::Config config;
      config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
      Serializer_mzML mzmlSerializer(config);
      mzmlSerializer.write(mzXMLOutFileP, newmsd);
    }
  else
    Rcpp::Rcerr << format << " format not supported! Please try mgf, mzML, mzXML or mz5." << std::endl;
  
}


// Adds information provided in the header and spectra data to the spectrumList
// content of the MSData.
// TODO: BIG QUESTION: what to use as spectrum ID? See issue #105
//       For now: use scan=acquisitionNum[i]. Possible problem: what if the
//       acquisitionNum has gaps? Are MSn spectra still linked correctly to
//       their precursor?
//       Alternative: scan=seqNum[i].
void RcppPwiz::addSpectrumList(MSData& msd,
			       Rcpp::DataFrame& spctr_header,
			       Rcpp::List& spctr_data,
			       bool rtime_seconds) {
  // Break the header down into its elements/columns:
  Rcpp::IntegerVector seqNum = spctr_header["seqNum"];
  Rcpp::IntegerVector acquisitionNum = spctr_header["acquisitionNum"];
  Rcpp::IntegerVector msLevel = spctr_header["msLevel"];
  Rcpp::IntegerVector polarity = spctr_header["polarity"];
  Rcpp::IntegerVector peaksCount = spctr_header["peaksCount"];
  Rcpp::IntegerVector totIonCurrent = spctr_header["totIonCurrent"];
  Rcpp::NumericVector retentionTime = spctr_header["retentionTime"];
  Rcpp::NumericVector basePeakMZ = spctr_header["basePeakMZ"];
  Rcpp::NumericVector basePeakIntensity = spctr_header["basePeakIntensity"];
  Rcpp::NumericVector collisionEnergy = spctr_header["collisionEnergy"];
  Rcpp::NumericVector ionisationEnergy = spctr_header["ionisationEnergy"];
  Rcpp::NumericVector lowMZ = spctr_header["lowMZ"];
  Rcpp::NumericVector highMZ = spctr_header["highMZ"];
  Rcpp::IntegerVector precursorScanNum = spctr_header["precursorScanNum"];
  Rcpp::NumericVector precursorMZ = spctr_header["precursorMZ"];
  Rcpp::IntegerVector precursorCharge = spctr_header["precursorCharge"];
  Rcpp::NumericVector precursorIntensity = spctr_header["precursorIntensity"];
  Rcpp::NumericVector mergedScan = spctr_header["mergedScan"];
  // Skipping mergedResultScanNum, mergedResultStartScanNum and mergedResultEndScanNum
  
  // From MSnbase::Spectrum        Column in the header
  // msLevel integer               $msLevel
  // peaksCount integer
  // rt numeric
  // acquisitionNum integer        $acquisitionNum
  // scanIndex integer             $seqNum
  // tic numeric                   $totIonCurrent
  // mz numeric                    peaks()[, 1]
  // intensity numeric             peaks()[, 2]
  // fromFile integer
  // centroided logical
  // smoothed logical
  // polarity integer              $polarity: 0 negative, 1 positive, -1 unknown
  // Spectrum2
  // merged numeric                $mergedScan
  // precScanNum integer           $precursorScanNum
  // precursorMz numeric           $precursorMz
  // precursorIntensity numeric    $precursorIntensity
  // precursorCharge integer       $precursorCharge
  // collisionEnergy numeric       $collisionEnergy
  
  // Now filling with new data
  shared_ptr<SpectrumListSimple> spectrumList(new SpectrumListSimple);
  msd.run.spectrumListPtr = spectrumList;
  // TODO add also eventual processings.
  for (int i = 0; i < spctr_data.size(); i++) {
    spectrumList->spectra.push_back(SpectrumPtr(new Spectrum));
    Spectrum& spct = *spectrumList->spectra[i];
    spct.set(MS_ms_level, msLevel[i]);
    // [X] polarity
    if (polarity[i] == 0)
      spct.set(MS_negative_scan);
    if (polarity[i] == 1)
      spct.set(MS_positive_scan);
    if (msLevel[i] == 1)
      spct.set(MS_MS1_spectrum);
    else
      spct.set(MS_MSn_spectrum);
    spct.set(MS_lowest_observed_m_z, lowMZ[i]);
    spct.set(MS_highest_observed_m_z, highMZ[i]);
    spct.set(MS_base_peak_m_z, basePeakMZ[i]);
    spct.set(MS_base_peak_intensity, basePeakIntensity[i]);
    spct.set(MS_total_ion_current, totIonCurrent[i]);
    // TODO:
    // [X] seqNum: number observed in file.
    spct.index = seqNum[i] - 1;
    // [X] acquisitionNum: number as reported (there might be gaps).
    spct.id = "scan=" + boost::lexical_cast<std::string>(acquisitionNum[i]);
    // [ ] peaksCount: no need to set this?
    // [X] retentionTime
    spct.scanList.scans.push_back(Scan());
    spct.scanList.set(MS_no_combination);
    Scan &spct_scan = spct.scanList.scans.back();
    if (rtime_seconds)
      spct_scan.set(MS_scan_start_time, retentionTime[i], UO_second);
    else
      spct_scan.set(MS_scan_start_time, retentionTime[i], UO_minute);
    // MSn - precursor:
    if (precursorScanNum[i] > 0 | precursorMZ[i] > 0) {
      spct.precursors.resize(1);
      Precursor& prec = spct.precursors.front();
      // assume we're linked to acquisitionNum (issue #105)
      prec.spectrumID =
	"scan=" + boost::lexical_cast<std::string>(precursorScanNum[i]);
      prec.activation.set(MS_collision_energy, collisionEnergy[i],
			  UO_electronvolt);
      prec.selectedIons.resize(1);
      prec.selectedIons[0].set(MS_selected_ion_m_z, precursorMZ[i], MS_m_z);
      prec.selectedIons[0].set(MS_peak_intensity, precursorIntensity[i],
			       MS_number_of_detector_counts);
      prec.selectedIons[0].set(MS_charge_state, precursorCharge[i]);
    }
    // [X] collisionEnergy
    // [ ] ionisationEnergy
    // [X] precursorScanNum
    // [X] precursorMZ
    // [X] precursorCharge
    // [X] precursorIntensity
    // [ ] mergedScan
    
    Rcpp::NumericMatrix spct_vals = spctr_data[i];
    // mz values
    Rcpp::NumericVector mz_vals = spct_vals( Rcpp::_, 0);
    BinaryDataArrayPtr spct_mz(new BinaryDataArray);
    spct_mz->set(MS_m_z_array, "", MS_m_z);
    spct_mz->data.resize(mz_vals.size());
    for (int j = 0; j < mz_vals.size(); j++)
      spct_mz->data[j] = mz_vals[j];
    spct.binaryDataArrayPtrs.push_back(spct_mz);
    // intensity values
    Rcpp::NumericVector ints_vals = spct_vals( Rcpp::_, 1);
    BinaryDataArrayPtr spct_ints(new BinaryDataArray);
    spct_ints->set(MS_intensity_array, "", MS_number_of_detector_counts);
    spct_ints->data.resize(ints_vals.size());
    for (int j = 0; j < ints_vals.size(); j++)
      spct_ints->data[j] = ints_vals[j];
    spct.binaryDataArrayPtrs.push_back(spct_ints);
    spct.defaultArrayLength = spct_mz->data.size();
  }
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

// How could we provide writing data? Have a look at src/pwiz/examples.cpp for
// details on how to generate data. But how could we provide/pass that from
// another package?
// Option a: read the original file, read it's content and write that to another
// file overwriting selected spectra.
