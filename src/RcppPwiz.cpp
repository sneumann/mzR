#include "RcppPwiz.h"

RcppPwiz::RcppPwiz()
{
  msd = NULL;
  nativeIdFormat = CVID_Unknown;
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

// void RcppPwiz::open(const string& fileName)
void RcppPwiz::open(Rcpp::StringVector fileName)
{

  filename = Rcpp::as<std::string>(fileName(0));
  msd = new MSDataFile(filename);
  // Better not to guess the native ID format. For mzML/mzXML all should be fine
  // with the default one.
  // nativeIdFormat = id::getDefaultNativeIDFormat(*msd);
}

/* Release all memory on close. */
void RcppPwiz::close()
{
  if (msd != NULL)
    {
      delete msd;
      msd = NULL;
      nativeIdFormat = CVID_Unknown;
      instrumentInfo = Rcpp::List::create();
      chromatogramsInfo = Rcpp::DataFrame::create();
      isInCacheInstrumentInfo = FALSE;
      allScanHeaderInfo = Rcpp::List::create();
      isInCacheAllScanHeaderInfo = FALSE;
    }
}

string RcppPwiz::getFilename() {
  return filename;
}

int RcppPwiz::getLastScan() const {
  if (msd != NULL) {
    SpectrumListPtr slp = msd->run.spectrumListPtr;
    return slp->size();
  }
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return -1;
}

int RcppPwiz::getLastChrom() const {
  if (msd != NULL) {
    ChromatogramListPtr clp = msd->run.chromatogramListPtr;
    return clp->size();
  }
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
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
	      std::string ionisation = "";
	      std::string analyzer = "";
	      std::string detector = "";
	      // Fix issue #113
	      // if (icp[0]->componentList.size() > 0)
	      // That does still not mean we have a ionisation available.
	      // Could be that either analyzer or detector or ionisation is
	      // defined.
	      // Have to use try-catch
	      try {
		ionisation = std::string(adapter.ionisation());		  
	      } catch(...) {}
	      try {
		analyzer = std::string(adapter.analyzer());		  
	      } catch(...) {}
	      try {
		detector = std::string(adapter.detector());		  
	      } catch(...) {}
	      instrumentInfo = Rcpp::List::create(
						  Rcpp::_["manufacturer"]  = std::string(adapter.manufacturer()),
						  Rcpp::_["model"]         = std::string(adapter.model()),
						  Rcpp::_["ionisation"]    = ionisation,
						  Rcpp::_["analyzer"]      = analyzer,
						  Rcpp::_["detector"]      = detector,
						  Rcpp::_["software"]      = (sp.size()>0?sp[0]->id + " " + sp[0]->version:""),
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
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return instrumentInfo;
}

// int RcppPwiz::getAcquisitionNumber(size_t index) const
// {
//   const SpectrumIdentity& si = msd->run.spectrumListPtr->spectrumIdentity(index);
//   string scanNumber = id::translateNativeIDToScanNumber(nativeIdFormat, si.id);
//   if (scanNumber.empty()) {
//     return static_cast<int>(index) + 1;
//   }
//   else
//     return lexical_cast<int>(scanNumber);
//   // return static_cast<int>(index) + 1;
// }

// Using this function instead of pwiz translateNativeIDToScanNumber because
// it randomly causes segfaults on macOS.
int RcppPwiz::getAcquisitionNumber(string id, size_t index) const
{
  if (id.find("controllerType") != std::string::npos) {
    if (id.find("controllerType=0 controllerNumber=1") == std::string::npos)
      return static_cast<int>(index) + 1;
  }
  string e;
  std::smatch match;
  if (id.find("scan=") != std::string::npos)
    e ="scan=(\\d+)";
  else if (id.find("index=") != std::string::npos)
    e = "index=(\\d+)";
  else if (id.find("spectrum=") != std::string::npos)
    e = "spectrum=(\\d+)";
  else if (id.find("scanId=") != std::string::npos)
    e = "scanId=(\\d+)";
  else return static_cast<int>(index) + 1;
  if (std::regex_search(id, match, std::regex(e)))
    return lexical_cast<int>(match[1]);
  else return static_cast<int>(index) + 1;
}

Rcpp::DataFrame RcppPwiz::getScanHeaderInfo (Rcpp::IntegerVector whichScan) {
  if (msd != NULL) {
    SpectrumListPtr slp = msd->run.spectrumListPtr;
    size_t N = slp->size();
    size_t N_scans = whichScan.size();
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
    Rcpp::IntegerVector mergedScan(N_scans);  /* only if MS level > 1 */
    Rcpp::IntegerVector mergedResultScanNum(N_scans); /* scan number of the resultant merged scan */
    Rcpp::IntegerVector mergedResultStartScanNum(N_scans); /* smallest scan number of the scanOrigin for merged scan */
    Rcpp::IntegerVector mergedResultEndScanNum(N_scans); /* largest scan number of the scanOrigin for merged scan */
    Rcpp::NumericVector ionInjectionTime(N_scans); /* The time spent filling an ion trapping device*/
    Rcpp::StringVector filterString(N_scans);
    Rcpp::StringVector spectrumId(N_scans);
    Rcpp::LogicalVector centroided(N_scans);
    Rcpp::NumericVector ionMobilityDriftTime(N_scans);
    Rcpp::NumericVector isolationWindowTargetMZ(N_scans);
    Rcpp::NumericVector isolationWindowLowerOffset(N_scans);
    Rcpp::NumericVector isolationWindowUpperOffset(N_scans);
    Rcpp::NumericVector scanWindowLowerLimit(N_scans);
    Rcpp::NumericVector scanWindowUpperLimit(N_scans);
    
    for (size_t i = 0; i < N_scans; i++) {
      int current_scan = whichScan[i];
      size_t current_index = static_cast<size_t>(current_scan - 1);
      // SpectrumPtr sp = slp->spectrum(current_index, false);
      SpectrumPtr sp = slp->spectrum(current_index, DetailLevel_FullMetadata);
      Scan dummy;
      Scan& scan = sp->scanList.scans.empty() ? dummy : sp->scanList.scans[0];
      if (scan.empty())
	Rprintf("Scan with index %d empty.\n", current_scan);
      // seqNum
      seqNum[i] = current_scan;
      acquisitionNum[i] = getAcquisitionNumber(sp->id, current_index);
      // spectrumId
      spectrumId[i] = Rcpp::String(sp->id);
      // msLevel
      msLevel[i] = sp->cvParam(MS_ms_level).valueAs<int>();
      // peaksCount
      peaksCount[i] = static_cast<int>(sp->defaultArrayLength);
      // totIonCurrent
      totIonCurrent[i] = sp->cvParam(MS_total_ion_current).valueAs<double>();
      // basePeakMZ
      basePeakMZ[i] = sp->cvParam(MS_base_peak_m_z).valueAs<double>();
      // basePeakIntensity
      basePeakIntensity[i] = sp->cvParam(MS_base_peak_intensity).valueAs<double>();
      // ionisationEnerty
      ionisationEnergy[i] = sp->cvParam(MS_ionization_energy_OBSOLETE).valueAs<double>();
      // lowMZ
      lowMZ[i] = sp->cvParam(MS_lowest_observed_m_z).valueAs<double>();
      // highMZ
      highMZ[i] = sp->cvParam(MS_highest_observed_m_z).valueAs<double>();
      // polarity
      CVParam param = sp->cvParamChild(MS_scan_polarity);
      polarity[i] = (param.cvid==MS_negative_scan ? 0 : (param.cvid==MS_positive_scan ? +1 : -1 ) );
      // centroided
      param = sp->cvParamChild(MS_spectrum_representation);
      centroided[i] = (param.cvid==MS_centroid_spectrum ? TRUE : (param.cvid==MS_profile_spectrum ? FALSE : NA_LOGICAL));      
      // retentionTime
      retentionTime[i] = scan.cvParam(MS_scan_start_time).timeInSeconds();
      // ionInjectionTime
      ionInjectionTime[i] = (scan.cvParam(MS_ion_injection_time).timeInSeconds() * 1000);
      // filterString
      filterString[i] = scan.cvParam(MS_filter_string).value.empty() ? NA_STRING : Rcpp::String(scan.cvParam(MS_filter_string).value);
      // ionMobilityDriftTime
      ionMobilityDriftTime[i] = scan.cvParam(MS_ion_mobility_drift_time).value.empty() ? NA_REAL : (scan.cvParam(MS_ion_mobility_drift_time).timeInSeconds() * 1000);
      // scanWindowLowerLimit and scanWindowUpperLimit
      if (!scan.scanWindows.empty()) {
	scanWindowLowerLimit[i] = scan.scanWindows[0].cvParam(MS_scan_window_lower_limit).valueAs<double>();
	scanWindowUpperLimit[i] = scan.scanWindows[0].cvParam(MS_scan_window_upper_limit).valueAs<double>();
      } else {
	scanWindowLowerLimit[i] = NA_REAL;
	scanWindowUpperLimit[i] = NA_REAL;
      }
      // mergedScan - also not supported by RAMPAdapter
      mergedScan[i] = NA_INTEGER;
      mergedResultScanNum[i] = NA_INTEGER;
      mergedResultStartScanNum[i] = NA_INTEGER;
      mergedResultEndScanNum[i] = NA_INTEGER;
      if (!sp->precursors.empty()) {
	const Precursor& precursor = sp->precursors[0];
	// collisionEnergy
	collisionEnergy[i] = precursor.activation.cvParam(MS_collision_energy).valueAs<double>();
	// precursorScanNum
	size_t precursorIndex = slp->find(precursor.spectrumID);
	if (precursorIndex < N) {
	  precursorScanNum[i] = getAcquisitionNumber(precursor.spectrumID, precursorIndex);
	} else {
	  precursorScanNum[i] = NA_INTEGER;
	}
	// precursorMZ, precursorCharge, precursorIntensity
	if (!precursor.selectedIons.empty()) {
	  precursorMZ[i] = precursor.selectedIons[0].cvParam(MS_selected_ion_m_z).value.empty() ? precursor.selectedIons[0].cvParam(MS_m_z).valueAs<double>() : precursor.selectedIons[0].cvParam(MS_selected_ion_m_z).valueAs<double>();
	  precursorCharge[i] = precursor.selectedIons[0].cvParam(MS_charge_state).valueAs<int>();
	  precursorIntensity[i] = precursor.selectedIons[0].cvParam(MS_peak_intensity).valueAs<double>();
	}
	// isolationWindowTargetMZ, ...
	IsolationWindow iwin = sp->precursors[0].isolationWindow;
	if (!iwin.empty()) {
	  isolationWindowTargetMZ[i] = iwin.cvParam(MS_isolation_window_target_m_z).value.empty() ? NA_REAL : iwin.cvParam(MS_isolation_window_target_m_z).valueAs<double>();
	  isolationWindowLowerOffset[i] = iwin.cvParam(MS_isolation_window_lower_offset).value.empty() ? NA_REAL : iwin.cvParam(MS_isolation_window_lower_offset).valueAs<double>();
	  isolationWindowUpperOffset[i] = iwin.cvParam(MS_isolation_window_upper_offset).value.empty() ? NA_REAL : iwin.cvParam(MS_isolation_window_upper_offset).valueAs<double>();
	} else {
	  isolationWindowTargetMZ[i] = NA_REAL;
	  isolationWindowLowerOffset[i] = NA_REAL;
	  isolationWindowUpperOffset[i] = NA_REAL;
	}
      } else {
	collisionEnergy[i] = NA_REAL;
	precursorScanNum[i] = NA_INTEGER;
	precursorMZ[i] = NA_REAL;
	precursorCharge[i] = NA_INTEGER;
	precursorIntensity[i] = NA_REAL;
	mergedScan[i] = NA_INTEGER;
	mergedResultScanNum[i] = NA_INTEGER;
	mergedResultStartScanNum[i] = NA_INTEGER;
	mergedResultEndScanNum[i] = NA_INTEGER;
	isolationWindowTargetMZ[i] = NA_REAL;
	isolationWindowLowerOffset[i] = NA_REAL;
	isolationWindowUpperOffset[i] = NA_REAL;
      }
    }
    
    Rcpp::List header(31);
    std::vector<std::string> names;
    size_t i = 0;
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
    names.push_back("injectionTime");
    header[i++] = Rcpp::wrap(ionInjectionTime);
    names.push_back("filterString");
    header[i++] = Rcpp::wrap(filterString);
    names.push_back("spectrumId");
    header[i++] = Rcpp::wrap(spectrumId);
    names.push_back("centroided");
    header[i++] = Rcpp::wrap(centroided);
    names.push_back("ionMobilityDriftTime");
    header[i++] = Rcpp::wrap(ionMobilityDriftTime);      
    names.push_back("isolationWindowTargetMZ");
    header[i++] = Rcpp::wrap(isolationWindowTargetMZ);      
    names.push_back("isolationWindowLowerOffset");
    header[i++] = Rcpp::wrap(isolationWindowLowerOffset);      
    names.push_back("isolationWindowUpperOffset");
    header[i++] = Rcpp::wrap(isolationWindowUpperOffset);      
    names.push_back("scanWindowLowerLimit");
    header[i++] = Rcpp::wrap(scanWindowLowerLimit);      
    names.push_back("scanWindowUpperLimit");
    header[i++] = Rcpp::wrap(scanWindowUpperLimit);      
    header.attr("names") = names;
    
    return header;
  }
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return Rcpp::DataFrame::create( );
}

Rcpp::DataFrame RcppPwiz::getAllScanHeaderInfo ( ) {
  if (msd != NULL) {
    if (!isInCacheAllScanHeaderInfo) {
      SpectrumListPtr slp = msd->run.spectrumListPtr;
      size_t N = slp->size();
      
      allScanHeaderInfo = getScanHeaderInfo(Rcpp::seq(1, N));
      isInCacheAllScanHeaderInfo = TRUE;	    
    }
    return allScanHeaderInfo ;
  }
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return Rcpp::DataFrame::create( );
}

Rcpp::List RcppPwiz::getPeakList(Rcpp::IntegerVector whichScan) {
  if (msd != NULL) {
    SpectrumListPtr slp = msd->run.spectrumListPtr;
    size_t n_scans = slp->size();
    size_t n_want = whichScan.size();
    int current_scan;
    SpectrumPtr sp;
    BinaryDataArrayPtr mzs,ints;
    std::vector<double> data;
    Rcpp::NumericVector data_matrix;
    Rcpp::List res(n_want);
    for (size_t i = 0; i < n_want; i++) {
      current_scan = whichScan[i];
      if (current_scan < 1 || current_scan > n_scans) {
	Rprintf("Index whichScan out of bounds [1 ... %d].\n", n_scans);
	return Rcpp::List::create( );
      }
      size_t current_index = static_cast<size_t>(current_scan - 1);
      // sp = slp->spectrum(current_index, true);
      sp = slp->spectrum(current_index, DetailLevel_FullData);
      mzs = sp->getMZArray();
      ints = sp->getIntensityArray();
      if (!mzs.get() || !ints.get()) {
	Rcpp::NumericMatrix pks(0, 2);
	res[i] = pks;
	continue;
      }
      if (mzs->data.size() != ints->data.size())
	Rcpp::Rcerr << "Sizes of mz and intensity arrays don't match." << std::endl;
      data = mzs->data;
      data.insert(data.end(), ints->data.begin(), ints->data.end());
      data_matrix = Rcpp::wrap(data);
      data_matrix.attr("dim") = Rcpp::Dimension(ints->data.size(), 2);
      res[i] = data_matrix;
    }
    return res;
  }
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return Rcpp::List::create();
}

/**
 * copyWriteMSFile copies (general) content from the originating MS file and
 * replaces the Spectrum list with the new data provided with arguments
 * spctr_header and spctr_data.
 * TODO: add strings that describe what processings have been done in R.
 * We're copying:
 * o fileDescription (adding also the originating MS file)
 * o softwareList (adding also additional info)
 * o instrumentConfigurationList
 * o dataProcessingList (adding also additional info)
 * o run: all "general" info from the run.
 * Potential additional parameters:
 * - centroided: whether spectra data is centroided: NA, TRUE, FALSE.
 * - processings: string/vector with all processing steps.
 * - software(s): name (and version?) of software.
 * software_processing: has to be a list of character vectors.
 **/
void RcppPwiz::copyWriteMSfile(const string& file, const string& format,
			       const string& originalFile,
			       Rcpp::DataFrame spctr_header,
			       Rcpp::List spctr_data,
			       bool rtime_seconds,
			       Rcpp::List software_processing) {
  MSDataFile *msd;
  msd = new MSDataFile(originalFile);
  MSData newmsd;
  newmsd.cvs = defaultCVList();

  Rcpp::IntegerVector msLevel = spctr_header["msLevel"];
  // Copy data from the original file.
  // o fileDescription with: fileContent, sourceFileList
  //   NOTE: don't copy the file description for mzXML export - somehow the
  //   spectra data will then not be written.
  if (format != "mzxml")
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

  // o paramGroupList
  if (format != "mzxml")
    newmsd.paramGroupPtrs = msd->paramGroupPtrs;
  // o sampleList
  newmsd.samplePtrs = msd->samplePtrs;
  // o instrumentConfigurationList
  newmsd.instrumentConfigurationPtrs = msd->instrumentConfigurationPtrs;
  // o softwareList
  newmsd.softwarePtrs = msd->softwarePtrs;
  // o dataProcessingList
  newmsd.dataProcessingPtrs = msd->dataProcessingPtrs;
  // Add new software and processings:
  if (software_processing.size() > 0) {
    for (int sp = 0; sp < software_processing.size(); sp++) {
      addDataProcessing(newmsd, Rcpp::as<Rcpp::StringVector>(software_processing(sp)));
    }
  }
  
  // o run
  // Initialize the run and fill with data from the original file.
  Run &original_run = msd->run;
  newmsd.run.id = original_run.id;
  if (format != "mzxml") {
    newmsd.run.defaultInstrumentConfigurationPtr =
      original_run.defaultInstrumentConfigurationPtr;
    newmsd.run.samplePtr = original_run.samplePtr;
    newmsd.run.startTimeStamp = original_run.startTimeStamp;
    newmsd.run.defaultSourceFilePtr = original_run.defaultSourceFilePtr;
  }
  // Now filling with new data
  addSpectrumList(newmsd, spctr_header, spctr_data, rtime_seconds);
  
  if (format == "mgf") {
    std::ofstream* mgfOutFileP = new std::ofstream(file.c_str());
    Serializer_MGF serializerMGF;
    serializerMGF.write(*mgfOutFileP, newmsd);
    mgfOutFileP->flush();
    mgfOutFileP->close();
  } else if (format == "mzxml") {
    std::ofstream mzXMLOutFileP(file.c_str());
    Serializer_mzXML::Config config;
    config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
    Serializer_mzXML serializerMzXML(config);
    serializerMzXML.write(mzXMLOutFileP, newmsd);
    mzXMLOutFileP.flush();
    mzXMLOutFileP.close();
  } else if (format == "mzml") {
    std::ofstream mzXMLOutFileP(file.c_str());
    Serializer_mzML::Config config;
    config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
    Serializer_mzML mzmlSerializer(config);
    mzmlSerializer.write(mzXMLOutFileP, newmsd);
    mzXMLOutFileP.flush();
    mzXMLOutFileP.close();
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
				 bool rtime_seconds,
				 Rcpp::List software_processing) {
  MSData newmsd;
  newmsd.cvs = defaultCVList();

  Rcpp::IntegerVector msLevel = spctr_header["msLevel"];
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

  // Add software_processing:
  if (software_processing.size() > 0) {
    for (int sp = 0; sp < software_processing.size(); sp++) {
      addDataProcessing(newmsd, Rcpp::as<Rcpp::StringVector>(software_processing(sp)));
    }
  }
  
  newmsd.run.id = "Experiment_1";

  // Now filling with new data
  addSpectrumList(newmsd, spctr_header, spctr_data, rtime_seconds);

  if (format == "mgf") {
    std::ofstream* mgfOutFileP = new std::ofstream(file.c_str());
    Serializer_MGF serializerMGF;
    serializerMGF.write(*mgfOutFileP, newmsd);
    mgfOutFileP->flush();
    mgfOutFileP->close();
  } else if (format == "mzxml") {
    std::ofstream mzXMLOutFileP(file.c_str());
    Serializer_mzXML::Config config;
    config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
    Serializer_mzXML serializerMzXML(config);
    serializerMzXML.write(mzXMLOutFileP, newmsd);
    mzXMLOutFileP.flush();
    mzXMLOutFileP.close();
  } else if (format == "mzml") {
    std::ofstream mzXMLOutFileP(file.c_str());
    Serializer_mzML::Config config;
    config.binaryDataEncoderConfig.compression = BinaryDataEncoder::Compression_Zlib;
    Serializer_mzML mzmlSerializer(config);
    mzmlSerializer.write(mzXMLOutFileP, newmsd);
    mzXMLOutFileP.flush();
    mzXMLOutFileP.close();
  }
  else
    Rcpp::Rcerr << format << " format not supported! Please try mgf, mzML, mzXML or mz5." << std::endl;
}

/*
 * o soft_proc: is supposed to be a character vector of length >= 2:
 *   soft_proc[0]: The software name (required).
 *   soft_proc[1]: The software version (required).
 *   soft_proc[2]: The CV ID of the software. Use "MS:-1" if not known, in 
 *                 which case we are NOT writing the corresponding CV element.
 *   soft_proc[3-length]: CV IDs of the processing steps (optional). 
 */
void RcppPwiz::addDataProcessing(MSData& msd, Rcpp::StringVector soft_proc) {
  SoftwarePtr new_soft(new Software);
  new_soft->id = soft_proc(0);
  new_soft->version = soft_proc(1);
  int soft_proc_size = soft_proc.size();
  if (soft_proc_size > 2) {
    if (soft_proc(2) != "MS:-1") {
      CVTermInfo cv_term = cvTermInfo(soft_proc(2));
      new_soft->set(cv_term.cvid);
    }
  }
  // Order: get the number of already present dataProcessingPtrs and
  // increment
  int order = msd.dataProcessingPtrs.size() + 1;
  DataProcessingPtr data_processing(new DataProcessing);
  std::ostringstream oss;
  oss << soft_proc[0] << "_processing";
  data_processing->id = oss.str();
  ProcessingMethod proc_meth;
  proc_meth.order = order;
  proc_meth.softwarePtr = new_soft;
  if (soft_proc_size > 3) {
    // Got also processing steps.
    for (int i = 3; i < soft_proc_size; i++) {
      CVTermInfo cv_term = cvTermInfo(soft_proc(i));
      proc_meth.set(cv_term.cvid);
    }
  }
  data_processing->processingMethods.push_back(proc_meth);
  msd.softwarePtrs.push_back(new_soft);
  msd.dataProcessingPtrs.push_back(data_processing);
}

/** Adds information provided in the header and spectra data to the spectrumList
 *  content of the MSData.
 *  TODO: OPEN QUESTION: what to use as spectrum ID? See issue #105
 *       For now: use scan=acquisitionNum[i]. According to the code of the
 *       RAMPAdapter.cpp this seems to be correct - the scan number (i.e.
 *       acquisitionNum) is extracted/guessed from the id of the spectrum.
 *       Need to test: what if the acquisitionNum has gaps? Are MSn spectra
 *       still linked correctly to their precursor?
 *       Alternative: scan=seqNum[i].
 **/
void RcppPwiz::addSpectrumList(MSData& msd,
			       Rcpp::DataFrame& spctr_header,
			       Rcpp::List& spctr_data,
			       bool rtime_seconds) {
  int precursor_idx;
  // Break the header down into its elements/columns:
  Rcpp::IntegerVector seqNum = spctr_header["seqNum"];
  Rcpp::IntegerVector acquisitionNum = spctr_header["acquisitionNum"];
  Rcpp::IntegerVector msLevel = spctr_header["msLevel"];
  Rcpp::IntegerVector polarity = spctr_header["polarity"];
  Rcpp::IntegerVector peaksCount = spctr_header["peaksCount"];
  Rcpp::NumericVector totIonCurrent = spctr_header["totIonCurrent"];
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
  Rcpp::IntegerVector mergedScan = spctr_header["mergedScan"];
  // Skipping mergedResultScanNum, mergedResultStartScanNum and mergedResultEndScanNum
  Rcpp::NumericVector ionInjectionTime = spctr_header["injectionTime"];
  Rcpp::StringVector filterString = spctr_header["filterString"];
  Rcpp::StringVector spectrumId = spctr_header["spectrumId"];
  Rcpp::LogicalVector centroided = spctr_header["centroided"];
  Rcpp::NumericVector ionMobilityDriftTime = spctr_header["ionMobilityDriftTime"];
  Rcpp::NumericVector isolationWindowTargetMZ = spctr_header["isolationWindowTargetMZ"];
  Rcpp::NumericVector isolationWindowLowerOffset = spctr_header["isolationWindowLowerOffset"];
  Rcpp::NumericVector isolationWindowUpperOffset = spctr_header["isolationWindowUpperOffset"];
  Rcpp::NumericVector scanWindowLowerLimit = spctr_header["scanWindowLowerLimit"];
  Rcpp::NumericVector scanWindowUpperLimit = spctr_header["scanWindowUpperLimit"];
  
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
  // Add the default Processing pointer (fix issue #151
  spectrumList->dp = msd.dataProcessingPtrs[(msd.dataProcessingPtrs.size() - 1)];
  msd.run.spectrumListPtr = spectrumList;
  // TODO add also eventual processings.
  for (int i = 0; i < spctr_data.size(); i++) {
    spectrumList->spectra.push_back(SpectrumPtr(new Spectrum));
    Spectrum& spct = *spectrumList->spectra[i];
    spct.set(MS_ms_level, msLevel[i]);
    // centroided
    if (centroided[i] != NA_LOGICAL && centroided[i] == TRUE)
      spct.set(MS_centroid_spectrum);
    if (centroided[i] != NA_LOGICAL && centroided[i] == FALSE)
      spct.set(MS_profile_spectrum);
    // [X] polarity
    if (polarity[i] == 0)
      spct.set(MS_negative_scan);
    if (polarity[i] == 1)
      spct.set(MS_positive_scan);
    if (msLevel[i] == 1)
      spct.set(MS_MS1_spectrum);
    else
      spct.set(MS_MSn_spectrum);
    spct.set(MS_lowest_observed_m_z, lowMZ[i], MS_m_z);
    spct.set(MS_highest_observed_m_z, highMZ[i], MS_m_z);
    spct.set(MS_base_peak_m_z, basePeakMZ[i], MS_m_z);
    spct.set(MS_base_peak_intensity, basePeakIntensity[i],
	     MS_number_of_detector_counts);
    spct.set(MS_total_ion_current, totIonCurrent[i]);
    // TODO:
    // [X] seqNum: number observed in file.
    spct.index = seqNum[i] - 1;	// Or just i?
    // [X] acquisitionNum: number as reported (there might be gaps).
    // spct.id = "scan=" + boost::lexical_cast<std::string>(acquisitionNum[i]);
    spct.id = spectrumId[i];	// Use the provided ID instead
    // [ ] peaksCount: no need to set this?
    // [X] retentionTime
    spct.scanList.scans.push_back(Scan());
    spct.scanList.set(MS_no_combination);
    Scan &spct_scan = spct.scanList.scans.back();
    if (rtime_seconds)
      spct_scan.set(MS_scan_start_time, retentionTime[i], UO_second);
    else
      spct_scan.set(MS_scan_start_time, retentionTime[i], UO_minute);
    if (ionInjectionTime[i] > 0)
      spct_scan.set(MS_ion_injection_time, ionInjectionTime[i], UO_millisecond);

    if (!Rcpp::StringVector::is_na(filterString[i]))
      spct_scan.set(MS_filter_string, filterString[i]);

    if (ionMobilityDriftTime[i] != NA_REAL)
      spct_scan.set(MS_ion_mobility_drift_time, ionMobilityDriftTime[i],
		    UO_millisecond);

    // scanWindow
    if (scanWindowLowerLimit[i] != NA_REAL && scanWindowUpperLimit[i] != NA_REAL) {
      spct_scan.scanWindows.push_back(ScanWindow(scanWindowLowerLimit[i], scanWindowUpperLimit[i], MS_m_z));
    }
    // MSn - precursor:
    if (precursorScanNum[i] > 0 | precursorMZ[i] > 0) {
      // Fill precursor data. This preserves the precursor data even if the
      // precursor scan is not available (e.g. after MS level filtering).
      spct.precursors.resize(1);
      Precursor& prec = spct.precursors.front();
      if (collisionEnergy[i] != 0) {
	prec.activation.set(MS_collision_induced_dissociation);
	prec.activation.set(MS_collision_energy, collisionEnergy[i],
			    UO_electronvolt);
      }
      prec.selectedIons.resize(1);
      prec.selectedIons[0].set(MS_selected_ion_m_z, precursorMZ[i], MS_m_z);
      prec.selectedIons[0].set(MS_peak_intensity, precursorIntensity[i],
			       MS_number_of_detector_counts);
      prec.selectedIons[0].set(MS_charge_state, precursorCharge[i]);
      // Get the spectrumId of the precursor. Assuming that precursorScanNum is
      // linked to the acquisitionNum of the precursor.
      // This seems to be correct, since both the acquisitionNum and the
      // precursorNum are extracted from the respective spectrum's ID.
      precursor_idx = -1;
      for (int j = 0; j < spctr_data.size(); j++) {
	if (precursorScanNum[i] == acquisitionNum[j]) {
	  precursor_idx = j;
	  break;
	}
      }
      if (precursor_idx >= 0) {
	prec.spectrumID = spectrumId[precursor_idx];
      }
      // isolation window
      if (isolationWindowTargetMZ[i] != NA_REAL) {
	prec.isolationWindow.set(MS_isolation_window_target_m_z, isolationWindowTargetMZ[i]);
	prec.isolationWindow.set(MS_isolation_window_lower_offset, isolationWindowLowerOffset[i]);
	prec.isolationWindow.set(MS_isolation_window_upper_offset, isolationWindowUpperOffset[i]);
      }
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
      Rf_warningcall(R_NilValue, "The direct support for chromatogram info is only available in mzML format.");
      return Rcpp::DataFrame::create();
    } else if (clp->size() == 0) {
      Rf_warningcall(R_NilValue, "No available chromatogram info.");
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
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return Rcpp::DataFrame::create( );
}

// get the header info for chromatograms.
Rcpp::DataFrame RcppPwiz::getChromatogramHeaderInfo (Rcpp::IntegerVector whichChrom)
{
  if (msd != NULL) {
    // CVID nativeIdFormat_ = id::getDefaultNativeIDFormat(*msd);
    ChromatogramListPtr clp = msd->run.chromatogramListPtr;
    if (clp.get() == 0) {
      Rf_warningcall(R_NilValue, "The direct support for chromatogram info is only available in mzML format.");
      return Rcpp::DataFrame::create();
    } else if (clp->size() == 0) {
      Rf_warningcall(R_NilValue, "No available chromatogram info.");
      return Rcpp::DataFrame::create();
    }

    int N = clp->size();  
    int N_chrom = whichChrom.size();

    Rcpp::StringVector chromatogramId(N_chrom); // the ID from the chrom
    Rcpp::IntegerVector chromatogramIndex(N_chrom);  // In contrast to the acquisitionNum we report here the index (1 based) of the chromatogram within the file.
    Rcpp::IntegerVector polarity(N_chrom);
    // MS:1000827: isolation window target m/z
    // MS:1000828: isolation window lower offset
    // MS:1000829: isolation window upper offset
    Rcpp::NumericVector precursorIsolationWindowTargetMZ(N_chrom);
    Rcpp::NumericVector precursorIsolationWindowLowerOffset(N_chrom);
    Rcpp::NumericVector precursorIsolationWindowUpperOffset(N_chrom);
    Rcpp::NumericVector precursorCollisionEnergy(N_chrom);
    Rcpp::NumericVector productIsolationWindowTargetMZ(N_chrom);
    Rcpp::NumericVector productIsolationWindowLowerOffset(N_chrom);
    Rcpp::NumericVector productIsolationWindowUpperOffset(N_chrom);
        
    for (int i = 0; i < N_chrom; i++) {
      int current_chrom = whichChrom[i];
      if (current_chrom < 0 || current_chrom > N) {
	Rf_warningcall(R_NilValue, "Provided index out of bounds.");
	Rcpp::Rcerr << "Provided index out of bounds" << std::endl;
      }
      ChromatogramPtr ch = clp->chromatogram(current_chrom - 1, false);
      chromatogramId[i] = ch->id;
      chromatogramIndex[i] = current_chrom;
      CVParam param = ch->cvParamChild(MS_scan_polarity);
      polarity[i] = (param.cvid==MS_negative_scan ? 0 : (param.cvid==MS_positive_scan ? +1 : -1 ) );
      if (!ch->precursor.empty()) {
	precursorIsolationWindowTargetMZ[i] = ch->precursor.isolationWindow.cvParam(MS_isolation_window_target_m_z).value.empty() ? NA_REAL : ch->precursor.isolationWindow.cvParam(MS_isolation_window_target_m_z).valueAs<double>();
	precursorIsolationWindowLowerOffset[i] = ch->precursor.isolationWindow.cvParam(MS_isolation_window_lower_offset).value.empty() ? NA_REAL : ch->precursor.isolationWindow.cvParam(MS_isolation_window_lower_offset).valueAs<double>();
	precursorIsolationWindowUpperOffset[i] = ch->precursor.isolationWindow.cvParam(MS_isolation_window_upper_offset).value.empty() ? NA_REAL : ch->precursor.isolationWindow.cvParam(MS_isolation_window_upper_offset).valueAs<double>();
	precursorCollisionEnergy[i] = ch->precursor.activation.cvParam(MS_collision_energy).value.empty() ? NA_REAL : ch->precursor.activation.cvParam(MS_collision_energy).valueAs<double>(); 
      } else {
	precursorIsolationWindowTargetMZ[i] = NA_REAL;
	precursorIsolationWindowLowerOffset[i] = NA_REAL;
	precursorIsolationWindowUpperOffset[i] = NA_REAL;
	precursorCollisionEnergy[i] = NA_REAL;
      }
      if (!ch->product.empty()) {
	productIsolationWindowTargetMZ[i] = ch->product.isolationWindow.cvParam(MS_isolation_window_target_m_z).value.empty() ? NA_REAL : ch->product.isolationWindow.cvParam(MS_isolation_window_target_m_z).valueAs<double>();
	productIsolationWindowLowerOffset[i] = ch->product.isolationWindow.cvParam(MS_isolation_window_lower_offset).value.empty() ? NA_REAL : ch->product.isolationWindow.cvParam(MS_isolation_window_lower_offset).valueAs<double>();
	productIsolationWindowUpperOffset[i] = ch->product.isolationWindow.cvParam(MS_isolation_window_upper_offset).value.empty() ? NA_REAL : ch->product.isolationWindow.cvParam(MS_isolation_window_upper_offset).valueAs<double>();
      } else {
	productIsolationWindowTargetMZ[i] = NA_REAL;
	productIsolationWindowLowerOffset[i] = NA_REAL;
	productIsolationWindowUpperOffset[i] = NA_REAL;
      }
    }
    Rcpp::List chromHeader(10);
    std::vector<std::string> names;
    int i = 0;
    names.push_back("chromatogramId");
    chromHeader[i++] = Rcpp::wrap(chromatogramId);
    names.push_back("chromatogramIndex");
    chromHeader[i++] = Rcpp::wrap(chromatogramIndex);
    names.push_back("polarity");
    chromHeader[i++] = Rcpp::wrap(polarity);
    names.push_back("precursorIsolationWindowTargetMZ");
    chromHeader[i++] = Rcpp::wrap(precursorIsolationWindowTargetMZ);
    names.push_back("precursorIsolationWindowLowerOffset");
    chromHeader[i++] = Rcpp::wrap(precursorIsolationWindowLowerOffset);
    names.push_back("precursorIsolationWindowUpperOffset");
    chromHeader[i++] = Rcpp::wrap(precursorIsolationWindowUpperOffset);
    names.push_back("precursorCollisionEnergy");
    chromHeader[i++] = Rcpp::wrap(precursorCollisionEnergy);
    names.push_back("productIsolationWindowTargetMZ");
    chromHeader[i++] = Rcpp::wrap(productIsolationWindowTargetMZ);
    names.push_back("productIsolationWindowLowerOffset");
    chromHeader[i++] = Rcpp::wrap(productIsolationWindowLowerOffset);
    names.push_back("productIsolationWindowUpperOffset");
    chromHeader[i++] = Rcpp::wrap(productIsolationWindowUpperOffset);
    
    chromHeader.attr("names") = names;
    return chromHeader;
  }
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return Rcpp::DataFrame::create( );
}

Rcpp::DataFrame RcppPwiz::getAllChromatogramHeaderInfo ( ) {
  if (msd != NULL) {
    ChromatogramListPtr clp = msd->run.chromatogramListPtr;
    if (clp.get() == 0) {
      Rf_warningcall(R_NilValue, "The direct support for chromatogram info is only available in mzML format.");
      return Rcpp::DataFrame::create();
    }
    int N = clp->size();
    if (N > 0) {
      return getChromatogramHeaderInfo(Rcpp::seq(1, N));
    } else {
      Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
    }
  }
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
      //Rprintf("%d\n",1);
      for (int i = 0; i < scanNumbers.size(); i++)
        {
	  SpectrumPtr s = slp->spectrum(scanNumbers[i] - 1, DetailLevel_FullMetadata);
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

  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return Rcpp::NumericMatrix(0,0);
}

string RcppPwiz::getRunStartTimeStamp() {
  if (msd != NULL) {
    return msd->run.startTimeStamp;
  }
  Rf_warningcall(R_NilValue, "pwiz not yet initialized.");
  return "";
}
