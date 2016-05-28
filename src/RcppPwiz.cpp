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

void RcppPwiz::open(const string& fileName)
{

  filename = fileName;
  msd = new MSDataFile(fileName);

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
    res[i++] = Rcpp::wrap(header.polarity);
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

Rcpp::DataFrame RcppPwiz::getChromatogramsInfo()
{
  if (msd != NULL)
  {
    ChromatogramListPtr clp = msd->run.chromatogramListPtr;
    if(clp.get() == 0)
    {
      Rcpp::Rcerr << "The direct support for chromatogram info is only available in mzML format." << std::endl;
      return Rcpp::DataFrame::create();
    }
    else if(clp->size() == 0)
    {
      Rcpp::Rcerr << "No available chromatogram info." << std::endl;
      return Rcpp::DataFrame::create();
    }
    else
    {
      std::vector<double> time;
      std::vector<double> intensity;
      ChromatogramPtr c = clp->chromatogram(0, true);
      vector<TimeIntensityPair> pairs;
      c->getTimeIntensityPairs (pairs);

      for(int i =0; i < pairs.size(); i++)
      {
        TimeIntensityPair p = pairs.at(i);
        time.push_back(p.time);
        intensity.push_back(p.intensity);
      }

      chromatogramsInfo = Rcpp::DataFrame::create(
          Rcpp::_["time"]	= time,
          Rcpp::_[c->id]	= intensity);

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
