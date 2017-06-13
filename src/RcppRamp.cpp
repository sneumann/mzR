#include "RcppRamp.h"

RcppRamp::RcppRamp()
{
    ramp = NULL;
    runInfo = Rcpp::List::create( );
    isInCacheRunInfo = FALSE;
    instrumentInfo = Rcpp::List::create( );
    isInCacheInstrumentInfo = FALSE;
    allScanHeaderInfo = Rcpp::DataFrame::create( );
    isInCacheAllScanHeaderInfo = FALSE;
    filename = Rcpp::StringVector::create( );
}

RcppRamp::~RcppRamp()
{
    RcppRamp::close();
}

void RcppRamp::open( const char* fileName, bool declaredScansOnly )
{
    RcppRamp::close();
    ramp = new cRamp(fileName, declaredScansOnly);
    if (ramp->OK())
    {
        filename = Rcpp::StringVector::create( fileName );
    }
    else
    {
        RcppRamp::close();
        Rprintf("Failed to open file.\n ");
    }
}

void RcppRamp::close()
{
    if (ramp != NULL)
    {
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


Rcpp::StringVector RcppRamp::getFilename (  )
{
    if (ramp != NULL)
    {
        return filename;
    }
    Rprintf("Warning: Ramp not yet initialized.\n ");
    return filename;
}

Rcpp::List RcppRamp::getRunInfo (  )
{
    if (ramp != NULL)
    {
        if (!isInCacheRunInfo)
        {
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
        }
        else
        {
            // Rprintf("Read from cache.\n ");
        }
        return runInfo;
    }
    Rprintf("Warning: Ramp not yet initialized.\n");
    return runInfo;
}

Rcpp::List RcppRamp::getInstrumentInfo ( )
{
    if (ramp != NULL)
    {
        if (!isInCacheInstrumentInfo)
        {
            // Rprintf("Read from disk.\n ");
            rampInstrumentInfo *info = ramp->getInstrumentInfo(); // NULL for mzData

            if (info != NULL)
            {
                InstrumentStruct * data = info->m_instrumentStructPtr;

                instrumentInfo = Rcpp::List::create(
                                     Rcpp::_["manufacturer"]  = std::string(data->manufacturer),
                                     Rcpp::_["model"]         = std::string(data->model),
                                     Rcpp::_["ionisation"]    = std::string(data->ionisation),
                                     Rcpp::_["analyzer"]      = std::string(data->analyzer),
                                     Rcpp::_["detector"]      = std::string(data->detector)
                                 ) ;
                delete info;
            }
            else
            {
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
        else
        {
            // Rprintf("Read from cache.\n ");
        }
        return(instrumentInfo);
    }
    Rprintf("Warning: Ramp not yet initialized.\n ");
    return instrumentInfo;
}

Rcpp::List RcppRamp::getScanHeaderInfo ( int whichScan  )
{
    if (ramp != NULL)
    {
        if ((whichScan <= 0) || (whichScan > ramp->getLastScan()))
        {
            Rprintf("Index whichScan out of bounds [1 ... %d].\n", ramp->getLastScan());
            return Rcpp::List::create( );
        }
        rampScanInfo *info = ramp->getScanHeaderInfo( whichScan );
        ScanHeaderStruct data = info->m_data;
        delete info;

        std::vector<std::string> names;
        Rcpp::List header(22);
        int i = 0;

        names.push_back("seqNum");
        header[i++] = Rcpp::wrap(data.seqNum);
        names.push_back("acquisitionNum");
        header[i++] = Rcpp::wrap(data.acquisitionNum);
        names.push_back("msLevel");
        header[i++] = Rcpp::wrap(data.msLevel);
        names.push_back("polarity");
        header[i++] = Rcpp::wrap(data.polarity);
        names.push_back("peaksCount");
        header[i++] = Rcpp::wrap(data.peaksCount);
        names.push_back("totIonCurrent");
        header[i++] = Rcpp::wrap(data.totIonCurrent);
        names.push_back("retentionTime");
        header[i++] = Rcpp::wrap(data.retentionTime);
        names.push_back("basePeakMZ");
        header[i++] = Rcpp::wrap(data.basePeakMZ);
        names.push_back("basePeakIntensity");
        header[i++] = Rcpp::wrap(data.basePeakIntensity);
        names.push_back("collisionEnergy");
        header[i++] = Rcpp::wrap(data.collisionEnergy);
        names.push_back("ionisationEnergy");
        header[i++] = Rcpp::wrap(data.ionisationEnergy);
        names.push_back("lowMZ");
        header[i++] = Rcpp::wrap(data.lowMZ);
        names.push_back("highMZ");
        header[i++] = Rcpp::wrap(data.highMZ);
        names.push_back("precursorScanNum");
        header[i++] = Rcpp::wrap(data.precursorScanNum);
        names.push_back("precursorMZ");
        header[i++] = Rcpp::wrap(data.precursorMZ);
        names.push_back("precursorCharge");
        header[i++] = Rcpp::wrap(data.precursorCharge);
        names.push_back("precursorIntensity");
        header[i++] = Rcpp::wrap(data.precursorIntensity);
        names.push_back("mergedScan");
        header[i++] = Rcpp::wrap(data.mergedScan);
        names.push_back("mergedResultScanNum");
        header[i++] = Rcpp::wrap(data.mergedResultScanNum);
        names.push_back("mergedResultStartScanNum");
        header[i++] = Rcpp::wrap(data.mergedResultStartScanNum);
        names.push_back("mergedResultEndScanNum");
        header[i++] = Rcpp::wrap(data.mergedResultEndScanNum);
        names.push_back("injectionTime");
        header[i++] = 0;

        header.attr("names") = names;

        return  header;
    }
    Rprintf("Warning: Ramp not yet initialized.\n ");
    return Rcpp::List::create( );
}

Rcpp::DataFrame RcppRamp::getAllScanHeaderInfo ( )
{
    if (ramp != NULL)
    {
        if (!isInCacheAllScanHeaderInfo)
        {
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
            Rcpp::IntegerVector mergedScan(N);  /* only if MS level > 1 */
            Rcpp::IntegerVector mergedResultScanNum(N); /* scan number of the resultant merged scan */
            Rcpp::IntegerVector mergedResultStartScanNum(N); /* smallest scan number of the scanOrigin for merged scan */
            Rcpp::IntegerVector mergedResultEndScanNum(N); /* largest scan number of the scanOrigin for merged scan */

            for (int whichScan=1; whichScan <= N; whichScan++)
            {
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

            Rcpp::List header(22);
            std::vector<std::string> names;
            int i = 0;

            names.push_back("seqNum");
            header[i++] =Rcpp::wrap(seqNum);
            names.push_back("acquisitionNum");
            header[i++] =Rcpp::wrap( acquisitionNum);
            names.push_back("msLevel");
            header[i++] =Rcpp::wrap(msLevel);
            names.push_back("polarity");
            header[i++] =Rcpp::wrap(polarity);
            names.push_back("peaksCount");
            header[i++] =Rcpp::wrap(peaksCount);
            names.push_back("totIonCurrent");
            header[i++] =Rcpp::wrap(totIonCurrent);
            names.push_back("retentionTime");
            header[i++] =Rcpp::wrap(retentionTime);
            names.push_back("basePeakMZ");
            header[i++] =Rcpp::wrap(basePeakMZ);
            names.push_back("basePeakIntensity");
            header[i++] =Rcpp::wrap(basePeakIntensity);
            names.push_back("collisionEnergy");
            header[i++] =Rcpp::wrap(collisionEnergy);
            names.push_back("ionisationEnergy");
            header[i++] =Rcpp::wrap(ionisationEnergy);
            names.push_back("lowMZ");
            header[i++] =Rcpp::wrap(lowMZ);
            names.push_back("highMZ");
            header[i++] =Rcpp::wrap(highMZ);
            names.push_back("precursorScanNum");
            header[i++] =Rcpp::wrap(precursorScanNum);
            names.push_back("precursorMZ");
            header[i++] =Rcpp::wrap(precursorMZ);
            names.push_back("precursorCharge");
            header[i++] =Rcpp::wrap(precursorCharge);
            names.push_back("precursorIntensity");
            header[i++] =Rcpp::wrap(precursorIntensity);
            names.push_back("mergedScan");
            header[i++] =Rcpp::wrap(mergedScan);
            names.push_back("mergedResultScanNum");
            header[i++] =Rcpp::wrap(mergedResultScanNum);
            names.push_back("mergedResultStartScanNum");
            header[i++] =Rcpp::wrap(mergedResultStartScanNum);
            names.push_back("mergedResultEndScanNum");
            header[i++] =Rcpp::wrap(mergedResultEndScanNum);
            names.push_back("injectionTime");
            header[i++] = 0;
            
			header.attr("names") = names;
			
            allScanHeaderInfo = header;
            isInCacheAllScanHeaderInfo = TRUE;
        }
        else
        {
            // Rprintf("Read from cache.\n ");
        }
        return(allScanHeaderInfo);
    }
    Rprintf("Warning: Ramp not yet initialized.\n ");
    return Rcpp::DataFrame::create( );
}

Rcpp::List RcppRamp::getPeakList ( int whichScan )
{
    if (ramp != NULL)
    {
        if ((whichScan <= 0) || (whichScan > ramp->getLastScan()))
        {
            Rprintf("Index whichScan out of bounds [1 ... %d].\n", ramp->getLastScan());
            return Rcpp::List::create( );
        }
        rampPeakList *pl = ramp->getPeakList( whichScan );
        int peaksCount = 0;
        if (pl != NULL)
        {
            peaksCount = pl->getPeakCount();
        }
        Rcpp::NumericMatrix peaks(peaksCount, 2);
        if (pl != NULL)
        {
            rampPeakInfoStruct *peak;
            peak = pl->getPeak(0);
            peaks(0,0) = peak->mz;
            peaks(0,1) = peak->intensity;
            for (int i=1; i < peaksCount; i++)
            {
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

Rcpp::NumericMatrix RcppRamp::get3DMap ( std::vector<int> scanNumbers, double whichMzLow, double whichMzHigh, double resMz )
{
    if (ramp != NULL)
    {
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
        // map3d = 0.0;
        int j=0;
        Rprintf("%d\n",1);
        for (int i = 0; i < scanNumbers.size(); i++)
        {
            rampPeakList *pl = ramp->getPeakList( scanNumbers[i] );
            int peaksCount = pl->getPeakCount();
            rampPeakInfoStruct *peak;
            peak = pl->getPeak(0);
            j = round(peak->mz * f) - low;
            if ((j >= 0) & (j < dmz))
            {
                if (peak->intensity > map3d(i,j))
                {
                    map3d(i,j) = peak->intensity;
                }
            }
            for (int k=1; k < peaksCount; k++)
            {
                peak++;
                j = round(peak->mz * f) - low;
                if ((j >= 0) & (j < dmz))
                {
                    if (peak->intensity > map3d(i,j))
                    {
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

int RcppRamp::getLastScan() const
{
    if (ramp != NULL)
    {
        return ramp->getLastScan();
    }
    Rprintf("Warning: Ramp not yet initialized.\n ");
    return -1;
}

bool RcppRamp::OK (  )
{
    if (ramp != NULL)
    {
        return ramp->OK();
    }
    // Rprintf("Warning: Ramp not yet initialized.\n ");
    return false;
}
