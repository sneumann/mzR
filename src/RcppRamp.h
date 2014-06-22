#ifndef _mzR_RCPP_RAMP_H
#define _mzR_RCPP_RAMP_H

#include "Rcpp.h"
// Taken from http://tolstoy.newcastle.edu.au/R/e2/devel/06/11/1242.html
// and http://stackoverflow.com/questions/11588765/using-rcpp-with-windows-specific-includes
// Undefine the Realloc macro, which is defined by both R and by Windows stuff
// Also need to undefine the Free macro
#if defined(__MINGW32__)
#undef Realloc
#undef Free
#endif

#include "cramp.h"

#if defined(__MINGW32__)
#include <windows.h>
#endif

class RcppRamp {

private:
  cRamp *ramp;
  Rcpp::List runInfo;
  bool isInCacheRunInfo;
  Rcpp::List instrumentInfo;
  bool isInCacheInstrumentInfo;
  Rcpp::DataFrame allScanHeaderInfo;
  bool isInCacheAllScanHeaderInfo;
  Rcpp::StringVector filename;

public:  
  RcppRamp();          // Constructor
  virtual ~RcppRamp(); // Destructor
  
  /** 
   * Opens a mass spec file (mzXML, mzData, etc.) and creates a cRamp object.
   * @param fileName: Name of the msxml file
   * @param declaredScansOnly: suppress RAMP's behavior of creating sparse 
   *        tables to accomodate unlisted scans
   */
  void open( const char* fileName, bool declaredScansOnly=false );

  /** 
   * Closes the mzXML file. Releases the memory of the cRamp object. 
   * This function allows memory management from R site.
   */
  void close();
  
  /**
   * Returns the filename.
   * @return the filename
   */
  Rcpp::StringVector getFilename (  );
  
  /**
   * Reads the run information from the mzXML header.
   * @return rampRunInfo* is parsed to a Rcpp::List
   */
  Rcpp::List getRunInfo (  );
  
  /**
   * Reads the instrument information from the mzXML header.
   * @return rampInstrumentInfo* is parsed to a Rcpp::List
   */
  Rcpp::List getInstrumentInfo();
  
  /**
   * Reads the scan header info. It reads for one mass spectrum the 
   * header information, like ms level, retention time, ion current, ... .
   * @return rampScanInfo* is parsed to a Rcpp::List
   */
  Rcpp::List getScanHeaderInfo ( int whichScan  );
  
  /**
   * Reads the scan header info. It reads for all mass spectra the header 
   * information, like ms level, retention time, ion current, ... .
   * @return rampScanInfo* is parsed to a Rcpp::List
   */
  Rcpp::DataFrame getAllScanHeaderInfo ( );
  
  /**
   * This function performs a non-sequential parsing operation on an indexed
   * mzXML file to obtain the peak list (i.e. the mass spectrum) for a numbered scan.
   * @param whichScan: Number of the scan we want to read from
   * @return rampPeakList* is parsed to a Rcpp::List. The first list element 
   *         contains the number of peaks, the second contains a n x 2 matrix 
   *         with the peak list. 
   */
  Rcpp::List getPeakList ( int whichScan );
  
  /**
   * This function reads all scans and returns them as a matrix. The number of 
   * rows is equal to the given number of scan numbers. The columns represent 
   * equidistant m/z values. mzXML file to obtain the peak list (i.e. the mass spectrum) 
   * for a numbered scan.
   * @param whichScan: The scan numbers we want to read from.
   * @param whichMzLow: The lowest m/z value to be returned. 
   * @param whichMzHigh: The highest m/z value to be returned.
   * @param resMz: The resolution in m/z direction.
   * @return The matrix is given back as a Rcpp::NumericMatrix.
   */
  Rcpp::NumericMatrix
  get3DMap ( std::vector<int> scanNumbers,
	     double whichMzLow,
	     double whichMzHigh,
	     double resMz );
  
  /**
   * Getting for the last scan number (not necessarily the number of scans 
   * because of missing scans).
   */
  int getLastScan() const;
  
  /**
   * checks the status of the object.
   */
  bool OK();
};

#endif
