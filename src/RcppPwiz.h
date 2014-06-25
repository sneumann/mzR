#ifndef _mzR_RCPP_PWIZ_H
#define _mzR_RCPP_PWIZ_H

#include "pwiz/data/msdata/RAMPAdapter.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/LegacyAdapter.hpp"
#include "pwiz/data/msdata/Serializer_mz5.hpp"
#include "pwiz/data/msdata/Serializer_mzML.hpp"
#include "pwiz/data/common/CVTranslator.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>

#include "Rcpp.h"

using namespace pwiz::cv;
using namespace pwiz::msdata;
using namespace pwiz::util;

class RcppPwiz {

private:
  MSDataFile *msd;
  Rcpp::List runInfo;
  bool isInCacheRunInfo;
  Rcpp::List instrumentInfo;
  Rcpp::List chromatogramsInfo;
  bool isInCacheInstrumentInfo;
  Rcpp::DataFrame allScanHeaderInfo;
  bool isInCacheAllScanHeaderInfo;
  string filename;

public: 

  RcppPwiz();

  void open(const string& fileNames);

  string getFilename();
  
  int getLastScan() const;
  
  Rcpp::List getInstrumentInfo();
  
  Rcpp::List getScanHeaderInfo(int whichScan);
  
  Rcpp::List getChromatogramsInfo();
  
  Rcpp::DataFrame getAllScanHeaderInfo();
  
  Rcpp::List getPeakList(int whichScan);
  
  Rcpp::NumericMatrix get3DMap(std::vector<int> scanNumbers, double whichMzLow, double whichMzHigh, double resMz);
    
};

#endif
