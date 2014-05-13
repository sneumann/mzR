#ifndef _mzR_RCPP_PWIZ_H
#define _mzR_RCPP_PWIZ_H

#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/LegacyAdapter.hpp"
#include "pwiz/data/common/CVTranslator.hpp"
#include "pwiz/utility/misc/Std.hpp"

#include "Rcpp.h"

using namespace pwiz::msdata;

class RcppPwiz {

private:
  MSDataFile *msd;
  Rcpp::List runInfo;
  bool isInCacheRunInfo;
  Rcpp::List instrumentInfo;
  bool isInCacheInstrumentInfo;
  Rcpp::DataFrame allScanHeaderInfo;
  bool isInCacheAllScanHeaderInfo;
  Rcpp::StringVector filename;

public: 

  RcppPwiz();

  void open(const string& fileNames);

  Rcpp::StringVector getFilename (  );
  
  int getLastScan() const;
  
  Rcpp::List getInstrumentInfo();
  
};

#endif
