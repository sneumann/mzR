#ifndef _mzR_RCPP_PWIZ_H
#define _mzR_RCPP_PWIZ_H

#include "Rcpp.h"

#if defined(__MINGW32__)
#undef Realloc
#undef Free
#endif

#include "pwiz/data/msdata/RAMPAdapter.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/LegacyAdapter.hpp"
#include "pwiz/data/msdata/Serializer_mzML.hpp"
#include "pwiz/data/msdata/Serializer_mzXML.hpp"
#include "pwiz/data/msdata/Serializer_MGF.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/common/CVTranslator.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/copy.hpp>

#include <fstream>
#include <string>
#include <iostream>

#if defined(__MINGW32__)
#include <windows.h>
#endif

using namespace pwiz::cv;
using namespace pwiz::msdata;
using namespace pwiz::util;


class RcppPwiz
{

 private:
  MSDataFile *msd;
  Rcpp::List instrumentInfo;
  Rcpp::DataFrame chromatogramsInfo;
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

  Rcpp::List getRunInfo();

  Rcpp::List getScanHeaderInfo(int whichScan);

  Rcpp::DataFrame getChromatogramsInfo();

  Rcpp::DataFrame getAllScanHeaderInfo();

  Rcpp::List getPeakList(int whichScan);

  Rcpp::NumericMatrix get3DMap(std::vector<int> scanNumbers, double whichMzLow, double whichMzHigh, double resMz);

};

#endif
