#ifndef _mzR_RCPP_IDENT_H
#define _mzR_RCPP_IDENT_H

#include "Rcpp.h"

#if defined(__MINGW32__)
#undef Realloc
#undef Free
#endif

#include "pwiz/data/identdata/IdentDataFile.hpp"
#include "pwiz/data/identdata/DefaultReaderList.hpp"
#include "pwiz/data/identdata/IO.hpp"
#include "pwiz/utility/misc/Std.hpp"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#if defined(__MINGW32__)
#include <windows.h>
#endif

using namespace pwiz::identdata;

class RcppIdent
{
 private:

  IdentDataFile *mzid;
  string date, filename, provider;

 public:

  RcppIdent();

  void open(const string& fileNames);

  Rcpp::List getIDInfo();

  Rcpp::DataFrame getPsmInfo();

  Rcpp::DataFrame getModInfo();

  Rcpp::DataFrame getSubInfo();

  Rcpp::DataFrame getScore();

  Rcpp::List getPara();

  Rcpp::DataFrame getDB();

  inline bool isNumber(const std::string& s)
  {
    std::istringstream iss( s );
    double dTestSink;
    iss >> dTestSink;
    return ( iss.rdbuf()->in_avail() == 0 );
  }

  inline bool isBool(std::string s)
  {

    boost::algorithm::to_lower(s);
    if (s=="true"||s=="false")
      return true;
    else
      return false;
  }

  inline bool toBool(std::string s)
  {
    boost::algorithm::to_lower(s);
    if (s=="true")
      return true;
    else
      return false;
  }

  inline std::string underscore(std::string text)
  {
    for(std::string::iterator it = text.begin(); it != text.end(); ++it)
    {
      if(*it == ' '||*it == ':'||*it == '-'||*it == '!')
      {
        *it = '_';
      }
    }
    return text;
  }
};

#endif
