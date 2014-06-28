#ifndef _mzR_RCPP_IDENT_H
#define _mzR_RCPP_IDENT_H

#include "pwiz/data/identdata/IdentDataFile.hpp"
#include "pwiz/data/identdata/DefaultReaderList.hpp"
#include "pwiz/data/identdata/IO.hpp"
#include "pwiz/utility/misc/Std.hpp"

#include <fstream>
#include <string>
#include <iostream>

#include "Rcpp.h"

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
};

#endif
