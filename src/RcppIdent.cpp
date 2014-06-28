#include "RcppIdent.h"

RcppIdent::RcppIdent()
{
    mzid = NULL;
}

void RcppIdent::open(const string& fileName)
{

    filename = fileName;
    mzid = new IdentDataFile(fileName);

}

string RcppIdent::getCreationDate(  )
{

    return mzid->creationDate;
}

