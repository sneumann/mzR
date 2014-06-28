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

Rcpp::List RcppIdent::getIDInfo(  )
{
	provider = (mzid->provider.contactRolePtr.get()!=0?mzid->provider.contactRolePtr.get()->name():"");
    date = mzid->creationDate;
    vector<AnalysisSoftwarePtr> as = mzid->analysisSoftwareList;
    int N = as.size();
    Rcpp::StringVector software(N);
    for (size_t i = 0; i < as.size(); i++)
    {
		software[i] = as[i]->name + " " + as[i]->version + " " + (as[i]->contactRolePtr.get()!=0?as[i]->contactRolePtr->contactPtr->name:"") ;
    }
	return Rcpp::List::create(
                   Rcpp::_["FileProvider"]	= provider,
                   Rcpp::_["CreationDate"]	= date,
                   Rcpp::_["software"]	= software
               );
    
}

