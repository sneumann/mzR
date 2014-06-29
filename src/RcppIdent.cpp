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
    vector<SearchDatabasePtr> sdb = mzid->dataCollection.inputs.searchDatabase;
    Rcpp::StringVector software(as.size());

    for (size_t i = 0; i < as.size(); i++)
    {
		software[i] = as[i]->name + " " + as[i]->version + " " + (as[i]->contactRolePtr.get()!=0?as[i]->contactRolePtr->contactPtr->name:"") ;
    }
    
    Rcpp::StringVector database(sdb.size());
    for (size_t i = 0; i < sdb.size(); i++)
    {
		database = sdb[i]->name + " (" + lexical_cast<string>(sdb[i]->numDatabaseSequences) + " sequences)";
    }
    
    vector<SpectrumIdentificationProtocolPtr> sip = mzid->analysisProtocolCollection.spectrumIdentificationProtocol;
    vector<SearchModificationPtr> sm = sip[0]->modificationParams;
    Rcpp::StringVector mod(sm.size());
    
    for (size_t i = 0; i < sm.size(); i++)
    {
		mod[i] = lexical_cast<string>(sm[i]->massDelta) + " (";
		for (size_t j = 0; j < sm[i]->residues.size(); j++)
        {
			mod[i] += sm[i]->residues[j];
        }        
        mod[i] += ")";
    }
    
	return Rcpp::List::create(
                   Rcpp::_["FileProvider"]	= provider,
                   Rcpp::_["CreationDate"]	= date,
                   Rcpp::_["software"]	= software,
                   Rcpp::_["database"]	= database,
                   Rcpp::_["modification"]	= mod
               );
    
}

