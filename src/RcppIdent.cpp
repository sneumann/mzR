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
    
    vector<EnzymePtr> enz = sip[0]->enzymes.enzymes;
    Rcpp::List enzymes;
    Rcpp::StringVector name(enz.size());
    Rcpp::StringVector nTermGain(enz.size());
    Rcpp::StringVector cTermGain(enz.size());
    Rcpp::StringVector minDistance(enz.size());
    Rcpp::StringVector missedCleavages(enz.size());
    
    for (size_t i = 0; i < enz.size(); i++)
    {
		name[i] = cvTermInfo(cleavageAgent(*enz[i].get())).name;
		nTermGain[i] = enz[i]->nTermGain;
		cTermGain[i] = enz[i]->cTermGain;
		minDistance[i] = enz[i]->minDistance;
        missedCleavages[i] = enz[i]->missedCleavages;
    }
    
    enzymes = Rcpp::List::create(
					Rcpp::_["name"]  = name,
					Rcpp::_["nTermGain"]  = nTermGain,
					Rcpp::_["cTermGain"]  = cTermGain,
					Rcpp::_["minDistance"]  = minDistance,
					Rcpp::_["missedCleavages"]  = missedCleavages
    );
    
    vector<SpectraDataPtr> sd = mzid->dataCollection.inputs.spectraData;
    Rcpp::StringVector format(sd.size());
    for (size_t i = 0; i < sd.size(); i++)
    {
        format[i] = sd[i]->fileFormat.name();
    }
    
	return Rcpp::List::create(
                   Rcpp::_["FileProvider"]	= provider,
                   Rcpp::_["CreationDate"]	= date,
                   Rcpp::_["software"]	= software,
                   Rcpp::_["database"]	= database,
                   Rcpp::_["enzymes"]	= enzymes,
                   Rcpp::_["SpectraDataFormat"]	= format
               );
    
}

Rcpp::DataFrame RcppIdent::getPepInfo(  )
{
	vector<PeptidePtr> pep = mzid->sequenceCollection.peptides;
    Rcpp::StringVector seq(pep.size());
    Rcpp::NumericVector modification(pep.size());
	
		
    for (size_t i = 0; i < pep.size(); i++) {
		seq[i] = pep[i]->peptideSequence;
		modification[i] = pep[i]->modification.size();

    }
    
	return Rcpp::DataFrame::create(
                   Rcpp::_["sequence"]	= seq,
                   Rcpp::_["modNum"]	= modification
               );
    
}

Rcpp::List RcppIdent::getModInfo(  )
{
	vector<PeptidePtr> pep = mzid->sequenceCollection.peptides;
	std::vector<std::string> seq;
	std::vector<std::string> modification;
    		
    for (size_t i = 0; i < pep.size(); i++) {
		
		if(pep[i]->modification.size() > 0){
			for(size_t j = 0 ; j < pep[i]->modification.size(); j++){
				seq.push_back(pep[i]->peptideSequence);
				modification.push_back(cvTermInfo(pep[i]->modification[j]->cvParams[0].cvid).name);
			}
		}
    }
    
	return Rcpp::List::create(
                   Rcpp::_["sequence"]	= seq,
                   Rcpp::_["name"]	= modification
               );
    
}

