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
	vector<PeptideEvidencePtr> peptideEvidence = mzid->sequenceCollection.peptideEvidence;
    Rcpp::StringVector seq(pep.size());
    Rcpp::NumericVector modification(pep.size());
    Rcpp::LogicalVector isDecoy(pep.size());
    Rcpp::StringVector post(pep.size());
    Rcpp::StringVector pre(pep.size());
	Rcpp::NumericVector start(pep.size());
	Rcpp::NumericVector end(pep.size());
	Rcpp::StringVector DBSequenceID(pep.size());
	
    for (size_t i = 0; i < pep.size(); i++) {
		seq[i] = pep[i]->peptideSequence;
		modification[i] = pep[i]->modification.size();
		isDecoy[i] = peptideEvidence[i]->isDecoy;
		post[i] =  string(1, peptideEvidence[i]->post);
		pre[i] = string(1, peptideEvidence[i]->pre);
		start[i] = peptideEvidence[i]->start;
		end[i] = peptideEvidence[i]->end;
		DBSequenceID[i] = peptideEvidence[i]->dbSequencePtr->id;
		
    }
    
	return Rcpp::DataFrame::create(
                   Rcpp::_["sequence"]	= seq,
                   Rcpp::_["modNum"]	= modification,
                   Rcpp::_["isDecoy"]	= isDecoy,
                   Rcpp::_["post"]		= post,
                   Rcpp::_["pre"]		= pre,
                   Rcpp::_["start"]		= start,
                   Rcpp::_["end"]		= end,
                   Rcpp::_["DatabaseID"]= DBSequenceID
               );
    
}

Rcpp::DataFrame RcppIdent::getModInfo(  )
{
	vector<PeptidePtr> pep = mzid->sequenceCollection.peptides;

	std::vector<std::string> seq;
	std::vector<std::string> name;
	std::vector<double> mass;
	std::vector<int> loc;
	
    		
    for (size_t i = 0; i < pep.size(); i++) {
		
		if(pep[i]->modification.size() > 0){
			for(size_t j = 0 ; j < pep[i]->modification.size(); j++){
				seq.push_back(pep[i]->peptideSequence);
				name.push_back(cvTermInfo(pep[i]->modification[j]->cvParams[0].cvid).name);
				mass.push_back(pep[i]->modification[j]->monoisotopicMassDelta);
				loc.push_back(pep[i]->modification[j]->location);
			}
		}
    }
    
	return Rcpp::DataFrame::create(
                   Rcpp::_["sequence"]	= seq,
                   Rcpp::_["name"]	= name,
                   Rcpp::_["mass"]	= mass,
                   Rcpp::_["location"]	= loc
               );
    
}

Rcpp::DataFrame RcppIdent::getSubInfo(  )
{
	vector<PeptidePtr> pep = mzid->sequenceCollection.peptides;
	std::vector<std::string> seq;
	std::vector<char> originalResidue;
	std::vector<char> replacementResidue;
	std::vector<int> loc;
	
    		
    for (size_t i = 0; i < pep.size(); i++) {
		
		if(pep[i]->substitutionModification.size() > 0){
			for(size_t j = 0 ; j < pep[i]->substitutionModification.size(); j++){
				seq.push_back(pep[i]->peptideSequence);
				originalResidue.push_back(pep[i]->substitutionModification[j]->originalResidue);
				replacementResidue.push_back(pep[i]->substitutionModification[j]->replacementResidue);
				loc.push_back(pep[i]->substitutionModification[j]->location);
			}
		}
    }
    
	return Rcpp::DataFrame::create(
                   Rcpp::_["sequence"]	= seq,
                   Rcpp::_["originalResidue"]	= originalResidue,
                   Rcpp::_["replacementResidue"]	= replacementResidue,
                   Rcpp::_["location"]	= loc
               );
    
}

