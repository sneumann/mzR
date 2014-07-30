#include <boost/lexical_cast.hpp>

#include "RcppIdent.h"
#include "ListBuilder.h"

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

    vector<SpectrumIdentificationProtocolPtr> sip = mzid->analysisProtocolCollection.spectrumIdentificationProtocol;
    string fragmentTolerance = "";
    string parentTolerance = "";
    if(!sip[0]->fragmentTolerance.empty())
    {
        fragmentTolerance = sip[0]->fragmentTolerance.cvParams[0].value + " " + sip[0]->fragmentTolerance.cvParam(MS_search_tolerance_plus_value).unitsName();
    }

    if(!sip[0]->parentTolerance.empty())
    {
        parentTolerance = sip[0]->parentTolerance.cvParams[0].value + " " + sip[0]->parentTolerance.cvParam(MS_search_tolerance_plus_value).unitsName();
    }

    vector<SearchModificationPtr> sm = sip[0]->modificationParams;

    Rcpp::StringVector mod(sm.size());
    for(size_t i = 0; i < sm.size(); i++)
    {
        mod[i] = cvTermInfo(sm[i]->cvParams[0].cvid).name;
    }

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
    Rcpp::StringVector spectra(sd.size());
    for (size_t i = 0; i < sd.size(); i++)
    {
        spectra[i] = sd[i]->location;
    }

    return Rcpp::List::create(
               Rcpp::_["FileProvider"]	= provider,
               Rcpp::_["CreationDate"]	= date,
               Rcpp::_["software"]	= software,
               Rcpp::_["ModificationSearched"]	= mod,
               Rcpp::_["FragmentTolerance"]	= fragmentTolerance,
               Rcpp::_["ParentTolerance"]	= parentTolerance,
               Rcpp::_["enzymes"]	= enzymes,
               Rcpp::_["SpectraSource"]	= spectra
           );

}

Rcpp::DataFrame RcppIdent::getPepInfo(  )
{
    ListBuilder res;

    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;

    Rcpp::StringVector spectrumID(spectrumIdResult.size());
    Rcpp::NumericVector chargeState(spectrumIdResult.size());
    Rcpp::NumericVector rank(spectrumIdResult.size());
    Rcpp::NumericVector experimentalMassToCharge(spectrumIdResult.size());
    Rcpp::NumericVector calculatedMassToCharge(spectrumIdResult.size());
    Rcpp::StringVector seq(spectrumIdResult.size());
    Rcpp::NumericVector modification(spectrumIdResult.size());
    Rcpp::LogicalVector isDecoy(spectrumIdResult.size());
    Rcpp::StringVector post(spectrumIdResult.size());
    Rcpp::StringVector pre(spectrumIdResult.size());
    Rcpp::NumericVector start(spectrumIdResult.size());
    Rcpp::NumericVector end(spectrumIdResult.size());
    Rcpp::StringVector DBSequenceID(spectrumIdResult.size());
    Rcpp::StringVector DBseq(spectrumIdResult.size());
    Rcpp::StringVector DBdesc(spectrumIdResult.size());

    for (size_t i = 0; i < spectrumIdResult.size(); i++)
    {

        spectrumID[i] = spectrumIdResult[i]->spectrumID;
        chargeState[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->chargeState;
        experimentalMassToCharge[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->experimentalMassToCharge;
        calculatedMassToCharge[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->calculatedMassToCharge;
        seq[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->peptideSequence;
        modification[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->modification.size();
        isDecoy[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->isDecoy;
        pre[i] = string(1, spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->pre);
        post[i] =  string(1, spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->post);
        start[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->start;
        end[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->end;
        DBSequenceID[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->accession;
        DBseq[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->seq;
        if(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->cvParams.size() > 0)
            DBdesc[i] = spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->cvParams[0].value;
    }

    res.add("spectrumID", Rcpp::wrap(spectrumID));
    res.add("chargeState", Rcpp::wrap(chargeState));
    res.add("experimentalM/Z", Rcpp::wrap(experimentalMassToCharge));
    res.add("calculatedM/Z", Rcpp::wrap(calculatedMassToCharge));
    res.add("sequence", Rcpp::wrap(seq));
    res.add("modNum", Rcpp::wrap(modification));
    res.add("isDecoy", Rcpp::wrap(isDecoy));
    res.add("post", Rcpp::wrap(post));
    res.add("pre", Rcpp::wrap(pre));
    res.add("start", Rcpp::wrap(start));
    res.add("end", Rcpp::wrap(end));
    res.add("DatabaseID", Rcpp::wrap(DBSequenceID));
    if(!spectrumIdResult[0]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->seq.empty())
        res.add("DatabaseSeq", Rcpp::wrap(DBseq));
    if(spectrumIdResult[0]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->cvParams.size() > 0)
        res.add("DatabaseDescription", Rcpp::wrap(DBdesc));

    return res;
}

Rcpp::DataFrame RcppIdent::getModInfo(  )
{
    ListBuilder res;

    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;
    vector<string> spectrumID;
    vector<string> seq;
    vector<string> name;
    vector<double> mass;
    vector<int> loc;
    vector<string> DBSequenceID;
    vector<string> DBseq;
    vector<string> DBdesc;

    for (size_t i = 0; i < spectrumIdResult.size(); i++)
    {
        if(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->modification.size()>0)
        {
            for(size_t j = 0 ; j < spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->modification.size(); j++)
            {
                spectrumID.push_back(spectrumIdResult[i]->spectrumID);
                seq.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->peptideSequence);
                name.push_back(cvTermInfo(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->modification[j]->cvParams[0].cvid).name);
                mass.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->modification[j]->monoisotopicMassDelta);
                loc.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->modification[j]->location);
                DBSequenceID.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->accession);
                DBseq.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->seq);
                if(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->cvParams.size() > 0)
                    DBdesc.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->cvParams[0].value);
            }
        }

    }

    res.add("spectrumID", Rcpp::wrap(spectrumID));
    res.add("sequence", Rcpp::wrap(seq));
    res.add("name", Rcpp::wrap(name));
    res.add("mass", Rcpp::wrap(mass));
    res.add("location", Rcpp::wrap(loc));
    res.add("DBSequenceID", Rcpp::wrap(DBSequenceID));
    if(!spectrumIdResult[0]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->seq.empty())
        res.add("DatabaseSeq", Rcpp::wrap(DBseq));
    if(spectrumIdResult[0]->spectrumIdentificationItem[0]->peptideEvidencePtr[0]->dbSequencePtr->cvParams.size() > 0)
        res.add("DatabaseDescription", Rcpp::wrap(DBdesc));
    return res;
}

Rcpp::DataFrame RcppIdent::getSubInfo(  )
{
    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;
    vector<string> spectrumID;
    std::vector<std::string> seq;
    std::vector<char> originalResidue;
    std::vector<char> replacementResidue;
    std::vector<int> loc;


    for (size_t i = 0; i < spectrumIdResult.size(); i++)
    {

        if(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->substitutionModification.size() > 0)
        {
            for(size_t j = 0 ; j < spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->substitutionModification.size(); j++)
            {
                spectrumID.push_back(spectrumIdResult[i]->spectrumID);
                seq.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->peptideSequence);
                originalResidue.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->substitutionModification[j]->originalResidue);
                replacementResidue.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->substitutionModification[j]->replacementResidue);
                loc.push_back(spectrumIdResult[i]->spectrumIdentificationItem[0]->peptidePtr->substitutionModification[j]->location);
            }
        }
    }

    return Rcpp::DataFrame::create(
               Rcpp::_["spectrumID"]	= spectrumID,
               Rcpp::_["sequence"]	= seq,
               Rcpp::_["originalResidue"]	= originalResidue,
               Rcpp::_["replacementResidue"]	= replacementResidue,
               Rcpp::_["location"]	= loc
           );

}

Rcpp::DataFrame RcppIdent::getScore(  )
{
    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;
    vector<string> spectrumID;
    vector<string> names;
    int count = 0;

    for(size_t i = 0; i < spectrumIdResult[0]->spectrumIdentificationItem[0]->cvParams.size(); i++)
    {
        if(!spectrumIdResult[0]->spectrumIdentificationItem[0]->cvParams[i].value.empty())
        {
            count++;
            names.push_back(cvTermInfo(spectrumIdResult[0]->spectrumIdentificationItem[0]->cvParams[i].cvid).name);
        }
    }
    if(count == 0)
    {
        Rcpp::Rcout << "No scoring information available" << std::endl;
        return Rcpp::DataFrame::create();
    }
    else
    {
        vector<vector<double> > score(count);

        for (size_t i = 0; i < spectrumIdResult.size(); i++)
        {
            spectrumID.push_back(spectrumIdResult[i]->spectrumID);
            count = 0;
            for(size_t j = 0; j < spectrumIdResult[i]->spectrumIdentificationItem[0]->cvParams.size(); j++)
            {
                if(!spectrumIdResult[i]->spectrumIdentificationItem[0]->cvParams[j].value.empty())
                {
                    score[count].push_back(lexical_cast<double>(spectrumIdResult[i]->spectrumIdentificationItem[0]->cvParams[j].value));
                    count++;
                }

            }
        }

        ListBuilder res;
        res.add("spectrumID", Rcpp::wrap(spectrumID));
        for(size_t i = 0; i < score.size(); i++)
        {
            res.add(names[i], Rcpp::wrap(score[i]));
        }

        return res;
    }


}

Rcpp::List RcppIdent::getPara(  )
{
	ListBuilder para;
    vector<SpectrumIdentificationProtocolPtr> sip = mzid->analysisProtocolCollection.spectrumIdentificationProtocol;

    para.add("searchType", Rcpp::wrap(cvTermInfo(sip[0]->searchType.cvid).name));

    for(int i = 0 ; i < sip[0]->additionalSearchParams.cvParams.size(); i++)
    {
        para.add(cvTermInfo(sip[0]->additionalSearchParams.cvParams[i].cvid).name, Rcpp::wrap((bool) 1));
    }

    for(int i = 0; i < sip[0]->additionalSearchParams.userParams.size(); i++)
    {
        if(sip[0]->additionalSearchParams.userParams[i].value.empty())
        {
            para.add(sip[0]->additionalSearchParams.userParams[i].name, Rcpp::wrap((bool) 1));
        }
        else
            para.add(sip[0]->additionalSearchParams.userParams[i].name, Rcpp::wrap(sip[0]->additionalSearchParams.userParams[i].value));
    }

    return para;
}

Rcpp::DataFrame RcppIdent::getDB(  )
{
    vector<SearchDatabasePtr> sdb = mzid->dataCollection.inputs.searchDatabase;

    std::vector<std::string> dbLocation;
    std::vector<std::string> dbID;
    std::vector<std::string> dbName;
    std::vector<std::string> dbVersion;
    std::vector<long> numDatabaseSequences;
    std::vector<long> numResidues;
    for (size_t i = 0; i < sdb.size(); i++)
    {
        dbLocation.push_back(sdb[i]->location);
        dbID.push_back(sdb[i]->id);
        dbName.push_back(sdb[i]->name);
        dbVersion.push_back(sdb[i]->version);
        numDatabaseSequences.push_back(sdb[i]->numDatabaseSequences);
        numResidues.push_back(sdb[i]->numResidues);
    }
    Rcpp::DataFrame database = Rcpp::List::create(
                                   Rcpp::_["location"]  = dbLocation,
                                   Rcpp::_["id"]  = dbID,
                                   Rcpp::_["name"]  = dbName,
                                   Rcpp::_["numDatabaseSequences"]  = numDatabaseSequences,
                                   Rcpp::_["numResidues"]  = numResidues,
                                   Rcpp::_["version"]  = dbVersion
                               );

    return database;
}
