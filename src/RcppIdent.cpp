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

    Rcpp::StringVector database(sdb.size());
    for (size_t i = 0; i < sdb.size(); i++)
    {
        database = sdb[i]->location + " (" + lexical_cast<string>(sdb[i]->numDatabaseSequences) + " sequences)";
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
               Rcpp::_["database"]	= database,
               Rcpp::_["ModificationSearched"]	= mod,
               Rcpp::_["FragmentTolerance"]	= fragmentTolerance,
               Rcpp::_["ParentTolerance"]	= parentTolerance,
               Rcpp::_["enzymes"]	= enzymes,
               Rcpp::_["SpectraSource"]	= spectra
           );

}

Rcpp::DataFrame RcppIdent::getPepInfo(  )
{
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

    }

    return Rcpp::DataFrame::create(
               Rcpp::_["spectrumID"]	= spectrumID,
               Rcpp::_["chargeState"]	= chargeState,
               Rcpp::_["experimentalM/Z"]	= experimentalMassToCharge,
               Rcpp::_["calculatedM/Z"]	= calculatedMassToCharge,
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
    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;
    vector<string> spectrumID;
    vector<string> seq;
    vector<string> name;
    vector<double> mass;
    vector<int> loc;

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
            }
        }

    }

    return Rcpp::DataFrame::create(
               Rcpp::_["spectrumID"]	= spectrumID,
               Rcpp::_["sequence"]	= seq,
               Rcpp::_["name"]	= name,
               Rcpp::_["mass"]	= mass,
               Rcpp::_["location"]	= loc
           );

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
    //names.push_back("spectrumID");
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

Rcpp::CharacterVector RcppIdent::getPara(  )
{
    vector<SpectrumIdentificationProtocolPtr> sip = mzid->analysisProtocolCollection.spectrumIdentificationProtocol;
    vector<SearchModificationPtr> sm = sip[0]->modificationParams;
    std::vector<std::string> para;
    para.push_back(cvTermInfo(sip[0]->searchType.cvid).name);

    for(int i = 0 ; i < sip[0]->additionalSearchParams.cvParams.size(); i++)
    {
        para.push_back(cvTermInfo(sip[0]->additionalSearchParams.cvParams[i].cvid).name);
    }

    for(int i = 0; i < sip[0]->additionalSearchParams.userParams.size(); i++)
    {
        if(sip[0]->additionalSearchParams.userParams[i].value.empty())
        {
            para.push_back(sip[0]->additionalSearchParams.userParams[i].name);
        }
        else
            para.push_back(sip[0]->additionalSearchParams.userParams[i].name + ": " + sip[0]->additionalSearchParams.userParams[i].value);
    }

    return wrap(para);
}

Rcpp::DataFrame RcppIdent::getDB(  )
{
    ListBuilder res;
    std::vector<std::string> access;
    std::vector<std::string> seq;
    std::vector<std::string> desc;
    for(int i = 0; i < mzid->sequenceCollection.dbSequences.size(); i++)
    {
        //Rcpp::Rcout << i << "\t" << mzid->sequenceCollection.dbSequences[i]->accession << std::endl;
        access.push_back(mzid->sequenceCollection.dbSequences[i]->accession);
        if(!mzid->sequenceCollection.dbSequences[i]->seq.empty())
            seq.push_back(mzid->sequenceCollection.dbSequences[i]->seq);
        if(mzid->sequenceCollection.dbSequences[i]->cvParams.size() > 0)
            desc.push_back(mzid->sequenceCollection.dbSequences[i]->cvParams[0].value);
    }
    res.add("accession", Rcpp::wrap(access));
    if(!mzid->sequenceCollection.dbSequences[0]->seq.empty())
        res.add("seq", Rcpp::wrap(seq));
    if(mzid->sequenceCollection.dbSequences[0]->cvParams.size() > 0 )
        res.add("description", Rcpp::wrap(desc));

    return res;
}
