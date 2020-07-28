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

Rcpp::DataFrame RcppIdent::getPsmInfo(  )
{
    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;

    std::vector<std::string> spectrumID;
    std::vector<int> chargeState;
    std::vector<int> rank;
    std::vector<double> experimentalMassToCharge;
    std::vector<double> calculatedMassToCharge;
    std::vector<std::string> seq;
    std::vector<std::string> peptide_ref;
    std::vector<int> modification;
    std::vector<bool> isDecoy;
    std::vector<bool> passThreshold;
    std::vector<std::string> post;
    std::vector<std::string> pre;
    std::vector<int> start;
    std::vector<int> end;
    std::vector<std::string> DBSequenceID;
    std::vector<std::string> DBseq;
    std::vector<int> DBSequenceLen;
    std::vector<std::string> DBdesc;

    for (size_t i = 0; i < spectrumIdResult.size(); i++)
    {
	for(size_t j = 0; j < spectrumIdResult[i]->spectrumIdentificationItem.size(); j++)
	{
	    for(size_t k = 0; k < spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr.size(); k++)
	    {

		spectrumID.push_back(spectrumIdResult[i]->spectrumID);
		chargeState.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->chargeState);
		rank.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->rank);
		passThreshold.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->passThreshold);
		experimentalMassToCharge.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->experimentalMassToCharge);
		calculatedMassToCharge.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->calculatedMassToCharge);
		seq.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptidePtr->peptideSequence);
		peptide_ref.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptidePtr->id);
		modification.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptidePtr->modification.size());
		isDecoy.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->isDecoy);
		pre.push_back(string(1, spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->pre));
		post.push_back(string(1, spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->post));
		start.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->start);
		end.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->end);
		if(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->dbSequencePtr.get()!=0)
		{
		    DBSequenceID.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->dbSequencePtr->accession);
		    DBSequenceLen.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->dbSequencePtr->length);
		    DBseq.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->dbSequencePtr->seq);
		    if(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->dbSequencePtr->cvParams.size() > 0)
		    {
			DBdesc.push_back(spectrumIdResult[i]->spectrumIdentificationItem[j]->peptideEvidencePtr[k]->dbSequencePtr->cvParams[0].value);
		    }
		    else
		    {
			DBdesc.push_back("");
		    }
		}
		else
		{
		    DBSequenceID.push_back("");
		    DBseq.push_back("");
		    DBdesc.push_back("");
		}

	    }
	}
    }



    return Rcpp::DataFrame::create(
	       Rcpp::_["spectrumID"]	= spectrumID,
	       Rcpp::_["chargeState"]	= chargeState,
	       Rcpp::_["rank"]	= rank,
	       Rcpp::_["passThreshold"]	= passThreshold,
	       Rcpp::_["experimentalMassToCharge"]	= experimentalMassToCharge,
	       Rcpp::_["calculatedMassToCharge"]	= calculatedMassToCharge,
	       Rcpp::_["sequence"]	= seq,
	       Rcpp::_["peptide_ref"]	= peptide_ref,
	       Rcpp::_["modNum"]	= modification,
	       Rcpp::_["isDecoy"]	= isDecoy,
	       Rcpp::_["post"]	= post,
	       Rcpp::_["pre"]	= pre,
	       Rcpp::_["start"]	= start,
	       Rcpp::_["end"]	= end,
	       Rcpp::_["DatabaseAccess"]	= DBSequenceID,
	       Rcpp::_["DBseqLength"]	= DBSequenceLen,
	       Rcpp::_["DatabaseSeq"]	= DBseq,
	       Rcpp::_["DatabaseDescription"]	= DBdesc
	   );
}

Rcpp::DataFrame RcppIdent::getModInfo(  )
{

    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;
    vector<string> spectrumID;
    vector<string> seq;
    vector<string> peptide_ref;
    vector<string> name;
    vector<double> mass;
    vector<int> loc;

    for (size_t i = 0; i < spectrumIdResult.size(); i++)
    {
	for(size_t k = 0; k < spectrumIdResult[i]->spectrumIdentificationItem.size() ; k++)
	{
	    if(spectrumIdResult[i]->spectrumIdentificationItem[k]->peptidePtr->modification.size()>0)
	    {
		for(size_t j = 0 ; j < spectrumIdResult[i]->spectrumIdentificationItem[k]->peptidePtr->modification.size(); j++)
		{
		    spectrumID.push_back(spectrumIdResult[i]->spectrumID);
		    seq.push_back(spectrumIdResult[i]->spectrumIdentificationItem[k]->peptidePtr->peptideSequence);
		    peptide_ref.push_back(spectrumIdResult[i]->spectrumIdentificationItem[k]->peptidePtr->id);
		    name.push_back(cvTermInfo(spectrumIdResult[i]->spectrumIdentificationItem[k]->peptidePtr->modification[j]->cvParams[0].cvid).name);
		    mass.push_back(spectrumIdResult[i]->spectrumIdentificationItem[k]->peptidePtr->modification[j]->monoisotopicMassDelta);
		    loc.push_back(spectrumIdResult[i]->spectrumIdentificationItem[k]->peptidePtr->modification[j]->location);
		}
	    }
	}
    }

    return Rcpp::DataFrame::create(
	       Rcpp::_["spectrumID"]	= spectrumID,
	       Rcpp::_["sequence"]	= seq,
	       Rcpp::_["peptide_ref"]	= peptide_ref,
	       Rcpp::_["name"]	= name,
	       Rcpp::_["mass"]	= mass,
	       Rcpp::_["location"]	= loc);

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

Rcpp::DataFrame RcppIdent::getScore(  ) {
  vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;
  vector<string> spectrumID;
  vector<string> names;
  int count = 0;
  int nCvParams = 0;

  for (size_t i = 0; i < spectrumIdResult[0]->spectrumIdentificationItem[0]->cvParams.size(); i++) {
    if (!spectrumIdResult[0]->spectrumIdentificationItem[0]->cvParams[i].value.empty()) {
      count++;
      nCvParams++;
      names.push_back(cvTermInfo(spectrumIdResult[0]->spectrumIdentificationItem[0]->cvParams[i].cvid).name);
    }
  }
  if(count == 0) {
    Rcpp::Rcout << "No scoring information available" << std::endl;
    return Rcpp::DataFrame::create();
  } else {
    vector<vector<double> > score(count);
    for (size_t i = 0; i < spectrumIdResult.size(); i++) {
      for (size_t k = 0; k < spectrumIdResult[i]->spectrumIdentificationItem.size(); k++) {
	for (size_t n = 0; n < spectrumIdResult[i]->spectrumIdentificationItem[k]->peptideEvidencePtr.size(); n++) {
	  spectrumID.push_back(spectrumIdResult[i]->spectrumID);
	  count = 0;

	  // The original loop iterated to j <
	  // spectrumIdResult[i]->spectrumIdentificationItem[k]->cvParams.size()
	  // which failed when some SpectrumIdentificationItem
	  // suddently have additional cvParams, such as in Mascot
	  // results - see https://github.com/sneumann/mzR/issues/136
	  for (size_t j = 0; j < nCvParams; j++) {
	    if (!spectrumIdResult[i]->spectrumIdentificationItem[k]->cvParams[j].value.empty()) {
	      score[count].push_back(lexical_cast<double>(spectrumIdResult[i]->spectrumIdentificationItem[k]->cvParams[j].value));
	      count++;
	    }
	  }
	}
      }
    }

    Rcpp::List res(score.size() + 1);
    names.insert(names.begin(), "spectrumID");
    res[0] = Rcpp::wrap(spectrumID);
    for(size_t i = 0; i < score.size(); i++) {
      res[i + 1] = Rcpp::wrap(score[i]);
    }

    res.attr("names") = names;
    Rcpp::DataFrame out(res);

    return out;
  }
}

Rcpp::DataFrame RcppIdent::getSpecParams(  )
{
    vector<SpectrumIdentificationResultPtr> spectrumIdResult = mzid->analysisCollection.spectrumIdentification[0]->spectrumIdentificationListPtr->spectrumIdentificationResult;
    vector<string> spectrumID;
    vector<string> names;
    int count = 0;

    for(size_t i = 0; i < spectrumIdResult[0]->cvParams.size(); i++)
    {
        if(!spectrumIdResult[0]->cvParams[i].value.empty())
        {
            count++;
            names.push_back(cvTermInfo(spectrumIdResult[0]->cvParams[i].cvid).name);
        }
    }
    if(count == 0)
    {
        Rcpp::Rcout << "No spectrum cvParams available" << std::endl;
        return Rcpp::DataFrame::create();
    }
    else
    {
        vector<vector<string> > score(count);

        for (size_t i = 0; i < spectrumIdResult.size(); i++)
        {
            spectrumID.push_back(spectrumIdResult[i]->spectrumID);
            count = 0;
            for(size_t j = 0; j < spectrumIdResult[i]->cvParams.size(); j++)
            {
                if(!spectrumIdResult[i]->cvParams[j].value.empty())
                {
                    score[count].push_back(lexical_cast<string>(spectrumIdResult[i]->cvParams[j].value));
                    count++;
                }
            }
        }

        Rcpp::List res(score.size() + 1);

        names.insert(names.begin(), "spectrumID");

        res[0] = Rcpp::wrap(spectrumID);

        for(size_t i = 0; i < score.size(); i++)
        {
            res[i + 1] = Rcpp::wrap(score[i]);
        }

        res.attr("names") = names;
        Rcpp::DataFrame out(res);

        return out;
    }
}

Rcpp::List RcppIdent::getPara(  )
{

    std::vector<SpectrumIdentificationProtocolPtr> sip = mzid->analysisProtocolCollection.spectrumIdentificationProtocol;
	std::vector<std::string> names, values;

	names.push_back("searchType");
	values.push_back(underscore(cvTermInfo(sip[0]->searchType.cvid).name));

    for(int i = 0 ; i < sip[0]->additionalSearchParams.cvParams.size(); i++)
    {
		names.push_back(underscore(cvTermInfo(sip[0]->additionalSearchParams.cvParams[i].cvid).name));
		values.push_back("true");
	}

    for(int i = 0; i < sip[0]->additionalSearchParams.userParams.size(); i++)
    {
		names.push_back(underscore(sip[0]->additionalSearchParams.userParams[i].name));
	if(sip[0]->additionalSearchParams.userParams[i].value.empty())
	{
			values.push_back("true");
	}
	else
	{
			values.push_back(sip[0]->additionalSearchParams.userParams[i].value);
	}
    }

    Rcpp::List res(names.size());

    for (size_t i = 0; i < names.size(); i++)
    {
		if (isNumber(values[i]))
		{
			res[i] = Rcpp::wrap(lexical_cast<double>(values[i]));
		}
		else if (isBool(values[i]))
		{
			res[i] = Rcpp::wrap(toBool(values[i]));
		}
		else
		{
			res[i] = Rcpp::wrap(values[i]);
		}
	}

    res.attr("names") = names;
    return res;
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
