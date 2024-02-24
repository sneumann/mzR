#include <Rcpp.h>
#include "RcppIdent.h"


RCPP_MODULE(Ident)
{

    using namespace Rcpp;

    class_<RcppIdent>( "Ident" )
    .constructor("Initialises a new Rccp ident object.")
    .method( "open", &RcppIdent::open, "Opens a mass spec file (mzXML, mzML, etc.) and creates a pwiz object" )
    .method( "getIDInfo", &RcppIdent::getIDInfo, "Basic information about this mzid files" )
    .method( "getPsmInfo", &RcppIdent::getPsmInfo, "Basic information about this mzid files" )
    .method( "getModInfo", &RcppIdent::getModInfo, "Modification information about this mzid files" )
    .method( "getSubInfo", &RcppIdent::getSubInfo, "Substitution information about this mzid files" )
    .method( "getScore", &RcppIdent::getScore, "Scoring information about this mzid files" )
    .method( "getPara", &RcppIdent::getPara, "Parameters used in identification." )
    .method( "getDB", &RcppIdent::getDB, "Database used in identification." )
    .method( "getSpecParams", &RcppIdent::getSpecParams, "SpectrumIdentificationResult cvParams" )
    ;
}
