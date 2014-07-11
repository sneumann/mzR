#include <Rcpp.h>
#include "RcppIdent.h"


RCPP_MODULE(Ident){
	
  using namespace Rcpp;

  class_<RcppIdent>( "Ident" )
    .constructor("Initialises a new Rccp ident object.")
    .method( "open", &RcppIdent::open, "Opens a mass spec file (mzXML, mzData, etc.) and creates a pwiz object" )
    .method( "getIDInfo", &RcppIdent::getIDInfo, "Basic information about this mzid files" )
    .method( "getPepInfo", &RcppIdent::getPepInfo, "Basic information about this mzid files" )
    .method( "getModInfo", &RcppIdent::getModInfo, "Modification information about this mzid files" )
    ;
}
