#include <Rcpp.h>
#include "RcppIdent.h"


RCPP_MODULE(Ident){
	
  using namespace Rcpp;

  class_<RcppIdent>( "Ident" )
    .constructor("Initialises a new Rccp ident object.")
    .method( "open", &RcppIdent::open, "Opens a mass spec file (mzXML, mzData, etc.) and creates a pwiz object" )
    .method( "getCreationDate", &RcppIdent::getCreationDate, "Opens a mass spec file (mzXML, mzData, etc.) and creates a pwiz object" )
    ;
}
