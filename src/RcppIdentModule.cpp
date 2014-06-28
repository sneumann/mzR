#include <Rcpp.h>
#include "RcppIdent.h"


RCPP_MODULE(Ident){
	
  using namespace Rcpp;

  class_<RcppIdent>( "Ident" )
    .constructor("Initialises a new Rccp ident object.")
    .method( "open", &RcppPwiz::open, "Opens a mass spec file (mzXML, mzData, etc.) and creates a pwiz object" )
    ;
}
