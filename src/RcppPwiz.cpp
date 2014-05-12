#include "RcppPwiz.h"

void RcppPwiz::open(const string& fileName) {

	filename = Rcpp::StringVector::create( fileName );
	msd = new MSDataFile(fileName);

}


Rcpp::StringVector RcppPwiz::getFilename (  ) {

  return filename;
}

int RcppPwiz::getLastScan() const {
  if (msd != NULL) {
	  SpectrumListPtr slp = msd->run.spectrumListPtr;
	  return slp->size();
  }
  Rprintf("Warning: Ramp not yet initialized.\n ");
  return -1;
}
