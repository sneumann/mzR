//
// Based on code sent by Kevin Ushey 
// to Rcpp-devel on Tue, 8 Jul 2014
//



#include <Rcpp.h>
using namespace Rcpp;

class ListBuilder {

public:

  ListBuilder() {};
  ~ListBuilder() {};

  inline ListBuilder& add(std::string name, SEXP x) {
    names.push_back(name);
    elements.push_back(x);
    return *this;
  }

  inline List get() const {
    return static_cast<List>(*this);
  }

  inline operator List() const {
    List result(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
      result[i] = elements[i];
    }
    result.attr("names") = wrap(names);
    return result;
  }

  inline operator DataFrame() const {
    List result = static_cast<List>(*this);
    result.attr("class") = "data.frame";
    result.attr("row.names") = IntegerVector::create(NA_INTEGER, XLENGTH(elements[0]));
    return result;
  }

private:

  std::vector<std::string> names;
  std::vector<SEXP> elements;

  ListBuilder(ListBuilder const&) {};

};

// [[Rcpp::export]]
DataFrame test_builder(SEXP x, SEXP y, SEXP z) {
  return ListBuilder()
    .add("foo", x)
    .add("bar", y)
    .add("baz", z);
}
