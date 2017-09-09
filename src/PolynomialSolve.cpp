#ifndef POLYNOMIALSOLVE_H_
#define POLYNOMIALSOLVE_H_
#include <iostream>
#include <boost/math/tools/polynomial.hpp>

#include "Methods/PolynomialDeflaction.h"


using namespace std;
using namespace boost::math;
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;


int main() {
	boost::array<double, 5> const d3a = {{-16, 2, 0, -8, 1}};
	polynomial<double> const a(d3a.begin(), d3a.end());
	polynomial<double> b{{-2.0, 1.0}};
	PolynomialDeflaction<double, 4> *pd;
	pd->deflate(a,-8,b);
	return 0;
}



#endif /* POLYNOMIALSOLVE_H_ */

