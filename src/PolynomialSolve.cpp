#ifndef POLYNOMIALSOLVE_H_
#define POLYNOMIALSOLVE_H_
#include <iostream>
#include <boost/math/tools/polynomial.hpp>
#include <complex>
#include "Methods/PolynomialDeflaction.h"


using namespace std;
using namespace boost::math;
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;


int main() {
	polynomial<double> b{{-2.0, 1.0}};
	PolynomialDeflaction<double, 10> *pd;

	//Prueba para reacionales
	boost::array<double, 10> const arr1 = {{-16, 2, 0, -8, 1}};
	polynomial<double> const pol(arr1.begin(), arr1.end());
	pd->deflate(pol,8,b);

	//Prueba para irracionales
	boost::array<double, 10> const arr2 = {{13, 59, -16, -1, 1}};
	polynomial<double> const pol2(arr2.begin(), arr2.end());
	complex<double> root(3,2);
	pd->deflate2(pol2,root,b);
	return 0;
}



#endif /* POLYNOMIALSOLVE_H_ */

