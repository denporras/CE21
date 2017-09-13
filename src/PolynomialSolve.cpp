//============================================================================
// Name        : Raices de Polinomios
// Authors     : Dennis Gerardo Porras Barrantes, David Eduardo Gomez Vargas and Kelvin Stanley Alfaro Vega
// Version     : 1.0
// Copyright   : Assignment for the course Numerical Analysis of the Costa Rica Institute of Tecnology
// Description : Program that solve polynomial by numerical methods.
//============================================================================



#ifndef POLYNOMIALSOLVE_H_
#define POLYNOMIALSOLVE_H_

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <iostream>
#include <boost/math/tools/polynomial.hpp>
#include <complex>
#include "Methods/LaGuerreMethod.h"
#include "Methods/MullerMethod.h"
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;
using namespace boost::math;
using namespace boost::math::tools;	// for polynomial
using boost::lexical_cast;

/**
 * @breif sign print
 * @param x
 * @return
 */
template<typename T>
string sign_str(T const &x) {
	return x < 0 ? "-" : "+";
}

/**
 * @breif print inner coefficient
 * @param x
 * @return
 */
template<typename T>
string inner_coefficient(T const &x) {
	string result(" " + sign_str(x) + " ");
	if (abs(x) != T(1))
		result += lexical_cast<string>(abs(x));
	return result;
}

/**
 * @brief print the polynomial
 * @param a polynomial
 * @return a string with the poly.
 */
template<typename T>
string formula_format(polynomial<T> const &a) {
	string result;
	if (a.size() == 0)
		result += lexical_cast<string>(T(0));
	else {
		// First one is a special case as it may need unary negate.
		unsigned i = a.size() - 1;
		if (a[i] < 0)
			result += "-";
		if (abs(a[i]) != T(1))
			result += lexical_cast<string>(abs(a[i]));

		if (i > 0) {
			result += "x";
			if (i > 1) {
				result += "^" + lexical_cast<string>(i);
				i--;
				for (; i != 1; i--)
					result += inner_coefficient(a[i]) + "x^"
							+ lexical_cast<string>(i);

				result += inner_coefficient(a[i]) + "x";
			}
			i--;

			result += " " + sign_str(a[i]) + " "
					+ lexical_cast<string>(abs(a[i]));
		}
	}
	return result;
} // string formula_format(polynomial<T> const &a)

/**
 * @brief prints the result
 * @param x an array of T numbers
 * @param degree of the polynomial
 */
template<typename T>
void printResult(T x[], int degree) {
	for (int i = 0; i < degree; i++) {
		cout << "x_" << i + 1 << " = " << x[i].real();
		if (x[i].imag() >= 0)
			cout << " + ";
		else
			cout << " ";
		cout << x[i].imag() << "j" << endl;
	}
}

/**
 * @brief Template class to call the methods.
 * @param d: degree of the polynomial.
 * @param polish: bool flag for polish val.
 */
template<class T>
void setUp(int d, bool polish) {

	if (d < 10 and d>0) {
		boost::array<T, 10> arr1 = { };
		for (int i = 0; i <= d; ++i) {
			cout << "Write a number for the degree " << i << endl;
			cin >> arr1[i];
		}
		polynomial<T> pol(arr1.begin(), arr1.end());
		//Metodo de LaGuerre
		LaGuerreMethod<T> *obj2;
		complex<T> x[10];
		obj2->solvePolynomial(pol, x, polish);

		//Metodo de Muller
		MullerMethod<T> *obj;
		complex<T> * roots = new complex<T> [10];
		roots = obj->solvePolynomial(pol, 5, 0.5, polish);

		//Resultado
		cout << "Polinomio: " << formula_format(pol) << endl << "\nSoluciones:"
				<< endl;
		cout << "Método de LaGuerre: " << endl;
		printResult(x, pol.degree());
		cout << "Método de Muller" << endl;
		printResult(roots, pol.degree());
	} else {
		cout << "The degree of the polynomial is very large" << endl;
	}
}

/**
 * @brief Main function.
 * @return 0 if everything is ok.
 */
int main() {
	int p, d;
	string po;
	cout << "Choose a number (1,2) to select the precision:" << endl
			<< "1. Float" << endl << "2. Double" << endl;
	cin >> p;
	cout << "Write a number to select the degree of the polynomial:" << endl;
	cin >> d;
	cout << "Choose a letter (y,n) to select if the polishing of roots is used:"
			<< endl << "y. Yes" << endl << "n. No" << endl;
	cin >> po;
	bool polish = false;
	if (po == "y"){
		polish = true;
	}
	if (p == 1) {
		setUp<float>(d, polish);
	} else if (p == 2) {
		setUp<double>(d, polish);
	} else {
		cout << "Invalid precision type" << endl;
	}
	return 0;
}

#endif /* POLYNOMIALSOLVE_H_ */

