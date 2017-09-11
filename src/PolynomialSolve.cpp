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
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;

//] [/polynomial_arithmetic_1]

template <typename T>
string sign_str(T const &x)
{
	return x < 0 ? "-" : "+";
}

template <typename T>
string inner_coefficient(T const &x)
{
	string result(" " + sign_str(x) + " ");
	if (abs(x) != T(1))
		result += lexical_cast<string>(abs(x));
	return result;
}

template <typename T>
string formula_format(polynomial<T> const &a)
{
	string result;
	if (a.size() == 0)
		result += lexical_cast<string>(T(0));
	else
	{
		// First one is a special case as it may need unary negate.
		unsigned i = a.size() - 1;
		if (a[i] < 0)
			result += "-";
		if (abs(a[i]) != T(1))
			result += lexical_cast<string>(abs(a[i]));

		if (i > 0)
		{
			result += "x";
			if (i > 1)
			{
				result += "^" + lexical_cast<string>(i);
				i--;
				for (; i != 1; i--)
					result += inner_coefficient(a[i]) + "x^" + lexical_cast<string>(i);

				result += inner_coefficient(a[i]) + "x";
			}
			i--;

			result += " " + sign_str(a[i]) + " " + lexical_cast<string>(abs(a[i]));
		}
	}
	return result;
} // string formula_format(polynomial<T> const &a)

int main() {
/**
	polynomial<double> b{{-1}};
	PolynomialDeflaction<double, 10> *pd;
	//Prueba para irracionales
	boost::array<double, 10> const arr2 = {{-4, 5, -1, 1}};
	polynomial<double> const pol2(arr2.begin(), arr2.end());
	complex<double> root(0,2);
	pd->deflate2(pol2,root,b);
	//Prueba para reacionales
	boost::array<double, 10> const arr1 = {{-16, 2, 0, -8, 1}};
	polynomial<double> const pol(arr1.begin(), arr1.end());
	pd->deflate(pol,8,b);
*/



	boost::array<double, 10> const arr1 = {{-8,-4,2,-1,1}}; //(X-2)(x^2+4)(x+1)
	polynomial<double> pol(arr1.begin(), arr1.end());

	//Metodo de Muller
	MullerMethod<double> *obj;
	complex<double> * roots = new complex<double>[pol.degree()];
	roots = obj->solvePolynomial(pol,5,0.5,true);

	for(unsigned int i = 0; i <= pol.degree()-1; i++){
		cout << "Real: " << roots[i].real() <<"\n";
		cout << "Imaginario: " << roots[i].imag() <<"\n";

	}


/**
	//Creando el polinomio
	//boost::array<double, 10> const arr1 = {{1, 2, 3}};
	//polynomial<double> pol(arr1.begin(), arr1.end());

	//Metodo de LaGuerre
	LaGuerreMethod<double> *obj2;
	const bool polish = true;
	complex<double> x[7];
	obj2->solvePolynomial(pol, x, polish);



	//Resultado
	cout << "Polinomio: " << formula_format(pol) << endl << "Soluciones:" << endl;
	for(int i=0; i < pol.degree(); i++){
		cout << "x_" << i+1 << " = " << x[i].real();
		if(x[i].imag()>=0)
			cout << " + ";
		else
			cout << " ";
		cout << x[i].imag() << "j" << endl;
	}

*/

	return 0;
}



#endif /* POLYNOMIALSOLVE_H_ */

