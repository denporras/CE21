/*
 * PolynomialDeflaction.h
 *
 *  Created on: 9 de set. de 2017
 *      Author: daedgomez
 */

#ifndef METHODS_POLYNOMIALDEFLACTION_H_
#define METHODS_POLYNOMIALDEFLACTION_H_

#include <boost/math/tools/polynomial.hpp>
#include <iostream> //Print values
#include <cmath>	//Math operators
#include <complex>


using namespace std;
using namespace boost::math;
using namespace boost::math::tools; // for polynomial
using boost::lexical_cast;

template<typename T, size_t N>
class PolynomialDeflaction {
public:
	PolynomialDeflaction();
	polynomial<T> deflate(const polynomial<T>& poly, const T& root,const polynomial<T>& residuo);
	polynomial<T> deflate2(const polynomial<T>& poly, const
			complex<T>& root,const polynomial<T>& residuo);
private:
	string sign_str(T const &x);
	string inner_coefficient(T const &x);
	string formula_format(polynomial<T> const &a);
};

template<typename T, size_t N>
PolynomialDeflaction<T, N>::PolynomialDeflaction() {
}

template<typename T, size_t N>
polynomial<T> PolynomialDeflaction<T, N>::deflate(const polynomial<T>& poly, const T& root,const polynomial<T>& residuo){
	int n = poly.size();
	boost::array<T, N> sol = {};
	T tmp[n];
	tmp[n-1] = poly[n-1];
	for(int i=n-2; i>=0; i--){
	      tmp[i] = poly[i] + tmp[i+1] * root;
	}
	for (int i=n-1; i>0; i--){
		sol[i-1] = tmp[i];
	}
	boost::array<T, 1> res = {};
	res[0] = tmp[0];
	polynomial<T> d3a(res.begin(), res.end());
	cout << "Residuo = " << formula_format(d3a) << "\n";

	polynomial<T> poly_sol(sol.begin(), sol.end());
	cout << "Cociente = " << formula_format(poly_sol) << "\n";
	return poly_sol;
}


template<typename T, size_t N>
polynomial<T> PolynomialDeflaction<T, N>::deflate2(const polynomial<T>& poly, const
		complex<T>& root,const polynomial<T>& residuo){
	int n = poly.size();
	boost::array<T, N> sol = {};
	complex<T> tmp[n];
	tmp[n-1] = poly[n-1];
	for(int i=n-2; i>=0; i--){
		tmp[i] = poly[i] + tmp[i+1] * root;
	}
	complex<T> tmp2[n];
	tmp2[n-1] = tmp[n-1];
	complex<T> root2 = conj(root);
	for(int i=n-2; i>=0; i--){
		tmp2[i] = tmp[i] + tmp2[i+1] * root2;
		}
	for (int i=n-1; i>0; i--){
		sol[i-1] = real(tmp2[i]);
	}
	boost::array<T, 1> res = {};
	res[0] = real(tmp2[0]);
	polynomial<T> d3a(res.begin(), res.end());
	cout << "Residuo = " << formula_format(d3a) << "\n";
	polynomial<T> poly_sol(sol.begin(), sol.end());
	cout << "Cociente = " << formula_format(poly_sol) << "\n";
	return poly_sol;
}


template <typename T, size_t N>
string PolynomialDeflaction<T, N>::sign_str(T const &x)
{
  return x < 0 ? "-" : "+";
}

template <typename T, size_t N>
string PolynomialDeflaction<T, N>::inner_coefficient(T const &x)
{
  string result(" " + sign_str(x) + " ");
  if (abs(x) != T(1))
      result += lexical_cast<string>(abs(x));
  return result;
}

/*! Output in formula format.
For example: from a polynomial in Boost container storage  [ 10, -6, -4, 3 ]
show as human-friendly formula notation: 3x^3 - 4x^2 - 6x + 10.
*/
template <typename T, size_t N>
string PolynomialDeflaction<T, N>::formula_format(polynomial<T> const &a)
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

    	if (i==0)
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



#endif /* METHODS_POLYNOMIALDEFLACTION_H_ */
