/**
 * @file PolynomialDeflaction.h
 * @brief Class that implements the polynomial deflaction to decrease a degree of the polynomial.
 * @author Daedgomez
 * @date 9 de sept. de 2017
 */
#ifndef METHODS_POLYNOMIALDEFLACTION_H_
#define METHODS_POLYNOMIALDEFLACTION_H_

#include <boost/math/tools/polynomial.hpp>
#include <cmath>	//Math operators
#include <complex>	//Math complex

using namespace std;
using namespace boost::math;
using namespace boost::math::tools;	// for polynomial

/**
 * @brief Class declaration.
 */
template<typename T, size_t N>
class PolynomialDeflaction {
public:
	PolynomialDeflaction();
	polynomial<T> deflate(const polynomial<T>& poly, const T& root,
			polynomial<T>& residuo);
	polynomial<T> deflate2(const polynomial<T>& poly, const complex<T>& root,
			polynomial<T>& residuo);
};
/**
 * @brief Constructor by default.
 */
template<typename T, size_t N>
PolynomialDeflaction<T, N>::PolynomialDeflaction() {
}

/**
 * @brief Deflate for real roots.
 * @param poly: A polynomial from boost::math::tools::polynomial.
 * @param root: A root that is going to be the divider.
 * @param residuo: A reference param for the residue.
 * @return a lower degree polynomial.
 */
template<typename T, size_t N>
polynomial<T> PolynomialDeflaction<T, N>::deflate(const polynomial<T>& poly,
		const T& root, polynomial<T>& residuo) {
	int n = poly.size();
	boost::array<T, N> sol = { };
	T tmp[n];
	tmp[n - 1] = poly[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		tmp[i] = poly[i] + tmp[i + 1] * root;
	}
	for (int i = n - 1; i > 0; i--) {
		sol[i - 1] = tmp[i];
	}
	residuo[0] = tmp[0];
	polynomial<T> poly_sol(sol.begin(), sol.end());
	return poly_sol;
}

/**
  * @brief Deflate for complex roots.
 * @param poly: A polynomial from boost::math::tools::polynomial.
 * @param root: A root that is going to be the divider.
 * @param residuo: A reference param for the residue.
 * @return a lower degree polynomial.
 */
template<typename T, size_t N>
polynomial<T> PolynomialDeflaction<T, N>::deflate2(const polynomial<T>& poly,
		const complex<T>& root, polynomial<T>& residuo) {
	int n = poly.size();
	boost::array<T, N> sol = { };
	complex<T> tmp[n];
	tmp[n - 1] = poly[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		tmp[i] = poly[i] + tmp[i + 1] * root;
	}
	complex<T> tmp2[n - 1];
	tmp2[n - 2] = tmp[n - 1];
	complex<T> root2 = conj(root);
	for (int i = n - 3; i >= 0; i--) {
		tmp2[i] = tmp[i + 1] + tmp2[i + 1] * root2;
	}
	for (int i = n - 2; i > 0; i--) {
		sol[i - 1] = real(tmp2[i]);
	}
	residuo[0] = real(tmp2[0]);
	polynomial<T> poly_sol(sol.begin(), sol.end());
	return poly_sol;
}

#endif /* METHODS_POLYNOMIALDEFLACTION_H_ */
