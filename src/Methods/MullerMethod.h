/**
 * @file MullerMethod.h
 * @brief Class that implements the Muller method to find the roots of a polynomial.
 * @author Kevin Alfaro
 * @date 11 de sept. de 2017
 */

#ifndef METHODS_MULLERMETHOD_H_
#define METHODS_MULLERMETHOD_H_

#define MAX_ITERATION 80
#define EPS_POW -3
#include <complex>
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

#include "PolynomialDeflaction.h"

template<typename T>
class MullerMethod {

private:
	complex<T> getRoot(polynomial<T> &poly, T xr, T h);
	complex<T> evaluatePolynomial(polynomial<T> &poly, complex<T> x);

public:
	MullerMethod();
	complex<T>* solvePolynomial(polynomial<T> &poly, T xr, T h,
			const bool &polish);
};

/**
 * @brief Constructor by default.
 */
template<typename T>
MullerMethod<T>::MullerMethod() {
}

/**
 * @brief Applies the Muller method to find all roots in the polynomial.
 * @param poly: Polynomial.
 * @param xr: Point of the parabole
 * @param h: h to get the other two points of the parabole.
 * @param polish: Activate or deactive root polish
 */
template<typename T>
complex<T>* MullerMethod<T>::solvePolynomial(polynomial<T> &poly, T xr, T h,
		const bool &polish) {

	const T EPS = pow(10, EPS_POW); //Estimated fractional roundoff error
	polynomial<T> temp_poly = poly;
	polynomial<T> aux_poly = poly;
	complex<T> * roots = new complex<T> [poly.degree()];

	for (unsigned int i = 0; i <= poly.degree() - 1; i++) { //Loop to find every root.
		roots[i] = this->getRoot(temp_poly, xr, h);

		if (polish) { //Polish the roots
			roots[i] = this->getRoot(temp_poly, roots[i].real(), h);
		}

		PolynomialDeflaction<T, 10> *pd; //Deflaction process to remove the root

		//Choose the correct deflate function depending if is a real or complex number
		if (abs(imag(roots[i])) <= (T(2) * EPS * abs(real(roots[i])))) {
			temp_poly = pd->deflate(temp_poly, real(roots[i]), aux_poly);
		} else {
			i++;
			std::cout << "MULLER: A complex root has been found!" << std::endl;
			roots[i] = complex<T>(real(roots[i - 1]), -1 * imag(roots[i - 1])); //Add the conjugate to roots
			temp_poly = pd->deflate2(temp_poly, roots[i - 1], aux_poly);
		}

	}

	return roots;
}

/**
 * @brief Applies the Muller method to find a root in the polynomial.
 * @param poly: Polynomial.
 * @param xr: Point of the parabole
 * @param h: h to get the other two points of the parabole.
 */
template<typename T>
complex<T> MullerMethod<T>::getRoot(polynomial<T> &poly, T xr, T h) {
	const T EPS = pow(10, EPS_POW); //Estimated fractional roundoff error

	//Reference x in complex format
	complex<T> x3 = complex<T>(xr, T(0));

	//Initial points of the parable
	complex<T> x2 = x3;
	complex<T> x1 = x3 + h;
	complex<T> x0 = x3 - h;
	//h and d coeficients
	complex<T> h0 = 0;
	complex<T> h1 = 0;
	complex<T> d0 = 0;
	complex<T> d1 = 0;
	//a b c coeficients
	complex<T> a = 0;
	complex<T> b = 0;
	complex<T> c = 0;

	//Aux variables
	complex<T> discriminant = 0;
	complex<T> rad = 0; //Radical used in the general equation
	complex<T> den = 0; //Denominator of the general eqution (b+rad or b-rad)
	complex<T> dx3 = std::numeric_limits<T>::max(); //Complete general ecuation, also equal to x3-x2
	int i = 0; //Iteration counter

	while ((abs(dx3) > abs(EPS * x3)) || (i > MAX_ITERATION)) {
		i++;
		h0 = x1 - x0;
		h1 = x2 - x1;
		d0 = (this->evaluatePolynomial(poly, x1)
				- this->evaluatePolynomial(poly, x0)) / h0; //(f(x1)-f(x0))/h0
		d1 = (this->evaluatePolynomial(poly, x2)
				- this->evaluatePolynomial(poly, x1)) / h1; //(f(x2)-f(x1))/h1
		a = (d1 - d0) / (h1 + h0);
		b = a * h1 + d1;
		c = this->evaluatePolynomial(poly, x2); //f(x2)

		discriminant = b * b - T(4) * a * c;
		rad = std::sqrt(discriminant);

		//Choose solution sign denominator
		if (abs(b + rad) > abs(b - rad)) {
			den = b + rad;
		} else {
			den = b - rad;
		}

		// Calculate the next (x3) root
		dx3 = (T(-2) * c) / den;
		x3 = x2 + dx3;

		// Secuential method
		x0 = x1;
		x1 = x2;
		x2 = x3;
	}

	//Check if the imaginary part is too little
	if (abs(imag(x3)) <= (T(2) * EPS * abs(real(x3)))) {
		x3 = complex<T>(real(x3), T(0));
	}

	//Check if the real part is too little
	if (abs(real(x3)) <= (T(2) * EPS * abs(imag(x3)))) {
		x3 = complex<T>(T(0), imag(x3));
	}

	return x3;
}

/**
 * @brief Evaluate a polynomial function as f(x)
 * @param poly: Polynomial.
 * @param x: Number to evaluate the function
 */
template<typename T>
complex<T> MullerMethod<T>::evaluatePolynomial(polynomial<T> &poly,
		complex<T> x) {
	complex<T> result = 0;

	for (int i = 0; i <= poly.degree(); i++) {
		result = result + poly[i] * pow(x, i); // coef*x^i
	}
	return result;
}

#endif /* METHODS_MULLERMETHOD_H_ */
