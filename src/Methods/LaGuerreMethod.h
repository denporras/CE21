/**
 * @file LaGuerreMethod.h
 * @brief Class that implements the LaGuerre method to find the roots of a polynomial.
 * @author Dennis Porras
 * @date 9 de sept. de 2017
 */


#ifndef METHODS_LAGUERREMETHOD_H_
#define METHODS_LAGUERREMETHOD_H_

#define MAX_ITERATION 80
#define MT 10
#define MR 8
#define MAX(x, y) (((x) > (y)) ? (x) : (y))


#include <complex>
#include <cmath>
#include <limits>

#include "PolynomialDeflaction.h"

template<typename T>
class LaGuerreMethod {
public:
	LaGuerreMethod();
	void solvePolynomial(polynomial<T> &a, complex<T> * roots, const bool &polish);
private:
	void process(polynomial<T> &a, complex<T> &x, int &iteration);
};

/**
 * @brief Constructor by default.
 */
template<typename T>
LaGuerreMethod<T>::LaGuerreMethod(){ }

/**
 * @brief Applies the LaGuerre method to find a root in the polynomial.
 * @param a: Polynomial.
 * @param x: An x to find or improve.
 * @param iteration: Amount of iterations.
 */
template<typename T>
void LaGuerreMethod<T>::process(polynomial<T> &a, complex<T> &x, int &iteration){
	const T frac[MR + 1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0}; //Fractions used to break a limit cycle.
	const T er_stop = numeric_limits<T>::epsilon(); //It's the estimated fractional roundoff error.
	int m = a.degree();
	complex<T> b, d, f, g, g2, h, sq, gp, gm, dx, x1; //LagGuerre's method uses this variable to compute the root.
	for(int i = 1; i < MAX_ITERATION; i++){ //Loop over iterations up to allowed maximum.
		iteration = i;
		b = a[m];
		T error = abs(b);
		d = f = T(0);
		T abx = abs(x);
		for(int j = m - 1; j >= 0; j--){ //Calculate the second derivative.
			f = x*f + d;
			d = x*d + b;
			b = x*b + a[j];
			error = abs(b) + abx * error;
		}
		error *= er_stop; //Estimate of roundoff error in evaluating polynomial
		if(abs(b) <= error) //It's on the root.
			return ;
		g = d/b; //Start LaGuerre Method.
		g2 = g*g;
		h = g2 - (T(2)*f/b);
		sq = sqrt(T(m-1)*(T(m)*h-g2));
		gp = g + sq;
		gm = g - sq;
		T abp = abs(gp);
		T abm = abs(gm);
		if(abp < abm)
			gp = gm;
		dx = MAX(abp,abm) > T(0) ? T(m)/gp : polar(1+abx, T(i));
		x1 = x - dx;
		if(x == x1) //Converged
			return ;
		if((i % 10) != 0)
			x = x1;
		else
			x -= frac[i/10]*dx; //Increase posibilities.
	}
}

/**
 * @brief Finds all the root in a polynomial applying the LaGuerre method and also deflacts the polinomial.
 * @param a: Polynomial.
 * @param roots: Array of complex roots.
 * @param polish: Boolean if polish.
 */
template<typename T>
void LaGuerreMethod<T>::solvePolynomial(polynomial<T> &a, complex<T> * roots, const bool &polish){
	const T EPS = pow(10,-14); //Estimated fractional roundoff error.
	int i, iterations;
	complex<T> x;
	int m = a.degree();
	polynomial<T> ad = a; //Copy the polynomial.
	for(int j = m-1; j >= 0; j--){ //Loop to find every root.
		x = T(0);//Zero to improve convergence.
		polynomial<T> ad_v = ad;
		this->process(ad_v, x, iterations);
		if(abs(imag(x)) <= (T(2)*EPS*abs(real(x))))
			x = complex<T>(real(x),T(0));
		roots[j]=x;
		PolynomialDeflaction<T, 10> *pd; //Se aplica deflacción
		ad = pd->deflate2(ad,x,ad);
	}

	if (polish){ //Polish the roots
		for (int j=0;j<m;j++){
			this->process(a,roots[j],iterations);
		}
	}
	for (int j=1;j<m;j++) { //Order the roots by the real part.
		x=roots[j];
		for (i=j-1;i>=0;i--) {
			if (real(roots[i]) <= real(x))
				break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}

	for(int i = 1; i <= m-1; i++){ //The method does not verify if the root are equal, which means the imaginary part changes sign
		if(abs(roots[i]-roots[i-1]) < numeric_limits<T>::epsilon()){
			T tmp = roots[i].imag() * -1;
			roots[i].imag(tmp);
		}
	}
}



#endif /* METHODS_LAGUERREMETHOD_H_ */
