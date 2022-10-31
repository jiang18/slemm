/*  == expint: Exponential Integral and Incomplete Gamma Function ==
 *
 *  Functions to compute the incomplete gamma function
 *
 *     G(a,x) = int_x^infty t^{a-1} exp(-t) dt
 *
 *  for a real and x >= 0. [This differs from 'pgamma' of base R in
 *  that negative values of 'a' are admitted.]
 *
 *  Copyright (C) 2016 Vincent Goulet
 *
 *  The code in part IMPLEMENTATION is derived from the GNU Scientific
 *  Library (GSL) v2.2.1 <https://www.gnu.org/software/gsl/>
 *
 *  Copyright (C) 2007 Brian Gough
 *  Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
 *
 *  The code in part R TO C INTERFACE is derived from R source code.
 *
 *  Copyright (C) 1995--1997 Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2016 The R Core Team.
 *  Copyright (C) 2003--2016 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301, USA.
 *
 *  AUTHOR for the GSL: G. Jungman
 *  AUTHOR for expint: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *                     with much indirect help from the R Core Team
 */

#include "gamma_inc.h"

#include <cmath>
#include <iostream>
#include <cfloat>
#include <tr1/cmath>


/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
 *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
 * See gamma_inc_Q_CF() below.
 *
 */
long double gamma_inc_F_CF(long double a, long double x)
{
    const int    nmax  =  5000;

    long double hn = 1.0;           /* convergent */
    long double Cn = 1.0 / LDBL_TINY;
    long double Dn = 1.0;
    int n;

    /* n == 1 has a_1, b_1, b_0 independent of a,x,
       so that has been done by hand                */
    for (n = 2 ; n < nmax ; n++)
    {
	long double an;
	long double delta;

	if (E1_IS_ODD(n))
	    an = 0.5 * (n - 1)/x;
	else
	    an = (0.5 * n - a)/x;

	Dn = 1.0 + an * Dn;
	if (fabs(Dn) < LDBL_TINY)
	    Dn = LDBL_TINY;
	Cn = 1.0 + an/Cn;
	if (fabs(Cn) < LDBL_TINY)
	    Cn = LDBL_TINY;
	Dn = 1.0/Dn;
	delta = Cn * Dn;
	hn *= delta;
	if (fabsl(delta-1.0) < LDBL_EPSILON)
	    break;
    }

    if (n == nmax)
	std::cout<<"maximum number of iterations reached in gamma_inc_F_CF"<<std::endl;

    return hn;
}

/* Useful for small a and x. Handles the subtraction analytically. */
long double gamma_inc_Q_series(long double a, long double x)
{
    const int nmax = 5000;
    long double term1;  /* 1 - x^a/Gamma(a+1) */
    long double sum;    /* 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ... */
    int n;
    long double term2;  /* a temporary variable used at the end */

    {
	/* Evaluate series for 1 - x^a/Gamma(a+1), small a */
	const long double pg21 = -2.404113806319188570799476;  /* PolyGamma[2,1] */
	const long double lnx = logl(x);
	const long double el = EULER_CNST + lnx;
	const long double c1 = -el;
	const long double c2 = M_PI*M_PI/12.0 - 0.5*el*el;
	const long double c3 = el*(M_PI*M_PI/12.0 - el*el/6.0) + pg21/6.0;
	const long double c4 = -0.04166666666666666667
	    * (-1.758243446661483480 + lnx)
	    * (-0.764428657272716373 + lnx)
	    * ( 0.723980571623507657 + lnx)
	    * ( 4.107554191916823640 + lnx);
	const long double c5 = -0.0083333333333333333
	    * (-2.06563396085715900 + lnx)
	    * (-1.28459889470864700 + lnx)
	    * (-0.27583535756454143 + lnx)
	    * ( 1.33677371336239618 + lnx)
	    * ( 5.17537282427561550 + lnx);
	const long double c6 = -0.0013888888888888889
	    * (-2.30814336454783200 + lnx)
	    * (-1.65846557706987300 + lnx)
	    * (-0.88768082560020400 + lnx)
	    * ( 0.17043847751371778 + lnx)
	    * ( 1.92135970115863890 + lnx)
	    * ( 6.22578557795474900 + lnx);
	const long double c7 = -0.00019841269841269841
	    * (-2.5078657901291800 + lnx)
	    * (-1.9478900888958200 + lnx)
	    * (-1.3194837322612730 + lnx)
	    * (-0.5281322700249279 + lnx)
	    * ( 0.5913834939078759 + lnx)
	    * ( 2.4876819633378140 + lnx)
	    * ( 7.2648160783762400 + lnx);
	const long double c8 = -0.00002480158730158730
	    * (-2.677341544966400 + lnx)
	    * (-2.182810448271700 + lnx)
	    * (-1.649350342277400 + lnx)
	    * (-1.014099048290790 + lnx)
	    * (-0.191366955370652 + lnx)
	    * ( 0.995403817918724 + lnx)
	    * ( 3.041323283529310 + lnx)
	    * ( 8.295966556941250 + lnx);
	const long double c9 = -2.75573192239859e-6
	    * (-2.8243487670469080 + lnx)
	    * (-2.3798494322701120 + lnx)
	    * (-1.9143674728689960 + lnx)
	    * (-1.3814529102920370 + lnx)
	    * (-0.7294312810261694 + lnx)
	    * ( 0.1299079285269565 + lnx)
	    * ( 1.3873333251885240 + lnx)
	    * ( 3.5857258865210760 + lnx)
	    * ( 9.3214237073814600 + lnx);
	const long double c10 = -2.75573192239859e-7
	    * (-2.9540329644556910 + lnx)
	    * (-2.5491366926991850 + lnx)
	    * (-2.1348279229279880 + lnx)
	    * (-1.6741881076349450 + lnx)
	    * (-1.1325949616098420 + lnx)
	    * (-0.4590034650618494 + lnx)
	    * ( 0.4399352987435699 + lnx)
	    * ( 1.7702236517651670 + lnx)
	    * ( 4.1231539047474080 + lnx)
	    * ( 10.342627908148680 + lnx);

	term1 = a*(c1+a*(c2+a*(c3+a*(c4+a*(c5+a*(c6+a*(c7+a*(c8+a*(c9+a*c10)))))))));
    }

    {
	/* Evaluate the sum */
	long double t = 1.0;
	sum = 1.0;

	for (n = 1; n < nmax; n++)
	{
	    t *= -x/(n+1.0);
	    sum += (a+1.0)/(a+n+1.0)*t;
	    if (fabsl(t/sum) < LDBL_EPSILON)
		break;
	}
    }

    term2 = (1.0 - term1) * a/(a + 1.0) * x * sum;

    if (n == nmax)
	std::cout<<"maximum number of iterations reached in gamma_inc_F_CF"<<std::endl;

    return term1 + term2;
}

/* Adapted from kfunc.c in samtools by Heng Li.

*/
/* The following computes regularized incomplete gamma functions.
 * Formulas are taken from Wiki, with additional input from Numerical
 * Recipes in C (for modified Lentz's algorithm) and AS245
 * (http://lib.stat.cmu.edu/apstat/245).
 *
 * A good online calculator is available at:
 *
 *   http://www.danielsoper.com/statcalc/calc23.aspx
 *
 * It calculates upper incomplete gamma function, which equals
 * kf_gammaq(s,z)*tgamma(s).
 */

// regularized lower incomplete gamma function, by series expansion
static long double _kf_gammap(long double s, long double z)
{
	long double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 1000; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < LDBL_EPSILON) break;
	}
	return expl(s * logl(z) - z - lgammal(s + 1.) + logl(sum));
}
// regularized upper incomplete gamma function, by continued fraction
static long double _kf_gammaq(long double s, long double z)
{
	int j;
	long double C, D, f;
	f = 1. + z - s; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	// See Numerical Recipes in C, 2nd edition, section 5.2
	for (j = 1; j < 1000; ++j) {
		long double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
		D = b + a * D;
		if (D < LDBL_TINY) D = LDBL_TINY;
		C = b + a / C;
		if (C < LDBL_TINY) C = LDBL_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabsl(d - 1.) < LDBL_EPSILON) break;
	}
	return expl(s * logl(z) - z - lgammal(s) - logl(f));
}

long double kf_gammap(long double s, long double z)
{
	return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);
}

long double kf_gammaq(long double s, long double z)
{
	return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z);
}

long double expint_E1(long double x)
{
	if(x<100) return (-std::tr1::expint(-(double)x));  // expintl would fail for a number close to 1
	else return (-std::tr1::expintl(-x));
}

/* Adapted from gamma_inc.c in the R package expint.

*/
long double gamma_inc(long double a, long double x)
{
    if (x < 0.0)
	return(NAN);
    else if (x == 0.0)
	return tgammal(a);
    else if (a == 0.0)
	return expint_E1(x);
    else if (a > 0.0)
	return tgammal(a) * kf_gammaq(a, x);
    else if (x > 0.25)
    {
	/* continued fraction seems to fail for x too small; otherwise
	   it is ok, independent of the value of |x/a|, because of the
	   non-oscillation in the expansion, i.e. the CF is
	   un-conditionally convergent for a < 0 and x > 0
	*/
	return expl((a - 1) * logl(x) - x) * gamma_inc_F_CF(a, x);
    }
    else if (fabsl(a) < 0.5)
    {
	return tgammal(a) * gamma_inc_Q_series(a, x);
    }
    else
    {
	/* a = fa + da; da >= 0 */
	const long double fa = floor(a);
	const long double da = a - fa;

	long double gax  = (da > 0.0 ? tgammal(da) * kf_gammaq(da, x)
		                : expint_E1(x) );
	long double alpha = da;

	/* Gamma(alpha-1,x) = 1/(alpha-1) (Gamma(a,x) - x^(alpha-1) e^-x) */
	do
	{
	    const long double shift = expl(-x + (alpha - 1.0L) * logl(x));
	    gax = (gax - shift)/(alpha - 1.0L);
	    alpha -= 1.0L;
	} while (alpha > a);

	return gax;
  }
}

