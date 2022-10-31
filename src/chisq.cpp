#include "chisq.h"
#include <iostream>

/*
	Adapted from chisqsum.R in the R package survey.
*/

double k0 (const double zeta, const Ref<const VectorXd>& lambda)
{
	return(-(1 - 2*zeta*lambda.array()).log().sum()/2.);
}

double kprime0 (const double zeta, const Ref<const VectorXd>& lambda)
{
	return( (lambda.array()/(1 - 2*zeta*lambda.array())).sum() );
}

double kpprime0 (const double zeta, const Ref<const VectorXd>& lambda)
{
	return( 2.*(lambda.array()/(1 - 2*zeta*lambda.array())).square().sum() );
}

inline double sign (double x)
{
	double value;
	if( x < 0.0) value = -1.0;
	else value = 1.0;
	return value;
}

// Use bisection method to find root of kprime0 - x =0.
// To implement Brent's method soon.
double uniroot(double a, double b, const Ref<const VectorXd>& lambda, const double x, const double tol)
{
	double dx;
	double r;
	double fr;
	double fa;
	double fb;
	
	fa = kprime0(a, lambda) - x;
	fb = kprime0(b, lambda) - x;
	if(sign(fa) == sign(fb)) {
		std::cerr << "\nError: f(a) and f(b) have the same sign.\n";
      	exit ( 1 );
	}
	
	while(1) {
		r = (a + b)/2.;
		fr = kprime0(r, lambda) - x;
		
		dx = 0.5 * std::abs( b-a );
		if(dx <= tol) break;
		
		if(sign(fr) == sign(fa)) {
			a = r;
			fa = fr;
		} else {
			b = r;
			fb = fr;
		}
	}
	return r;
}

double pchisqsum_saddle(double x, VectorXd lambda)
{
	// all lambda's must be >= 0
	if( (lambda.array() < 0).any() ) {
		std::cerr << "\nAll coefficients must be >= 0 in the weighted sum of chisq distributions.\n";
		exit ( 1 );
	}
	double pval;
	const double tol = 1e-8;
	
	double d = lambda.maxCoeff();
	lambda /= d;
	x /= d;
	int n = lambda.size();

	// special treatment for one chisq
	if(n == 1) {
		return getOneDfChisqPval(x);
	}
	
	double lmin;
	if(x > lambda.sum())
		lmin = -0.01;
	else
		lmin = -n/(2*x);
	double lmax = 1/2. * 0.99999;
	
	double hatzeta = uniroot(lmin, lmax, lambda, x, tol);
	double w = sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta,lambda)));
	double v = hatzeta*sqrt(kpprime0(hatzeta,lambda));
	
	if (std::abs(hatzeta)<1e-4) {
		double tr = lambda.mean();
		double tr2 = lambda.array().square().mean()/(tr*tr);
		double scale = tr*tr2;
		double df = n/tr2;
		pval = pchisq(df, x/scale);
	} else {
		pval = getStNormalPval( w+log(v/w)/w );
	}
	return pval;
}
