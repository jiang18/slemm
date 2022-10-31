#ifndef CHISQ_H_
#define CHISQ_H_

#include <cmath>

#ifndef EIGEN_LIB
#define EIGEN_LIB
#include "Eigen/Dense"
using namespace Eigen;
#endif

#include "gamma_inc.h"

inline double getOneDfChisqPval (const double x)
{
	if(x <= 0) return(1.0);
	return erfc(sqrt(x) * M_SQRT1_2);
}

inline double getStNormalPval (const double x)
{
	return erfc(x * M_SQRT1_2)/2.;
}

// upper tail probability
inline double pchisq(const double df, const double x)
{
	return (double)kf_gammaq(df/2.,x/2.);
}

double pchisqsum_saddle(double x, VectorXd lambda);

#endif
