#ifndef WLS_H_
#define WLS_H_
#include "slemm.h"

void weighted_least_squares(
	const Ref<VectorXf>& rvec,  // error variance weight R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<MatrixXf>& X,     // covariate matrix
	const Ref<VectorXf>& y,     // phenotype vector
	const int subset_size,      // 
	Ref<VectorXf> snp_blup
);

#endif
