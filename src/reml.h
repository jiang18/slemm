#ifndef REML_H_
#define REML_H_
#include "slemm.h"

void SLDF_REML(
	const Ref<VectorXf>& rvec,  // error variance weight R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<MatrixXf>& X,     // covariate matrix
	const Ref<VectorXf>& y,     // phenotype vector
	const Ref<VectorXf>& persnp,// per-SNP heritability weight
	const int subset_size,      // 
	const int rng_seed,         // random seed
	const int n_V,              // number of random probes
	const float tol_L,          // relative lanczos tolerance
	Ref<VectorXf> Py,
	Ref<VectorXf> snp_blup,
	Ref<VectorXf> blue,
	float& vg,
	float& ve,
	float& llr,
	const float s2max = .7000,  // maximal h2 value
	const float s2min = .0001,  // minimal h2 value
	const float tol_VC = 1e-4,  // abs. tolerance for h2 estimation
	const bool verbose = true,  // verbose output
	const int p_freq = 5        // print frequency
);

void LMC_REML(
	const Ref<VectorXf>& rvec,  // error variance weight R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<MatrixXf>& X,     // covariate matrix
	const Ref<VectorXf>& y,     // phenotype vector
	const Ref<VectorXf>& persnp,// per-SNP heritability weight
	const int subset_size,      // 
	const int rng_seed,         // random seed
	const int n_MC,             // number of MC samples
	const float tol_L,          // relative lanczos tolerance
	const bool do_lmm,
	const int num_qf_markers,
	const bool fake_geno,
	const int window_size,
	Ref<VectorXf> Py,
	Ref<VectorXf> snp_blup,
	Ref<VectorXf> blue,
	float& vg,
	float& ve,
	std::vector<gstat>& lmm_out,
	const float s2max = .7000,  // maximal h2 value
	const float s2min = .0001,  // minimal h2 value
	const bool verbose = true,  // verbose output
	const int p_freq = 5        // print frequency
);

#endif
