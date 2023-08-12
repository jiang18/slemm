#ifndef SLEMM_H_
#define SLEMM_H_

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifndef EIGEN_MKL_LIB
#define EIGEN_MKL_LIB
#ifndef EIGEN_LIB
#define EIGEN_LIB
#include "Eigen/Dense"
#endif
#include <mkl.h>
#endif

#include <omp.h>

#include "StrFunc.h"

using namespace Eigen;

inline bool file_check (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}
// input data
std::map<std::string, std::pair<float, float> > read_phenotype_file(std::string phenotype_file, std::string trait_name, std::string error_variance_weight_name);
std::map<std::string, VectorXf> read_covariate_file(std::string covariate_file, std::vector<std::string>& covariate_names);
void read_marker_info_file(
	std::string marker_info_file, 
	std::string group_header, 
	std::string weight_name, 
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight );
void get_dim_from_plink(
	const std::string binary_genotype_file_prefix, 
	const std::map<std::string, std::pair<float, float> >& indi2pheno_weight, 
	const std::map<std::string, VectorXf >& indi2covar, 
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight, 
	unsigned& indi_keep_num,
	unsigned& marker_num, 
	unsigned& group_keep_num );
void read_plink_bin_geno_file_and_generate_pxvb_data(
	const std::string binary_genotype_file_prefix, 
	std::map<std::string, std::pair<float, float> >& indi2pheno_weight, 
	std::map<std::string, VectorXf >& indi2covar, 
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight, 
	const std::pair<float, float>& beta_pdf_params,
	const float min_maf,
	std::map<std::string, float>& marker2maf,
	const float min_hwep,
	const int hwe_midp,
	std::map<std::string, float>& marker2hwep,
	std::map<std::string, std::pair<std::string, int> >& marker2pos,
	std::map<std::string, std::pair<std::string, std::string> >& marker2alleles,
	std::vector<std::string>& indi_keep,
	std::vector<std::string>& marker_keep,
	std::vector<std::string>& group_keep,
	Ref<VectorXi> group_sindex,
	Ref<VectorXi> group_ncols,
	std::vector<VectorXf>& gdiag,
	Ref<MatrixXf> xmat,
	std::vector<std::vector<bool> >& kmat,
	Ref<VectorXf> yvec,
	Ref<VectorXf> rvec);
// marker weights
void cal_maf_from_genobits(const std::vector<std::vector<bool> >& geno, Ref<VectorXf> maf);
void update_weight_by_maf (
	const std::pair<float, float>& beta_pdf_params,
	const std::map<std::string, float>& marker2maf,
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight );
// prediction of total genetic values
void cal_dgv(
	const std::string marker_estimate_file,
	const std::string binary_genotype_file_prefix,
	const std::string output_file_prefix);
// output data
void write_reml_into_files(
	const std::string output_file_prefix,
	const std::vector<std::string>& covar_names,
	const std::vector<std::string>& indi_keep,
	const Ref<VectorXf> rvec,
	const Ref<MatrixXf> xmat,
	std::map<std::string, float>& marker2maf,
	std::map<std::string, float>& marker2hwep,
	std::map<std::string, std::pair<std::string, int> >& marker2pos,
	std::map<std::string, std::pair<std::string, std::string> >& marker2alleles,
	const Ref<VectorXi> ncols,
	const std::vector<std::string>& group_keep,
	const std::vector<std::string>& marker_keep,
	const Ref<VectorXf> Py,
	const Ref<VectorXf> snp_blup,
	const Ref<VectorXf> blue,
	const Ref<MatrixXf> var_blue,
	const float& vg,
	const float& ve,
	const float& lrt );

typedef struct {int vi; float xPy; float xPx; float xBx; float xSx;} gstat;
void write_gstat_into_files(
	const std::string output_file_prefix,
	const int pop_size,
	const int window_size,
	const float vg,
	const std::map<std::string, float>& marker2maf,
	const std::map<std::string, float>& marker2hwep,
	std::map<std::string, std::pair<std::string, int> >& marker2pos,
	const std::vector<std::string>& marker_keep,
	const std::vector<gstat>& lmm_out );

#endif
