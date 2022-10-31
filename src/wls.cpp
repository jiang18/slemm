#include "wls.h"
#include <tuple>

void weighted_least_squares(
	const Ref<VectorXf>& rvec,  // error variance weight R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<MatrixXf>& X,     // covariate matrix
	const Ref<VectorXf>& y,     // phenotype vector
	const int subset_size,      // 
	Ref<VectorXf> snp_blup
)
{
	// prep scaling for error variance weight
	VectorXf r_rsqrt = rvec.array().rsqrt();
	
	std::cout<<"\nStarted weighted least squares.\n";
	// extract needed dimensions
	const int n = y.size();
	const int m = Z.size();
	
	// qr decomposition of covariate matrix
	MatrixXf Q = r_rsqrt.asDiagonal() * X;
	ColPivHouseholderQR<MatrixXf> qr(Q);
	Q = qr.householderQ() * MatrixXf::Identity(n, qr.rank());
	
	// apply projection to phenotype vector
	VectorXf y_proj = r_rsqrt.cwiseProduct(y);
	y_proj = y_proj - Q * (Q.transpose() * y_proj);
	std::cout<<"Projected covariates out of phenotypes.\n";
	
	// fit suitable subsets into memory
	int num_subsets = m / subset_size;
	int num_left_cols = m % subset_size;
	MatrixXf mat1 = MatrixXf::Zero(n, subset_size);
	MatrixXf mat2 = MatrixXf::Zero(n, subset_size);
	
	std::cout<<"Dividing LHS into "<<num_subsets+(num_left_cols>0)<<"x"<<num_subsets+(num_left_cols>0)<<" blocks.\n";
	
	// genotype mean and weight
	VectorXf gm(m); gm.setZero();
	#pragma omp parallel for
	for(int i=0; i<m; ++i) {
		int nz_bits = 0;
		for(int j=0; j<n*2; ++j) {
			if(Z[i][j]) nz_bits ++;
		}
		gm(i) = nz_bits/float(n);
	}
	
	MatrixXf XtX(m, m);
	VectorXf Xty(m);
	
	for(int i=0; i<num_subsets; ++i) {
		#pragma omp parallel for
		for(int j=0; j<subset_size; ++j) {
			float mg = gm(i*subset_size+j);
			for(int k=0, kk=0; kk<n; k+=2, ++kk) {
				mat1(kk,j) = (Z[i*subset_size+j][k] + Z[i*subset_size+j][k+1]) - mg;
			}
		}
		mat1.array().colwise() *= r_rsqrt.array();
		mat1 = mat1 - Q * (Q.transpose() * mat1);
		
		XtX.block(i*subset_size,i*subset_size,subset_size,subset_size) = mat1.transpose() * mat1;
		Xty.segment(i*subset_size, subset_size) = mat1.transpose() * y_proj;
		
		std::cout<<"Computed block "<<i<<","<<i<<"\n";
		
		mat2.resize(n, subset_size);
		for(int ip=i+1; ip<num_subsets; ++ip) {
			#pragma omp parallel for
			for(int j=0; j<subset_size; ++j) {
				float mg = gm(ip*subset_size+j);
				for(int k=0, kk=0; kk<n; k+=2, ++kk) {
					mat2(kk,j) = (Z[ip*subset_size+j][k] + Z[ip*subset_size+j][k+1]) - mg;
				}
			}
			mat2.array().colwise() *= r_rsqrt.array();
			XtX.block(ip*subset_size,i*subset_size,subset_size,subset_size) = mat2.transpose() * mat1;
			
			std::cout<<"Computed block "<<ip<<","<<i<<"\n";
		}
		
		if(num_left_cols != 0) {
			mat2.resize(n, num_left_cols);
			#pragma omp parallel for
			for(int j=0; j<num_left_cols; ++j) {
				float mg = gm(num_subsets*subset_size+j);
				for(int k=0, kk=0; kk<n; k+=2, ++kk) {
					mat2(kk,j) = (Z[num_subsets*subset_size+j][k] + Z[num_subsets*subset_size+j][k+1]) - mg;
				}
			}
			mat2.array().colwise() *= r_rsqrt.array();
			XtX.block(num_subsets*subset_size,i*subset_size,num_left_cols,subset_size) = mat2.transpose() * mat1;
			
			std::cout<<"Computed block "<<num_subsets<<","<<i<<"\n";
		}
	}
	
	if(num_left_cols != 0) {
		mat1.resize(n, num_left_cols);
		#pragma omp parallel for
		for(int j=0; j<num_left_cols; ++j) {
			float mg = gm(num_subsets*subset_size+j);
			for(int k=0, kk=0; kk<n; k+=2, ++kk) {
				mat1(kk,j) = (Z[num_subsets*subset_size+j][k] + Z[num_subsets*subset_size+j][k+1]) - mg;
			}
		}
		mat1.array().colwise() *= r_rsqrt.array();
		mat1 = mat1 - Q * (Q.transpose() * mat1);
		
		XtX.bottomRightCorner(num_left_cols,num_left_cols) = mat1.transpose() * mat1;
		Xty.tail(num_left_cols) = mat1.transpose() * y_proj;
		
		std::cout<<"Computed block "<<num_subsets<<","<<num_subsets<<"\n";
	}
	
	// free memory
	mat1.resize(0,0);
	mat2.resize(0,0);
	
	XtX.diagonal().array() += 1e-4*XtX.diagonal().mean();
	std::cout<<"Computing Cholesky decomposition...\n";
	LLT<Ref<MatrixXf> > llt(XtX);
	snp_blup = llt.solve(Xty);
	
	std::cout<<"\nCompleted weighted least squares."<<std::endl;
}

