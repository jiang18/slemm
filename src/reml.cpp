#include "reml.h"
#include <bits/stdc++.h>
#include <tuple>


// linear operator functions

void H0_quadform_mv(
	Ref<MatrixXf> mat1,     // two temporary matrices 
	Ref<MatrixXf> mat2,     // to speed up char-float conversion
	const Ref<VectorXf>& rr,// reciprocal square root of R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& m, // genotype mean
	const Ref<VectorXf>& w, // weight for grm
	const float tau0,       // (1-h2)/h2
	const Ref<MatrixXf>& Q, // qr decomposition of covariate matrix
	const Ref<MatrixXf>& v, // input matrix
	Ref<MatrixXf> loo       // linear operator output
)
{
	const int subset_size = mat1.cols();
	const int num_left_cols = mat2.cols();
	const int num_subsets = Z.size()/subset_size;
	const int nn = rr.size();
	
	MatrixXf v_adj = rr.asDiagonal() * (v - Q * (Q.transpose() * v));
	MatrixXf ZPv = MatrixXf::Zero(v.rows(), v.cols());
	for(int i=0; i<num_subsets; ++i) {
		#pragma omp parallel for
		for(int j=0; j<subset_size; ++j) {
			float mg = m(i*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat1(kk,j) = (Z[i*subset_size+j][k] + Z[i*subset_size+j][k+1]) - mg;
			}
		}
		ZPv.noalias() +=  mat1 * (w.segment(i*subset_size,subset_size).asDiagonal() * (mat1.transpose() * v_adj));
	}
	if(num_left_cols != 0) {
		#pragma omp parallel for
		for(int j=0; j<num_left_cols; ++j) {
			float mg = m(num_subsets*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat2(kk,j) = (Z[num_subsets*subset_size+j][k] + Z[num_subsets*subset_size+j][k+1]) - mg;
			}
		}
		ZPv.noalias() +=  mat2 * (w.tail(num_left_cols).asDiagonal() * (mat2.transpose() * v_adj));
	}
	ZPv.array().colwise() *= rr.array();
	loo.noalias() = ZPv + tau0*v - Q * (Q.transpose() * ZPv);
}

void H0_ldet_mv(
	Ref<MatrixXf> mat1,     // two temporary matrices 
	Ref<MatrixXf> mat2,     // to speed up char-float conversion
	const Ref<VectorXf>& rr,// reciprocal square root of R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& m, // genotype mean
	const Ref<VectorXf>& w, // weight for grm
	const float tau0,       // (1-h2)/h2
	const Ref<MatrixXf>& v, // input matrix
	Ref<MatrixXf> loo       // linear operator output
)
{
	const int subset_size = mat1.cols();
	const int num_left_cols = mat2.cols();
	const int num_subsets = Z.size()/subset_size;
	const int nn = rr.size();

	loo.setZero();
	MatrixXf rv = rr.asDiagonal() * v;
	for(int i=0; i<num_subsets; ++i) {
		#pragma omp parallel for
		for(int j=0; j<subset_size; ++j) {
			float mg = m(i*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat1(kk,j) = (Z[i*subset_size+j][k] + Z[i*subset_size+j][k+1]) - mg;
			}
		}
		loo.noalias() +=  mat1 * (w.segment(i*subset_size,subset_size).asDiagonal() * (mat1.transpose() * rv));
	}
	if(num_left_cols != 0) {
		#pragma omp parallel for
		for(int j=0; j<num_left_cols; ++j) {
			float mg = m(num_subsets*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat2(kk,j) = (Z[num_subsets*subset_size+j][k] + Z[num_subsets*subset_size+j][k+1]) - mg;
			}
		}
		loo.noalias() +=  mat2 * (w.tail(num_left_cols).asDiagonal() * (mat2.transpose() * rv));
	}
	loo.array().colwise() *= rr.array();
	
	loo.noalias() += tau0 * v;
}

void lnr_op_select(
	Ref<MatrixXf> tmp_mat1, // two temporary matrices 
	Ref<MatrixXf> tmp_mat2, // to speed up char-float conversion
	const std::string lnr,  // linear operator name
	const Ref<VectorXf>& rr,// reciprocal square root of R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& m, // genotype mean
	const Ref<VectorXf>& w, // weight for grm
	const float tau0,       // (1-h2)/h2
	const Ref<MatrixXf>& Q, // qr decomposition of covariate matrix
	const Ref<MatrixXf>& v, // input matrix
	Ref<MatrixXf> loo       // linear operator output
)
{
	if(lnr == "quadform") {
		H0_quadform_mv(tmp_mat1, tmp_mat2, rr, Z, m, w, tau0, Q, v, loo);
	} else if(lnr == "ldet") {
		H0_ldet_mv(tmp_mat1, tmp_mat2, rr, Z, m, w, tau0, v, loo);
	} else {
		std::cout<<"Unknown linear operator in lnr_op_select()"<<std::endl;
		exit(1);
	}
}

// return W(Z-m)'v
MatrixXf weighted_rmatvec(
	Ref<MatrixXf> mat1,     // two temporary matrices 
	Ref<MatrixXf> mat2,     // to speed up char-float conversion
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& m, // genotype mean
	const Ref<VectorXf>& w, // weight for grm
	const Ref<MatrixXf>& v  // input matrix
)
{
	const int subset_size = mat1.cols();
	const int num_left_cols = mat2.cols();
	const int num_subsets = Z.size()/subset_size;
	const int nn = v.rows();

	MatrixXf loo(m.size(),v.cols());
	for(int i=0; i<num_subsets; ++i) {
		#pragma omp parallel for
		for(int j=0; j<subset_size; ++j) {
			float mg = m(i*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat1(kk,j) = (Z[i*subset_size+j][k] + Z[i*subset_size+j][k+1]) - mg;
			}
		}
		loo.middleRows(i*subset_size,subset_size).noalias() =  w.segment(i*subset_size,subset_size).asDiagonal() * (mat1.transpose() * v);
	}
	if(num_left_cols != 0) {
		#pragma omp parallel for
		for(int j=0; j<num_left_cols; ++j) {
			float mg = m(num_subsets*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat2(kk,j) = (Z[num_subsets*subset_size+j][k] + Z[num_subsets*subset_size+j][k+1]) - mg;
			}
		}
		loo.bottomRows(num_left_cols).noalias() =  w.tail(num_left_cols).asDiagonal() * (mat2.transpose() * v);
	}
	return loo;
}

// return (Z-m)'v
MatrixXf rmatvec(
	Ref<MatrixXf> mat1,     // two temporary matrices 
	Ref<MatrixXf> mat2,     // to speed up char-float conversion
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& m, // genotype mean
	const Ref<MatrixXf>& v  // input matrix
)
{
	const int subset_size = mat1.cols();
	const int num_left_cols = mat2.cols();
	const int num_subsets = Z.size()/subset_size;
	const int nn = v.rows();

	MatrixXf loo(m.size(),v.cols());
	for(int i=0; i<num_subsets; ++i) {
		#pragma omp parallel for
		for(int j=0; j<subset_size; ++j) {
			float mg = m(i*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat1(kk,j) = (Z[i*subset_size+j][k] + Z[i*subset_size+j][k+1]) - mg;
			}
		}
		loo.middleRows(i*subset_size,subset_size).noalias() = mat1.transpose() * v;
	}
	if(num_left_cols != 0) {
		#pragma omp parallel for
		for(int j=0; j<num_left_cols; ++j) {
			float mg = m(num_subsets*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat2(kk,j) = (Z[num_subsets*subset_size+j][k] + Z[num_subsets*subset_size+j][k+1]) - mg;
			}
		}
		loo.bottomRows(num_left_cols).noalias() =  mat2.transpose() * v;
	}
	return loo;
}

// return (Z-m)v
MatrixXf matvec(
	Ref<MatrixXf> mat1,     // two temporary matrices 
	Ref<MatrixXf> mat2,     // to speed up char-float conversion
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& m, // genotype mean
	const Ref<MatrixXf>& v  // input matrix
)
{
	const int subset_size = mat1.cols();
	const int num_left_cols = mat2.cols();
	const int num_subsets = Z.size()/subset_size;
	const int nn = Z[0].size()/2;

	MatrixXf loo = MatrixXf::Zero(nn,v.cols());
	for(int i=0; i<num_subsets; ++i) {
		#pragma omp parallel for
		for(int j=0; j<subset_size; ++j) {
			float mg = m(i*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat1(kk,j) = (Z[i*subset_size+j][k] + Z[i*subset_size+j][k+1]) - mg;
			}
		}
		loo.noalias() += mat1 * v.middleRows(i*subset_size,subset_size);
	}
	if(num_left_cols != 0) {
		#pragma omp parallel for
		for(int j=0; j<num_left_cols; ++j) {
			float mg = m(num_subsets*subset_size+j);
			for(int k=0, kk=0; kk<nn; k+=2, ++kk) {
				mat2(kk,j) = (Z[num_subsets*subset_size+j][k] + Z[num_subsets*subset_size+j][k+1]) - mg;
			}
		}
		loo.noalias() +=  mat2 * v.bottomRows(num_left_cols);
	}
	return loo;
}

//######################## L_Seed ##########################
//## constructs bases for Krylov subspaces:               ##
//## B, (A+sI)B, (A+sI)²B, ...                            ##
//##########################################################

typedef struct {
	std::vector<MatrixXf> U;
	MatrixXf beta;
	MatrixXf delta;
} lseed;

lseed L_Seed (
	Ref<MatrixXf> tmp_mat1, // two temporary matrices 
	Ref<MatrixXf> tmp_mat2, // to speed up char-float conversion
// start of linear operator for A
	const std::string lnr,  // linear operator name
	const Ref<VectorXf>& rr,// reciprocal square root of R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& m, // genotype mean
	const Ref<VectorXf>& w, // weight for grm
	const float tau0,       // (1-h2)/h2
	const Ref<MatrixXf>& Q, // qr decomposition of covariate matrix
// end of linear operator for A
	const Ref<MatrixXf>& B, // RHS
	const float tol = 5e-4, // relative tolerance to B
	const int p_freq = 5,   // print frequency
	const int maxit = 600,  // maximum iteration count
	bool verbose = true     // print extra information
)
{
	lseed out;
	
	int n = Q.rows();
	int t = B.cols();
	RowVectorXf B_norm = B.colwise().norm();
	
	MatrixXf R = MatrixXf::Zero(n,t); // CG residuals
	out.U.clear();                    // orthonormal bases for K-subspaces
	out.U.push_back(MatrixXf::Zero(n,t));
	
	// coefficients
	MatrixXf rho = MatrixXf::Zero(t,maxit);
	MatrixXf beta = MatrixXf::Zero(t,maxit);
	MatrixXf omega = MatrixXf::Zero(t,maxit);
	MatrixXf gamma = MatrixXf::Ones(t,maxit);
	MatrixXf delta = MatrixXf::Zero(t,maxit);
	
	// initial values
	rho.col(0) = B_norm;
	beta.col(0) = rho.col(0);
	R = B;
	out.U[0] = R.array().rowwise() / beta.col(0).array().transpose();
	
	int j = 0;
	bool cnvg = false;
	
	MatrixXf tmp(n, t);
	
	while(j < maxit-1) {
		out.U.push_back(MatrixXf::Zero(n,t));
		// Lanczos iteration
		lnr_op_select(tmp_mat1, tmp_mat2, lnr, rr, Z, m, w, tau0, Q, out.U[j], tmp);
		if(j >0) tmp.array() -= out.U[j-1].array().rowwise() * beta.col(j).array().transpose();
		delta.col(j) = (out.U[j].transpose() * tmp).diagonal();
		tmp.array() -= out.U[j].array().rowwise() * delta.col(j).array().transpose();
		beta.col(j+1) = tmp.colwise().norm();
		out.U[j+1] = tmp.array().rowwise() / beta.col(j+1).array().transpose();
		
		// CG coefficents update
		if(j == 0) gamma.col(j) = delta.col(j).array().inverse();
		else gamma.col(j) = (delta.col(j).array() - omega.col(j-1).array()/gamma.col(j-1).array()).inverse();
		omega.col(j) = (beta.col(j+1).array() * gamma.col(j).array()).square();
		rho.col(j+1) = -beta.col(j+1).array() * gamma.col(j).array() * rho.col(j).array();
		
		// CG vectors update
		R = out.U[j+1].array().rowwise() * rho.col(j+1).array().transpose();
		
		++j;
		
		float res_norm = R.colwise().norm().cwiseQuotient(B_norm).maxCoeff();
		
		if(verbose and j % p_freq ==0)
			std::cout<<"Relative error at step "<<j<<" is "<<res_norm<<std::endl;
		if(res_norm <= tol) {
			cnvg = true;
			break;
		}
	}
	if(verbose and cnvg) std::cout<<"Converged after "<<j<<" iterations."<<std::endl;
	else if(verbose) std::cout<<"Failed to converge after "<<maxit<<" iterations."<<std::endl;
	
	R.resize(0,0);
	tmp.resize(0,0);
	omega.resize(0,0);
	gamma.resize(0,0);
	rho.resize(0,0);
	
	// Jacobi coefficients for Lanczos polys
	out.beta = beta.leftCols(j+1);
	beta.resize(0,0);
	out.delta = delta.leftCols(j);
	delta.resize(0,0);
	
	return out;
}

//######################## L_Solve #########################
//## solves (A + sigma I)X = B using results from L_Seed) ##
//##########################################################
void solve_symtridiag_inplace(
	const Ref<VectorXf>& a,    // LHS subdiagonal
	const Ref<VectorXf>& b,    // LHS diagonal
	Ref<VectorXf> d            // RHS and sol
)
{
	int n = b.size();
	VectorXf c = a.segment(1, n-1);
	
	n--;
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i]*c[i-1];
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
	}

	d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

	for (int i = n; i-- > 0;) {
		d[i] -= c[i]*d[i+1];
	}
}

MatrixXf L_Solve(
	const lseed& seedSystem, // ouput from L_Seed
	const Ref<MatrixXf>& B,  // RHS
	const float sigma,       // shift >=0
	const float tol = 1e-6,  // relative tolerance for CG (fixed, sufficiently small)
	bool verbose= false,     // verbose output
	const int p_freq = 5     // print frequency (if verbose)
)
{
	RowVectorXf B_norm = B.colwise().norm();
	// copy delta only
	MatrixXf delta = seedSystem.delta.array() + sigma;
	
	int n = seedSystem.U[0].rows();
	int t = seedSystem.U[0].cols();
	int maxit = seedSystem.U.size();
	
	MatrixXf X(n,t);     // approximate soln
	MatrixXf P(n,t);     // search directions
	MatrixXf R(n,t);     // residuals
	
	// coefficients
	MatrixXf omega = MatrixXf::Zero(t,maxit);
	MatrixXf gamma = MatrixXf::Ones(t,maxit);
	MatrixXf rho = MatrixXf::Zero(t,maxit);
	
	// initial values
	rho.col(0) = B_norm;
	R = B;
	P = B;
	X.setZero();
	
	int j=0;
	bool cnvg = false;
	
	while(j<maxit-1) {
		// CG coefficents update
		if(j == 0) gamma.col(j) = delta.col(j).array().inverse();
		else gamma.col(j) = (delta.col(j).array() - omega.col(j-1).array()/gamma.col(j-1).array()).inverse();
		omega.col(j) = (seedSystem.beta.col(j+1).array() * gamma.col(j).array()).square();
		rho.col(j+1) = -seedSystem.beta.col(j+1).array() * gamma.col(j).array() * rho.col(j).array();
		
		// CG vectors update
		X.array() += P.array().rowwise() * gamma.col(j).array().transpose();
		R = seedSystem.U[j+1].array().rowwise() * rho.col(j+1).array().transpose();
		P = R.array() + P.array().rowwise() * omega.col(j).array().transpose();
		
		++j;
		
		float res_norm = R.colwise().norm().cwiseQuotient(B_norm).maxCoeff();
		
		if(verbose and j % p_freq ==0)
			std::cout<<"Relative error at step "<<j<<" is "<<res_norm<<std::endl;
		if(res_norm <= tol) {
			cnvg = true;
			break;
		}
	}
	if(verbose and cnvg) std::cout<<"Converged after "<<j<<" iterations."<<std::endl;
	else if(verbose) std::cout<<"Failed to converge after "<<maxit<<" iterations."<<std::endl;
	
	return X;
}

//######################## SLQ_LDet #################################
//## returns approximate log det (A + sI) given spectral           ##
//## decompositions of Jacobi matrices from Lanczos decompositions ##
//## of seed Krylov subspaces for probes                           ##
//###################################################################
// log det A = tr log A

float SLQ_LDet(
	const Ref<MatrixXf>& D,     // eigenvalues of Jacobi matrices
	const Ref<MatrixXf>& W,     // eigenvectors of Jacobi matrices
	const int n,                // dimension of A
	const int n_V,              // number of probing vectors
	float sigma = 0             // shift
)
{
	float soln = ( W.array().square() * (D.array() + sigma).log() ).sum();
	
	return(soln*n/n_V);         // estimate of tr log (A + sI)
}

//############################# SLDF_REML ##############################
//## zero order REML estimation using (shifted) Lanczos conjugate     ##
//## gradients and (shifted) stochastic Lanczos quadrature            ##
//######################################################################
float REML_criterion(
	const float ss,
	const float tau0,
	const Ref<VectorXf>& y_proj,
	const lseed& seed_y,
	const Ref<MatrixXf>& X,
	const lseed& seed_X,
	const Ref<MatrixXf>& D_V,
	const Ref<MatrixXf>& W_V,
	const float ySy,
	float& v_g,
	float& v_e,
	float& lrt
)
{
	const int n = X.rows();
	const int c = X.cols();
	
	float tau = (1-ss)/ss;
	float sigma = tau - tau0;
	float yPy = (y_proj.transpose() * L_Solve(seed_y, y_proj, sigma))(0);
	
	v_g = yPy/(n-c);
	v_e = v_g*tau;
	
	float ldet_H = SLQ_LDet(D_V, W_V, n, D_V.rows(), sigma);
	
	lrt = n*log(ySy/yPy) - ldet_H;
	
	float slogLL = (n-c)*log(ySy/yPy) - (X.transpose() * L_Solve(seed_X, X, sigma)).colPivHouseholderQr().logAbsDeterminant() - ldet_H;
	// add log det(X'X) to calculate LRT statistic by REML
	return(slogLL);
}

// golden section search
// note that we always have 0 < a < b.
std::tuple<float,float,float> REML_GSS (
	float a,
	float b,
	const float tau0,
	const Ref<VectorXf>& y_proj,
	const lseed& seed_y,
	const Ref<MatrixXf>& X,
	const lseed& seed_X,
	const Ref<MatrixXf>& D_V,
	const Ref<MatrixXf>& W_V,
	const float tol )
{
	const double invphi = 0.618033988749894848204586834365638117720309179805762862135;
	const double invphi2 = 0.381966011250105151795413165634361882279690820194237137864;
	
	float ve, vg, lrt, fm;
	float ySy = y_proj.squaredNorm();
	
	std::cout<<"\nDoing golden section search in ("<<a<<", "<<b<<")"<<std::endl;
	float h = b - a;
	int nitr;
	if(h < tol) nitr = 1;
	else nitr = ceil(log(tol/h)/log(invphi));
	
	std::cout<<nitr<<" iterations:"<<std::endl;
	
	float c, d, fc, fd;
	c = a + invphi2 * h;
	d = a + invphi * h;
	fc = REML_criterion(c, tau0, y_proj, seed_y, X, seed_X, D_V, W_V, ySy, vg, ve, lrt);
	fd = REML_criterion(d, tau0, y_proj, seed_y, X, seed_X, D_V, W_V, ySy, vg, ve, lrt);
	
	for(int i=0; i<nitr; ++i) {
		if(fc > fd) {
			b = d;
			d = c;
			fd = fc;
			h *= invphi;
			c = a + invphi2 * h;
			fc = REML_criterion(c, tau0, y_proj, seed_y, X, seed_X, D_V, W_V, ySy, vg, ve, lrt);
			fm = fc;
		} else {
			a = c;
			c = d;
			fc = fd;
			h *= invphi;
			d = a + invphi * h;
			fd = REML_criterion(d, tau0, y_proj, seed_y, X, seed_X, D_V, W_V, ySy, vg, ve, lrt);
			fm = fd;
		}
		std::cout<<"Itr "<<i+1<<": h2 in ("<<a<<", "<<b<<") and 2logLL = "<<fm<<std::endl;
	}
	std::cout<<"\nCompleted golden section search."<<std::endl;
	return std::make_tuple(vg, ve, lrt);
}

void SLDF_REML(
	const Ref<VectorXf>& rvec,  // error variance weight R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<MatrixXf>& X,     // covariate matrix
	const Ref<VectorXf>& y,     // phenotype vector
	const Ref<VectorXf>& persnp,// per-SNP heritability weight
	const int subset_size,
	const int rng_seed,         // random seed
	const int n_V,              // number of random probes
	const float tol_L,          // relative lanczos tolerance
	Ref<VectorXf> Py,
	Ref<VectorXf> snp_blup,
	Ref<VectorXf> blue,
	float& vg,
	float& ve,
	float& lrt,
	const float s2max,          // maximal h2 value
	const float s2min,          // minimal h2 value
	const float tol_VC,         // abs. tolerance for h2 estimation
	const bool verbose,         // verbose output
	const int p_freq            // print frequency
)
{
	// prep scaling for error variance weight
	VectorXf r_rsqrt = rvec.array().rsqrt();
	
	std::cout<<"\nStarted stochastic Lanczos DF-REML."<<std::endl;
	// extract needed dimensions
	const int n = y.size();
	const int m = Z.size();
	
	// qr decomposition of covariate matrix
	// use Q as X
	MatrixXf Q = r_rsqrt.asDiagonal() * X;
	ColPivHouseholderQR<MatrixXf> qr(Q);
	Q = qr.householderQ() * MatrixXf::Identity(n, qr.rank());
	MatrixXf qrR = qr.matrixR().topLeftCorner(qr.rank(), qr.rank()).template triangularView<Upper>();
	if(Q.cols() != X.cols()) throw("\nFixed-effects design matrix is not of full rank.\n");
	
	// apply projection to phenotype vector
	VectorXf y_proj = r_rsqrt.cwiseProduct(y);
	y_proj = y_proj - Q * (Q.transpose() * y_proj);
	std::cout<<"Projected covariates out of phenotypes."<<std::endl;

	// sample normalized Rademacher probing vetors
	std::cout<<"Generating "<<n_V<<" Rademacher probing vetors."<<std::endl;
	MatrixXf V = MatrixXf::Ones(n, n_V);
	std::mt19937 rng(rng_seed+54321);
	std::uniform_real_distribution<float> runif(0.0, 1.0);
	for(int j=0; j<n_V; ++j) {
		for(int i=0; i<n; ++i)
			if(runif(rng) < 0.5) V(i,j) = -1.0;
	}
	RowVectorXf V_norm = V.colwise().norm();
	V = V.array().rowwise() / V_norm.array();
	
	// fit suitable subsets into memory
	int num_left_cols = m % subset_size;
	MatrixXf tmp_mat1 = MatrixXf::Zero(n, subset_size);
	MatrixXf tmp_mat2 = MatrixXf::Zero(n, num_left_cols);
	
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
	VectorXf gw = (gm.array() * (2.0-gm.array())/2.).inverse()/float(m) * ( persnp.array()/persnp.mean() );
	
	// Lanczos decompositions of seed systems
	float tau0 = (1-s2max)/s2max;
	std::cout<<"\nPerforming Lanczos decomposition for phenotypes:"<<std::endl;
	lseed seed_y = L_Seed(tmp_mat1, tmp_mat2, "quadform", r_rsqrt, Z, gm, gw, tau0, Q, y_proj, tol_L, p_freq);
	std::cout<<"\nPerforming Lanczos decomposition for covariates:"<<std::endl;
	lseed seed_Q = L_Seed(tmp_mat1, tmp_mat2, "ldet", r_rsqrt, Z, gm, gw, tau0, Q, Q, tol_L, p_freq);
	std::cout<<"\nPerforming Lanczos decomposition for Rademacher vectors:"<<std::endl;
	lseed seed_V = L_Seed(tmp_mat1, tmp_mat2, "ldet", r_rsqrt, Z, gm, gw, tau0, Q, V, tol_L, p_freq);
	
	// free memory
	V.resize(0,0);
	seed_V.U.clear();
	
	// decompose Jacobi matrices for SLQ
	int nld = seed_V.delta.cols();
	MatrixXf W_V = MatrixXf::Zero(n_V, nld);  // eigenvalues
	MatrixXf D_V = MatrixXf::Zero(n_V, nld);  // first elements of eigenvectors
	SelfAdjointEigenSolver<MatrixXf> eigs(nld);
	for(int i=0; i<n_V; ++i) {
		eigs.computeFromTridiagonal(seed_V.delta.row(i), seed_V.beta.row(i).segment(1,nld-1));
		D_V.row(i) = eigs.eigenvalues();
		W_V.row(i) = eigs.eigenvectors().row(0);
	}
	
	// free memory
	seed_V.beta.resize(0,0);
	seed_V.delta.resize(0,0);
	
	std::tie(vg, ve, lrt) = REML_GSS(s2min, s2max, tau0, y_proj, seed_y, Q, seed_Q, D_V, W_V, tol_VC );
	
	// free memory
	W_V.resize(0,0);
	D_V.resize(0,0);
	seed_Q.U.clear();
	seed_Q.beta.resize(0,0);
	seed_Q.delta.resize(0,0);
	
	float tau1 = ve/vg;
	Py = L_Solve(seed_y, y_proj, tau1-tau0).col(0);
	Py = Py - Q*(Q.transpose() * Py);
	Py = r_rsqrt.cwiseProduct(Py);
	snp_blup = weighted_rmatvec(tmp_mat1, tmp_mat2, Z, gm, gw, Py).col(0);
	std::cout<<"\nVg,"<<vg<<std::endl;
	std::cout<<"Ve,"<<ve<<std::endl;
	std::cout<<"h2,"<<vg/(vg+ve)<<std::endl;
	std::cout<<"LRT,"<<lrt<<std::endl;
	
	// compute blue of fixed effects
	y_proj = matvec(tmp_mat1, tmp_mat2, Z, gm, snp_blup);
	y_proj += tau1 * rvec.cwiseProduct(Py);
	y_proj = y - y_proj;
	blue = (X.transpose() * X).llt().solve(X.transpose() * y_proj);
	
	// free memory
	y_proj.resize(0);
	tmp_mat1.resize(0,0);
	tmp_mat2.resize(0,0);
	
	std::cout<<"\nCompleted stochastic Lanczos DF-REML."<<std::endl;
	
	if(s2max - vg/(vg+ve) < tol_VC) {
		std::cout<<"\nWarning: the h2 estimate is virtually the same as the max you specified.\n";
		std::cout<<"Rerun with a higher --max_heritability."<<std::endl;
	}
}

//############################ L_FOMC_REML #################################
//## extension of BOLT-LMM algorithm to recycle Krylov subspace bases     ##
//## involved in solving linear systems                                   ##
//##########################################################################
float REML_first_order(
	Ref<MatrixXf> tmp_mat1, // two temporary matrices 
	Ref<MatrixXf> tmp_mat2, // to speed up char-float conversion
	const float ss,
	const float tau0,
	const Ref<VectorXf>& rr,// reciprocal square root of R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<VectorXf>& gm,
	const Ref<VectorXf>& gw,
	const Ref<MatrixXf>& Q,
	const Ref<MatrixXf>& V,
	const lseed& seed_V,
	Ref<MatrixXf> Py       // returned value is exactly Py with V = ZWZ'+R(varE/varA).
)
{
	const int n_MC = Py.cols()-1;
	
	float tau = (1-ss)/ss;
	float sigma = tau - tau0;
	MatrixXf PV = L_Solve(seed_V, V, sigma);
	Py.rightCols(n_MC) = PV.middleCols(1,n_MC) + sqrt(tau) * PV.rightCols(n_MC);
	Py.col(0) = PV.col(0);
	Py = Py - Q * (Q.transpose() * Py);
	
	// scaled BLUP
	MatrixXf e_BLUP = Py;
	Py.array().colwise() *= rr.array();
	MatrixXf u_BLUP = rmatvec(tmp_mat1, tmp_mat2, Z, gm, Py);
	u_BLUP = u_BLUP.array().colwise() * gw.array().sqrt();
	
	float f = log((u_BLUP.leftCols(1).squaredNorm()/e_BLUP.leftCols(1).squaredNorm()) / (u_BLUP.rightCols(n_MC).squaredNorm()/e_BLUP.rightCols(n_MC).squaredNorm()));
	return(f);
}

void LMC_REML(
	const Ref<VectorXf>& rvec,  // error variance weight R
	const std::vector<std::vector<bool> >& Z, // genotype bits
	const Ref<MatrixXf>& X,     // covariate matrix
	const Ref<VectorXf>& y,     // phenotype vector
	const Ref<VectorXf>& persnp,// per-SNP heritability weight
	const int subset_size,
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
	const float s2max,          // maximal h2 value
	const float s2min,          // minimal h2 value
	const bool verbose,         // verbose output
	const int p_freq            // print frequency
)
{
	/// prep scaling for error variance weight
	VectorXf r_rsqrt = rvec.array().rsqrt();
	
	std::cout<<"\nStarted Lanczos MC-REML."<<std::endl;
	// extract needed dimensions
	const int n = y.size();
	const int m = Z.size();
	
	// qr decomposition of covariate matrix
	MatrixXf Q = r_rsqrt.asDiagonal() * X;
	ColPivHouseholderQR<MatrixXf> qr(Q);
	const int c = qr.rank();
	Q = qr.householderQ() * MatrixXf::Identity(n, c);
	if(Q.cols() != X.cols()) throw("\nFixed-effects design matrix is not of full rank.\n");
	
	// genotype means and weights
	VectorXf gm(m); gm.setZero();
	#pragma omp parallel for
	for(int i=0; i<m; ++i) {
		int nz_bits = 0;
		for(int j=0; j<n*2; ++j) {
			if(Z[i][j]) nz_bits ++;
		}
		gm(i) = nz_bits/float(n);
	}
	VectorXf gw = (gm.array() * (2.0-gm.array())/2.).inverse()/float(m) * ( persnp.array()/persnp.mean() );
	
	// fit suitable subsets into memory
	int num_left_cols = m % subset_size;
	MatrixXf tmp_mat1 = MatrixXf::Zero(n, subset_size);
	MatrixXf tmp_mat2 = MatrixXf::Zero(n, num_left_cols);
	
	// 1+2*n_MC columns for Lanczos decomposition
	MatrixXf V = MatrixXf(n, 1+2*n_MC);
	// y is the first column
	V.col(0) = y;
	// MC samples for genetic effects
	std::cout<<"Generating "<<n_MC<<" MC samples for genetic effects."<<std::endl;
	MatrixXf u_MC = MatrixXf(m, n_MC);
	std::mt19937 rng(rng_seed+54321);
	std::normal_distribution<float> rnorm(0.0, 1.0);
	for(int j=0; j<n_MC; ++j) {
		for(int i=0; i<m; ++i)
			u_MC(i,j) = rnorm(rng);
	}
	u_MC = u_MC.array().colwise() * gw.array().sqrt();
	V.middleCols(1,n_MC) = matvec(tmp_mat1, tmp_mat2, Z, gm, u_MC);
	u_MC.resize(0,0);
	// MC samples for residuals
	std::cout<<"Generating "<<n_MC<<" MC samples for residuals."<<std::endl;
	for(int j=1+n_MC; j<1+2*n_MC; ++j) {
		for(int i=0; i<n; ++i)
			V(i,j) = rnorm(rng);
	}
	// project covariates out
	std::cout<<"Projecting covariates out of phenotypes and MC samples."<<std::endl;
	V.leftCols(1+n_MC).array().colwise() *= r_rsqrt.array();
	V = V - Q * (Q.transpose() * V);
	
	// Lanczos decompositions of seed systems
	std::cout<<"\nPerforming Lanczos decomposition:"<<std::endl;
	float tau0 = (1-s2max)/s2max;
	// follow BOLT to relax the CG tolerance 
	lseed seed_V = L_Seed(tmp_mat1, tmp_mat2, "quadform", r_rsqrt, Z, gm, gw, tau0, Q, V, 10*tol_L, p_freq);
	
	// secant iteration
	std::cout<<"\nSecant iterations:"<<std::endl;
	MatrixXf PyMC(n, 1+n_MC);
	float s2 = s2max/2;
	float log_tau1 = log((1-s2)/s2);
	float f1 = REML_first_order(tmp_mat1, tmp_mat2, s2, tau0, r_rsqrt, Z, gm, gw, Q, V, seed_V, PyMC);
	printf("  Iter 1: h2=%.4f log((1-h2)/h2)=%.4f funct=%.4f\n", s2, log_tau1, f1);
	fflush(stdout);
	
	if(f1 < 0) s2 = s2max/4;
	else s2 = s2max;
	float log_tau2 = log((1-s2)/s2);
	float f2 = REML_first_order(tmp_mat1, tmp_mat2, s2, tau0, r_rsqrt, Z, gm, gw, Q, V, seed_V, PyMC);
	printf("  Iter 2: h2=%.4f log((1-h2)/h2)=%.4f funct=%.4f\n", s2, log_tau2, f2);
	fflush(stdout);
	
	for(int si = 3; si <= 9; ++si) {
		float log_tau_si = (log_tau1*f2 - log_tau2*f1) / (f2 - f1);
		if(log_tau_si > 9) log_tau_si = 9;          // constrain h2 to 0.0001-0.9999
		else if(log_tau_si < -9) log_tau_si = -9;   // added on Mar 18, 2022
		
		if(std::abs(log_tau_si - log_tau2) < 0.002) break; // TODO: decide
		
		s2 = 1/(exp(log_tau_si) + 1);
		f1 = f2;
		f2 = REML_first_order(tmp_mat1, tmp_mat2, s2, tau0, r_rsqrt, Z, gm, gw, Q, V, seed_V, PyMC);
		log_tau1 = log_tau2;
		log_tau2 = log_tau_si;
		printf("  Iter %d: h2=%.4f log((1-h2)/h2)=%.4f funct=%.4f\n", si, s2, log_tau_si, f2);
		fflush(stdout);
	}
	// added on Mar 18, 2022
	if(std::abs(f2) > 0.005) std::cout<<"\nWarning: secant iterations not converged. h2 close to 0 or small sample size?\n";
	
	Py = PyMC.col(0);
	vg = y.dot(Py)/float(n-c);
	float tau1 = (1-s2)/s2;
	ve = tau1 * vg;
	
	// free memory
	PyMC.resize(0,0);
	seed_V.U.clear();
	seed_V.beta.resize(0,0);
	seed_V.delta.resize(0,0);
	
	snp_blup = weighted_rmatvec(tmp_mat1, tmp_mat2, Z, gm, gw, Py).col(0);
	std::cout<<"\nVg,"<<vg<<std::endl;
	std::cout<<"Ve,"<<ve<<std::endl;
	std::cout<<"h2,"<<s2<<std::endl;
	
	// compute blue of fixed effects
	V.col(0) = matvec(tmp_mat1, tmp_mat2, Z, gm, snp_blup);
	V.col(0) += tau1 * rvec.cwiseProduct(Py);
	V.col(0) = y - V.col(0);
	blue = (X.transpose() * X).llt().solve(X.transpose() * V.col(0));
	
	std::cout<<"\nCompleted Lanczos MC-REML."<<std::endl;
	
	if(s2 > s2max) {
		std::cout<<"\nWarning: the h2 estimate is larger than the max you specified.\n";
		std::cout<<"Can rerun with a higher --max_heritability to gain accuracy."<<std::endl;
	}
	
	if(do_lmm) {
		std::cout<<"\nStarted approximation of LMM associations."<<std::endl;
		if(fake_geno) {
			std::cout<<"Generating "<<num_qf_markers<<" fake markers (MAF>0.01) to compute quadratic forms."<<std::endl;
			std::normal_distribution<float> rnorm(0.0, 1.0);
			std::uniform_real_distribution<float> runif(0.01,0.99);
			
			std::vector<int> variants;
			V.resize(n, num_qf_markers);
			for(int j=0; j<num_qf_markers; ++j) {
				variants.push_back(-(j+1));
				float a1_freq = runif(rng);
				std::binomial_distribution<int> rbiom(2, a1_freq);
				for(int i=0; i<n; ++i)
					V(i,j) = rbiom(rng);
			}
			
			V = V.rowwise() - V.colwise().mean();
			VectorXf VtPy = V.transpose() * Py;
			V.array().colwise() *= r_rsqrt.array();  // added January 11, 2022 for using reliability in GWA
			V = V - Q * (Q.transpose() * V);
			seed_V = L_Seed(tmp_mat1, tmp_mat2, "quadform", r_rsqrt, Z, gm, gw, tau1, Q, V, tol_L, p_freq);
			MatrixXf VtPV = V.transpose() * L_Solve(seed_V, V, 0.0);
			RowVectorXf vtSv = V.colwise().squaredNorm();
			
			lmm_out.clear();
			gstat elem;
			for(int i=0; i<V.cols(); ++i) {
				elem.vi = variants[i];
				elem.xPy = VtPy(i);
				elem.xPx = VtPV(i,i);
				elem.xBx = 0;
				elem.xSx = vtSv(i);
				lmm_out.push_back(elem);
			}
		} else {
			std::cout<<"Sampling "<<num_qf_markers<<" markers (MAC>20) to compute quadratic forms."<<std::endl;
			std::vector<int> variants;
			for(int i=0; i<m; ++i)
				if(gm[i]*n >20 && (2-gm[i])*n >20) variants.push_back(i);
			std::random_shuffle(variants.begin(), variants.end());
			
			V.resize(n, std::min(num_qf_markers, static_cast<int>(variants.size())));
			#pragma omp parallel for
			for(int i=0; i<V.cols(); ++i) {
				int j = variants[i];
				float mg = gm[j];
				for(int k=0, kk=0; kk<n; k+=2, ++kk) {
					V(kk,i) = (Z[j][k] + Z[j][k+1]) - mg;
				}
			}
			VectorXf VtPy = V.transpose() * Py;
			
			V.array().colwise() *= r_rsqrt.array();  // added January 11, 2022 for using reliability in GWA
			V = V - Q * (Q.transpose() * V);
			seed_V = L_Seed(tmp_mat1, tmp_mat2, "quadform", r_rsqrt, Z, gm, gw, tau1, Q, V, tol_L, p_freq);
			MatrixXf VtPV = V.transpose() * L_Solve(seed_V, V, 0.0);
			RowVectorXf vtSv = V.colwise().squaredNorm();
			
			VectorXf vtBv(V.cols());
			MatrixXf Z_wd(n, window_size);
			MatrixXf ZtZ_wd(window_size, window_size);
			int wds;
			for(int i=0; i<V.cols(); ++i) {
				if(variants[i]<window_size/2) wds = 0;
				else if(variants[i]>m-window_size/2) wds=m-window_size;
				else wds=variants[i]-window_size/2;
				
				#pragma omp parallel for
				for(int j=wds; j<wds+window_size; ++j) {
					float mg = gm[j];
					for(int k=0, kk=0; kk<n; k+=2, ++kk) {
						Z_wd(kk,j-wds) = (Z[j][k] + Z[j][k+1]) - mg;
					}
				}
				Z_wd.array().colwise() *= r_rsqrt.array();  // added January 11, 2022 for using reliability in GWA
				Z_wd = Z_wd - Q * (Q.transpose() * Z_wd);
				ZtZ_wd = Z_wd.transpose() * Z_wd;
				ZtZ_wd.diagonal() += tau1 * gw.segment(wds,window_size).cwiseInverse();
				VectorXf Zv = Z_wd.transpose() * V.col(i);
				
				vtBv(i) = ( vtSv(i) - Zv.dot(ZtZ_wd.llt().solve(Zv)) )/tau1;
			}
			Z_wd.resize(0,0);
			ZtZ_wd.resize(0,0);
			
			lmm_out.clear();
			gstat elem;
			for(int i=0; i<V.cols(); ++i) {
				elem.vi = variants[i];
				elem.xPy = VtPy(i);
				elem.xPx = VtPV(i,i);
				elem.xBx = vtBv(i);
				elem.xSx = vtSv(i);
				lmm_out.push_back(elem);
			}
		}
		std::cout<<"\nCompleted approximation of LMM associations."<<std::endl;
	}
	
	// free memory
	tmp_mat1.resize(0,0);
	tmp_mat2.resize(0,0);
	seed_V.U.clear();
	seed_V.beta.resize(0,0);
	seed_V.delta.resize(0,0);
	V.resize(0,0);
}

