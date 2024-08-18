
#include <getopt.h>
#include <chrono>
#include <ctime>
#include "slemm.h"
#include "reml.h"
#include "wls.h"

void print_help();
void option(int option_num, char **option_str);

int main(int argc, char **argv)
{   
	std::cout<<"*******************************************************************"<<std::endl;
	std::cout<<"* SLEMM by Jicai Jiang"<<std::endl;
	std::cout<<"* Stochastic-Lanczos-Expedited Mixed Models"<<std::endl;
	std::cout<<"* Version 0.90.1 (August 18, 2024)"<<std::endl;
	std::cout<<"* (C) 2021-present, Jicai Jiang, NC State University"<<std::endl;
	std::cout<<"*******************************************************************"<<std::endl;
	
	std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	std::cout<<"Analysis started: "<<std::ctime(&start_time)<<std::endl;
	
	try{ option(argc, argv); }
	catch(const std::string &err_msg){ std::cerr<<"\n"<<err_msg<<std::endl; }
	catch(const char *err_msg){ std::cerr<<"\n"<<err_msg<<std::endl; }
	
	std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout<<"\nAnalysis finished: "<<std::ctime(&end_time);
	std::chrono::duration<double> elapsed_seconds = end - start;
	long time_used = static_cast<long>(elapsed_seconds.count());
	std::cout<<"Computational time: "<<time_used/3600<<":"<<(time_used%3600)/60<<":"<<time_used%60<<std::endl;

	return 0;
}

void option(int option_num, char **option_str) {
	if(option_num < 3) {
		print_help();
		exit(1);
	}
	// variables
	bool pred_flag = false;
	
	bool reml_flag = false;
	bool lmm_flag = false;
	bool wls_flag = false;
	bool lrt_flag = false;
	
	bool iter_weighting = false;
	float sig_chisq = 1e6;
	float indep_r2 = 0.5;
	
	int rng_seed = 0;
	int num_rand_prob = 30;
	bool using_nr_opt = false;
	float max_h2 = 0.7;
	int window_size = 1000;
	bool using_ws_opt = false;
	int num_qf_markers = 30;
	bool fake_geno = false;
	
	std::string marker_estimate_file = "";
	std::string phen_file="", trait_name="", error_variance_weight_name="";
	std::string covar_file="", str_covar_names="";
	std::vector<std::string> covar_names;
	std::string marker_info_file="", marker_group="", marker_effvar_weight_name="";
	std::string binary_genotype_file="";
	std::pair<float, float> beta_pdf_params = std::make_pair(1, 1);
	float min_maf=0;
	int num_threads=1;
	std::string output_file="";
	unsigned subset_size = 1000;
	float rel_tol = 5e-4;
	float min_hwep = 0;
	int hwe_midp = 0;
	
	// parse command-line arguments
	static struct option long_options[] =
	{
		{"min_hwe_pval",  required_argument, 0, 'i'},
		{"hwe_midp",  no_argument, 0, 'H'},
		{"prediction",  no_argument, 0, 'P'},
		{"snp_estimate_file",  required_argument, 0, 'a'},
		{"beta_weight_parameters",  required_argument, 0, 'w'},
		{"min_maf",  required_argument, 0, 'j'},
		{"phenotype_file",  required_argument, 0, 'p'},
		{"trait_name",  required_argument, 0, 't'},
		{"error_weight_name",  required_argument, 0, 'e'},
		{"covariate_file",  required_argument, 0, 'c'},
		{"covariate_names",  required_argument, 0, 'v'},
		{"snp_info_file",  required_argument, 0, 'm'},
		{"snp_set_name",  required_argument, 0, 'x'},
		{"snp_weight_name",  required_argument, 0, 'y'},
		{"binary_genotype_file",  required_argument, 0, 'g'},
		{"bfile",  required_argument, 0, 'b'},
		{"num_threads",    required_argument, 0, 'u'},
		{"output_file",    required_argument, 0, 'o'},
		{"subset_size",    required_argument, 0, 'z'},
		{"rel_tol",    required_argument, 0, 'r'},
		{"reml",  no_argument, 0, 'R'},
		{"num_random_vectors",  required_argument, 0, 'k'},
		{"lmm",  no_argument, 0, 'L'},
		{"num_qf_markers",  required_argument, 0, 'f'},
		{"window_size",  required_argument, 0, 'q'},
		{"max_heritability",  required_argument, 0, 'h'},
		{"fake_geno",  no_argument, 0, 'F'},
		{"seed",  required_argument, 0, 'd'},
		{"wls",  no_argument, 0, 'S'},
		{"lrt",  no_argument, 0, 'T'},
		{"iter_weighting",  no_argument, 0, 'I'},
		{"sig_chisq",  required_argument, 0, 's'},
		{"indep_r2",  required_argument, 0, 'n'},
		{0, 0, 0, 0}
	};
	int option_index = 0;
	int opt;
	while ((opt = getopt_long(option_num, option_str, "HPRLFSTIi:a:w:j:p:t:e:c:v:m:x:y:g:b:u:o:z:r:k:f:q:h:d:s:n:",long_options, &option_index)) != -1)
	{
		switch (opt) {
			case 'i':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				min_hwep = std::stod(optarg);
				break;
			case 'H':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				hwe_midp = 1;
				break;
			case 'w':
			{
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				const std::vector<std::string> beta_optarg = StrFunc::split(optarg,',');
				if(beta_optarg.size() != 2)
					throw("\nError: incorrect option for beta weight.\n");
				beta_pdf_params = std::make_pair(std::stod(beta_optarg[0]), std::stod(beta_optarg[1]));
				break;
			}
			case 'j':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				min_maf = std::stod(optarg);
				break;
			case 'P':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				pred_flag = true;
				break;
			case 'a':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_estimate_file = optarg;
				break;
			case 'p':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				phen_file = optarg;
				break;
			case 't':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				trait_name = optarg;
				break;
			case 'e':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				error_variance_weight_name = optarg;
				break;
			case 'c':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				covar_file = optarg;
				break;
			case 'v':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				str_covar_names = optarg;
				break;
			case 'm':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_info_file = optarg;
				break;
			case 'x':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_group = optarg;
				break;
			case 'y':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				marker_effvar_weight_name = optarg;
				break;
			case 'g':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				binary_genotype_file = optarg;
				break;
			case 'b':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				binary_genotype_file = optarg;
				break;
			case 'u':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				num_threads = std::stoi(optarg);
				break;
			case 'o':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				output_file = optarg;
				break;
			case 'z':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				subset_size = std::stoi(optarg);
				break;
			case 'r':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				rel_tol = std::stod(optarg);
				break;
			case 'R':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				reml_flag = true;
				break;
			case 'k':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				num_rand_prob = std::stoi(optarg);
				using_nr_opt = true;
				break;
			case 'L':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				lmm_flag = true;
				break;
			case 'f':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				num_qf_markers = std::stoi(optarg);
				break;
			case 'q':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				window_size = std::stof(optarg);
				using_ws_opt = true;
				break;
			case 'h':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				max_h2 = std::stof(optarg);
				break;
			case 'F':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				fake_geno = true;
				break;
			case 'd':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				rng_seed = std::stoi(optarg);
				break;
			case 'S':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				wls_flag = true;
				break;
			case 'T':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				lrt_flag = true;
				break;
			case 'I':
				printf("OPTION %s initialized\n", long_options[option_index].name);
				iter_weighting = true;
				break;
			case 's':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				sig_chisq = std::stof(optarg);
				break;
			case 'n':
				printf("OPTION %s with ARG %s\n", long_options[option_index].name,optarg);
				indep_r2 = std::stof(optarg);
				break;
			default:
				print_help();
				exit(1);
		}
	}
	
	if(!(pred_flag || reml_flag || lmm_flag || wls_flag)) throw("\nError: --lmm OR --reml OR --pred OR --wls is missing.\n");
	
	// predict total genetic values
	if(pred_flag) {
		std::cout<<"Calculating total genetic values."<<std::endl;
		if(binary_genotype_file.empty()) throw("\nError: --binary_genotype_file is required for analysis.\n");
		if(marker_estimate_file.empty()) throw("\nError: --snp_estimate_file is required for analysis.\n");
		if(output_file.empty()) throw("\nError: --output_file is required for analysis.\n");
		cal_dgv(marker_estimate_file, binary_genotype_file, output_file);
		std::cout<<"Calculation of total genetic values finished."<<std::endl;
		return;
	}
	
	// check options and associated options
	if(phen_file.empty()) throw("\nError: --phenotype_file is required for analysis.\n");
	if(trait_name.empty()) throw("\nError: --trait_name is required for analysis.\n");
	if(!str_covar_names.compare("all") && covar_file.empty()) 
		throw("\nError: --covariate_file is required when --covariate_names all is given.\n");
	if(covar_file.empty()) covar_file = phen_file;
	if(!str_covar_names.compare("all")) 
		covar_names = StrFunc::get_token_names(covar_file);
	else if(!str_covar_names.empty())
		covar_names = StrFunc::split(str_covar_names,',');
	if(marker_info_file.empty()) throw("\nError: --snp_info_file is required for analysis.\n");
	if(binary_genotype_file.empty()) throw("\nError: --binary_genotype_file is required for analysis.\n");
	if(min_maf < 0 || min_maf >= 0.5) throw("\nError: --min_maf is out of [0, 0.5).\n");
	if(min_hwep < 0 || min_hwep >= 1) throw("\nError: --min_hwe_pval is out of [0, 1).\n");
	if(output_file.empty()) throw("\nError: --output_file is required for analysis.\n");
	if(num_threads > 100 || num_threads < 1) throw("\nError: --num_threads is out of [1,100].\n");
	if(subset_size<10 || subset_size>5000) throw("\nError: --subset_size is out of [10, 5000].\n");
	
	if(num_rand_prob<3 || num_rand_prob>3000) throw("\nError: --num_rand_prob is out of [3, 3000].\n");
	if(num_qf_markers<3 || num_qf_markers>300) throw("\nError: --num_qf_markers is out of [3, 300].\n");
	
	// change the default window size if not LMM (added Oct 30, 2022)
	if((! lmm_flag) && (! using_ws_opt)) window_size = 20;
	if(window_size % 2 != 0) {
		std::cout<<"\nWarning: --window_size is reset to its clostest even number.\n";
		window_size = window_size/2*2;
	}
	if(lmm_flag) {
		if(window_size<10 || window_size>5000) throw("\nError: --window_size is out of [10, 5000] for --lmm.\n");
	} else if(window_size<0) throw("\nError: --window_size is negative.\n");
	
	if((reml_flag || lmm_flag) && !marker_group.empty()) {
		std::cout<<"\nWarning: --snp_set_name is temporarily disabled in --reml or --lmm.\n";
		marker_group = "";
	}
	// set threads
	omp_set_num_threads(num_threads);
	// input data
	// read phenotype file
	std::map<std::string, std::pair<float, float> > indi2pheno_weight = read_phenotype_file(phen_file, trait_name, error_variance_weight_name);
	// read covariate file
	std::map<std::string, VectorXf > indi2covar = read_covariate_file(covar_file, covar_names);
	// read marker prior file
	std::map<std::string, std::pair<std::string, float> > marker2group_weight;
	read_marker_info_file(marker_info_file, marker_group, marker_effvar_weight_name, marker2group_weight);
	// get dimensions for analysis
	unsigned indi_keep_num, marker_num, group_keep_num;
	if(file_check(binary_genotype_file + ".fam") && file_check(binary_genotype_file + ".bim") && file_check(binary_genotype_file + ".bed"))
		get_dim_from_plink(binary_genotype_file, indi2pheno_weight, indi2covar, marker2group_weight, indi_keep_num, marker_num, group_keep_num );
	else throw("\nError: bed/bim/fam fileset NOT found.\n");
	std::map<std::string, float> marker2maf;
	std::map<std::string, float> marker2hwep;
	std::map<std::string, std::pair<std::string, int> > marker2pos;
	std::map<std::string, std::pair<std::string, std::string> > marker2alleles;
	std::vector<std::string> indi_keep;
	std::vector<std::string> marker_keep;
	std::vector<std::string> group_keep;
	VectorXi group_sindex;
	VectorXi group_ncols;
	std::vector<VectorXf> gdiag;
	MatrixXf xmat;
	VectorXf yvec;
	VectorXf rvec;
	yvec.resize(indi_keep_num);
	rvec.resize(indi_keep_num);
	xmat.resize(indi_keep_num, covar_names.size());
	std::vector<std::vector<bool> > kmat(marker_num, std::vector<bool>(indi_keep_num*2, true));
	group_sindex.resize(group_keep_num);
	group_ncols.resize(group_keep_num);
	
	read_plink_bin_geno_file_and_generate_pxvb_data(binary_genotype_file, indi2pheno_weight, indi2covar, marker2group_weight,
	  beta_pdf_params, min_maf, marker2maf, min_hwep, hwe_midp, marker2hwep, marker2pos, marker2alleles, 
	  indi_keep, marker_keep, group_keep, group_sindex, group_ncols, gdiag, xmat, kmat, yvec, rvec );
	
	group_keep_num = group_keep.size();
	group_sindex.conservativeResize(group_keep_num);
	group_ncols.conservativeResize(group_keep_num);
	
	if(marker_keep.empty()) throw("\nError: no SNPs are kept in analysis.\n");
	
	// WLS
	if(wls_flag) {
		VectorXf Py(indi_keep_num); Py.setZero();
		VectorXf snp_blup(marker_keep.size());
		VectorXf blue(xmat.cols());  // not used by --wls
		MatrixXf var_blue(0, 0);     // not used by --wls
		float vg=0, ve=0, llr=0;
		
		weighted_least_squares(rvec, kmat, xmat, yvec, subset_size, snp_blup);
		
		write_reml_into_files(output_file, covar_names, indi_keep, rvec, xmat, marker2maf, marker2hwep, marker2pos, marker2alleles, group_ncols, group_keep, marker_keep,
		  Py, snp_blup, blue, var_blue, vg, ve, llr);
		
		return;
	}
	
	// LMM assoc using L-FOMC-REML
	if(lmm_flag) {
		if(marker_keep.size() < window_size+1) throw("\nError: --window_size >= # of SNPs kept for --lmm.\n");
		
		int n_MC = std::max(std::min(int(4e9/indi_keep_num/indi_keep_num), 30), 3);
		
		VectorXf persnp = gdiag[0];
		gdiag[0].resize(0);
		
		VectorXf Py(indi_keep_num);
		VectorXf snp_blup(marker_keep.size());
		VectorXf blue(xmat.cols());
		MatrixXf var_blue(0, 0);     // not used by --lmm
		float vg, ve, llr=0;
		std::vector<gstat> lmm_out;
		LMC_REML(rvec, kmat, xmat, yvec, persnp, subset_size, rng_seed, n_MC, rel_tol, 
		  lmm_flag, num_qf_markers, fake_geno, window_size, Py, snp_blup, blue, vg, ve, lmm_out, max_h2);
		
		write_reml_into_files(output_file, covar_names, indi_keep, rvec, xmat, marker2maf, marker2hwep, marker2pos, marker2alleles, group_ncols, group_keep, marker_keep,
		  Py, snp_blup, blue, var_blue, vg, ve, llr);
		write_gstat_into_files(output_file, indi_keep_num, window_size, vg, marker2maf, marker2hwep, marker2pos, marker_keep, lmm_out);
		
		return;
	}
	
	// REML using SLDF (for --lrt) or L-FOMC
	if(reml_flag) {
		VectorXf persnp = gdiag[0];
		gdiag[0].resize(0);
		
		VectorXf Py(indi_keep_num);
		VectorXf snp_blup(marker_keep.size());
		VectorXf blue(xmat.cols());
		MatrixXf var_blue;
		float vg, ve, llr=0;
		
		// L-FOMC-specific parameters
		int n_MC = std::max(std::min(int(4e9/indi_keep_num/indi_keep_num), 30), 3);
		std::vector<gstat> lmm_out; // not used
		
		if(iter_weighting && (persnp.size() < window_size+1)) {
			std::cout<<"\nWarning: --window_size >= total # of SNPs in --iter_weighting.\n";
			std::cout<<"--iter_weighting is thus disabled.\n";
			iter_weighting = false;
		}
		
		if(! iter_weighting) {
			persnp /= persnp.mean();
			std::cout<<"\nSNP weights for heritability: "<<persnp.minCoeff()<<" - "<<persnp.maxCoeff()<<"\n";
			
			if(lrt_flag) {
				var_blue.resize(xmat.cols(), xmat.cols());
				SLDF_REML(rvec, kmat, xmat, yvec, persnp, subset_size, rng_seed, num_rand_prob, rel_tol, 
					Py, snp_blup, blue, var_blue, vg, ve, llr, max_h2);
			} else {
				if(using_nr_opt) n_MC = num_rand_prob;
				LMC_REML(rvec, kmat, xmat, yvec, persnp, subset_size, rng_seed, n_MC, rel_tol, 
					false, num_qf_markers, fake_geno, window_size, Py, snp_blup, blue, vg, ve, lmm_out, max_h2);
			}
		} else {
			std::cout<<"\nIterative SNP weighting with window size = "<<window_size<<"\n";
			VectorXf snpvar = VectorXf::Ones(persnp.size()); 
			VectorXf tmp_persnp = VectorXf::Ones(persnp.size()); 
			
			std::cout<<"\n==========================================================\n";
			std::cout<<"\nEqual-weight GREML started.\n";
			std::cout<<"\nSNP weights for heritability: "<<tmp_persnp.minCoeff()<<" - "<<tmp_persnp.maxCoeff()<<"\n";
			
			// --lrt not used during window weighting for large samples (TODO: how large?)
			if(lrt_flag && indi_keep.size() < 5000) {
				var_blue.resize(xmat.cols(), xmat.cols());
				SLDF_REML(rvec, kmat, xmat, yvec, tmp_persnp, subset_size, rng_seed, num_rand_prob, rel_tol, 
					Py, snp_blup, blue, var_blue, vg, ve, llr, max_h2);
			} else {
				LMC_REML(rvec, kmat, xmat, yvec, tmp_persnp, subset_size, rng_seed, n_MC, rel_tol, 
					false, num_qf_markers, fake_geno, window_size, Py, snp_blup, blue, vg, ve, lmm_out, max_h2);
			}
			
			std::cout<<"\nEqual-weight GREML completed.\n";
			std::cout<<"\n==========================================================\n";
			
			for(unsigned i=0; i<persnp.size(); ++i) {
				float maf = marker2maf[marker_keep[i]];
				snpvar(i) = 2.0 * maf * (1.0 - maf) * snp_blup(i) * snp_blup(i);
			}
			std::cout<<"\nComputed each SNP's genetic variance by SNP effects.\n";
			
			snpvar /= snpvar.mean();
			tmp_persnp = snpvar;
			
			std::vector<unsigned> sigsnp_idx;
			MatrixXf sigsnp_X(indi_keep_num, 1);
			VectorXf sigsnp_y(indi_keep_num);
			VectorXf sigsnp_y_res(indi_keep_num);
			unsigned max_idx;
			while(1) {
				if(tmp_persnp.maxCoeff(&max_idx) > sig_chisq) { // key parameter 1 in thining/clumping
					tmp_persnp(max_idx) = 0;
					for(int k=0, kk=0; kk<indi_keep_num; k+=2, ++kk) 
						sigsnp_y(kk) = (kmat[max_idx][k] + kmat[max_idx][k+1]);
					sigsnp_y.array() -= sigsnp_y.mean();
					
					if(sigsnp_idx.empty()) {
						sigsnp_idx.push_back(max_idx);
						sigsnp_X.col(0) = sigsnp_y;
					} else {
						sigsnp_y_res = sigsnp_y - sigsnp_X * (sigsnp_X.transpose()*sigsnp_X).llt().solve(sigsnp_X.transpose()*sigsnp_y);
						if(sigsnp_y_res.squaredNorm()/sigsnp_y.squaredNorm() > indep_r2) { // key parameter 2 in thining/clumping
							sigsnp_idx.push_back(max_idx);
							sigsnp_X.conservativeResize(NoChange, sigsnp_idx.size());
							sigsnp_X.col(sigsnp_idx.size()-1) = sigsnp_y;
						}
					}
				} else {
					break;
				}
			}
			sigsnp_y.resize(0);
			sigsnp_y_res.resize(0);
			
			if(! sigsnp_idx.empty()) {
				xmat.conservativeResize(NoChange, covar_names.size()+sigsnp_idx.size() );
				xmat.rightCols(sigsnp_idx.size()) = sigsnp_X;
			}
			sigsnp_X.resize(0,0);
			
			std::cout<<"Identified independent significant SNPs.\n";
			
			for(unsigned i=0; i<persnp.size(); ++i) {
				if(i < window_size/2) tmp_persnp(i) = snpvar.head(window_size+1).mean();
				else if(i> persnp.size() - window_size/2 -1) tmp_persnp(i) = snpvar.tail(window_size+1).mean(); 
				else tmp_persnp(i) = snpvar.segment(i - window_size/2, window_size+1).mean();
			}
			// Added August 15, 2024
			tmp_persnp = 0.001 / (0.001 + 0.999 * (tmp_persnp * (-0.5)).array().exp());
			
			std::cout<<"Computed window-based SNP weights.\n";
			std::cout<<"\n==========================================================\n";
			
			std::cout<<"\nWindow-weighted GREML started.\n";
			
			if(sigsnp_idx.empty()) 
				std::cout<<"\nNo SNPs will be fitted as fixed effects.\n";
			else {
				std::cout<<"\nThe following SNPs will be fitted as fixed effects:\n";
				for(unsigned i=0; i<sigsnp_idx.size(); ++i) {
					tmp_persnp(sigsnp_idx[i]) = 0;
					std::cout<<marker_keep[sigsnp_idx[i]]<<"\n";
				}
			}
			
			persnp.array() *= tmp_persnp.array();
			persnp /= persnp.mean();
			std::cout<<"\nSNP weights for heritability: "<<persnp.minCoeff()<<" - "<<persnp.maxCoeff()<<"\n";
			
			blue.resize(xmat.cols());
			if(lrt_flag) {
				var_blue.resize(xmat.cols(), xmat.cols());
				SLDF_REML(rvec, kmat, xmat, yvec, persnp, subset_size, rng_seed, num_rand_prob, rel_tol, 
					Py, snp_blup, blue, var_blue, vg, ve, llr, max_h2);
			} else {
				if(using_nr_opt) n_MC = num_rand_prob;
				LMC_REML(rvec, kmat, xmat, yvec, persnp, subset_size, rng_seed, n_MC, rel_tol, 
					false, num_qf_markers, fake_geno, window_size, Py, snp_blup, blue, vg, ve, lmm_out, max_h2);
			}
			for(unsigned i=0; i<sigsnp_idx.size(); ++i) snp_blup(sigsnp_idx[i]) = blue(i+covar_names.size());
			
			std::cout<<"\nWindow-weighted GREML completed.\n";
		}
		
		xmat.conservativeResize(NoChange, covar_names.size());
		blue.conservativeResize(covar_names.size());
		if(lrt_flag) var_blue.conservativeResize(covar_names.size(), covar_names.size());
		write_reml_into_files(output_file, covar_names, indi_keep, rvec, xmat, marker2maf, marker2hwep, marker2pos, marker2alleles, group_ncols, group_keep, marker_keep,
		  Py, snp_blup, blue, var_blue, vg, ve, llr);
		
		return;
	}
}

void print_help() {
	printf("\nPlease see online documentation at https://github.com/jiang18/slemm\n");
}
