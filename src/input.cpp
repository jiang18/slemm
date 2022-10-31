#include <bitset>
#include "slemm.h"
#include "snphwe.h"

const std::string default_group ("NULL");
const float default_marker_effvar_weight = 1.0;
const float default_error_variance_weight = 1.0;
const char delimiter = ',';

void cal_hwep_from_genobits(const std::vector<std::vector<bool> >& geno, const int midp, Ref<VectorXf> hwep)
{
	#pragma omp parallel for
	for(unsigned i=0; i<geno.size(); ++i) {
		int n0 = 0, n1 = 0, n2 = 0;
		for(unsigned j=0; j<geno[i].size(); j+=2) {
			if(geno[i][j] & geno[i][j+1]) n2++;
			else if(geno[i][j] | geno[i][j+1]) n1++;
			else n0++;
		}
		hwep(i) = SNPHWE2(n1, n0, n2, midp);
	}
}

// Remember to take advantage of Return Value Optimization in the main function.
std::map<std::string, std::pair<float, float> > read_phenotype_file(std::string phenotype_file, std::string trait_name, std::string error_variance_weight_name)
{
	std::map<std::string, std::pair<float, float> > indi2pheno_weight;
	std::ifstream ifs;
	std::string line;
	
	ifs.open(phenotype_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<phenotype_file<<std::endl;
		exit(1);
	}
	
	std::getline(ifs, line);
	if(line[line.size()-1] == '\r') line.erase(line.size()-1);
	std::vector<std::string> header = StrFunc::split(line, delimiter);
	
	int trait_col = StrFunc::find_token_in_header(trait_name, header, phenotype_file);
	int error_variance_weight_col = StrFunc::find_token_in_header(error_variance_weight_name, header, phenotype_file);
	// Beware of the situation where the last column is missing.
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::vector<std::string> cols = StrFunc::split(line, delimiter);
		if(!(header.size() == cols.size() || (line.back() == delimiter && header.size() == cols.size()+1))) {
			std::cout<<"\nError: "<<line<<" in "<<phenotype_file<<std::endl;
			exit(1);
		}
		if(error_variance_weight_col == na_col && cols.size() > trait_col && !cols[trait_col].empty())
			indi2pheno_weight[cols[0]] = std::make_pair(std::stod(cols[trait_col]), default_error_variance_weight);
		else if(cols.size() > error_variance_weight_col && !cols[error_variance_weight_col].empty() && cols.size() > trait_col && !cols[trait_col].empty())
			indi2pheno_weight[cols[0]] = std::make_pair(std::stod(cols[trait_col]), std::stod(cols[error_variance_weight_col]));
	}
	ifs.close();
	
	return indi2pheno_weight;
}

std::map<std::string, VectorXf> read_covariate_file(std::string covariate_file, std::vector<std::string>& covariate_names)
{
	std::map<std::string, VectorXf > indi2covar;
	int i;
	std::ifstream ifs;
	std::string line;
	
	ifs.open(covariate_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<covariate_file<<std::endl;
		exit(1);
	}
	
	std::map<std::string, int> name2col; 
	std::getline(ifs, line);
	if(line[line.size()-1] == '\r') line.erase(line.size()-1);
	std::vector<std::string> header = StrFunc::split(line, delimiter);
	std::cout<<"Header in covariate file: "<<std::endl;
	for (i=1; i<header.size(); ++i) {
		std::cout<<" "<<header[i];
		name2col[header[i]] = i;
	}
	std::cout<<std::endl;
	std::vector<int> covariate_cols;
	std::cout<<"Covariates in command line option: "<<std::endl;
	for (i=0; i<covariate_names.size(); ++i) {
		std::cout<<" "<<covariate_names[i];
		if(name2col.find(covariate_names[i]) == name2col.end()) {
			std::cout<<"\n"<<covariate_names[i]<<" not in covariate file"<<std::endl;
			exit(1);
		} else covariate_cols.push_back(name2col[covariate_names[i]]);
	}
	std::cout<<std::endl;
	VectorXf covars (covariate_cols.size());
	// Beware of the situation where the last column is missing.
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::vector<std::string> cols = StrFunc::split(line, delimiter);
		if(!(header.size() == cols.size() || (line.back() == delimiter && header.size() == cols.size()+1))) {
			std::cout<<"\nError: "<<line<<" in "<<covariate_file<<std::endl;
			exit(1);
		}
		if(covariate_cols.empty()) indi2covar[cols[0]] = VectorXf::Ones(1);
		else {
			for (i=0; i<covariate_cols.size(); ++i) {
				if(cols.size() > covariate_cols[i] && !cols[covariate_cols[i]].empty()) covars[i] = std::stod(cols[covariate_cols[i]]);
				else goto nextline;
			}
			indi2covar[cols[0]] = covars;
		}
		nextline:;
	}
	ifs.close();
	if(covariate_cols.empty()) covariate_names.push_back("intercept");
	return indi2covar;
}

void read_marker_info_file(
	std::string marker_info_file, 
	std::string group_header, 
	std::string weight_name, 
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight )
{
	marker2group_weight.clear();
	
	std::ifstream ifs;
	std::string line;
	
	ifs.open(marker_info_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<marker_info_file<<std::endl;
		exit(1);
	}
	
	std::getline(ifs, line);
	if(line[line.size()-1] == '\r') line.erase(line.size()-1);
	std::vector<std::string> header = StrFunc::split(line, delimiter);
	
	int group_col = StrFunc::find_token_in_header(group_header, header, marker_info_file);
	if(!group_header.empty() && group_col == na_col) {
		std::cout<<"Marker group column "<<group_header<<" not found"<<std::endl;
		exit(1);
	}
	int weight_col = StrFunc::find_token_in_header(weight_name, header, marker_info_file);
	if(!weight_name.empty() && weight_col == na_col) {
		std::cout<<"Marker weight column "<<weight_name<<" not found"<<std::endl;
		exit(1);
	}
	
	while (std::getline(ifs, line)) {
		if(line[line.size()-1] == '\r') line.erase(line.size()-1);
		std::vector<std::string> cols = StrFunc::split(line, delimiter);
		if(!(header.size() == cols.size() || (line.back() == delimiter && header.size() == cols.size()+1))) {
			std::cout<<"\nError: "<<line<<" in "<<marker_info_file<<std::endl;
			exit(1);
		}
		if(group_col == na_col) {
			if(weight_col == na_col)
				marker2group_weight[cols[0]] = std::make_pair(default_group, default_marker_effvar_weight);
			else if(cols.size() > weight_col && !cols[weight_col].empty())
				marker2group_weight[cols[0]] = std::make_pair(default_group, std::stod(cols[weight_col]));
			else continue;
		} else if(cols.size() > group_col && !cols[group_col].empty()) {
			if(weight_col == na_col)
				marker2group_weight[cols[0]] = std::make_pair(cols[group_col], default_marker_effvar_weight);
			else if(cols.size() > weight_col && !cols[weight_col].empty())
				marker2group_weight[cols[0]] = std::make_pair(cols[group_col], std::stod(cols[weight_col]));
			else continue;
		}
	}
	ifs.close();
}

void get_dim_from_plink(
	const std::string binary_genotype_file_prefix, 
	const std::map<std::string, std::pair<float, float> >& indi2pheno_weight, 
	const std::map<std::string, VectorXf >& indi2covar, 
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight, 
	unsigned& indi_keep_num,
	unsigned& marker_num, 
	unsigned& group_keep_num )
{
	std::ifstream ifs;
	// read individual list file
	std::string indi_file = binary_genotype_file_prefix + ".fam";
	ifs.open(indi_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<indi_file<<std::endl;
		exit(1);
	}
	indi_keep_num=0;

	std::stringstream ss;
	std::string line;
	std::string pid;
	while( std::getline(ifs, line) ){
	// line content: Family_ID Within-family_ID Within-family_ID_of_father Within-family_ID_of_mother Sex Case-control_phenotype
		if(line.size() < 3) continue;
		ss.str(line);
		ss >> pid >> pid;
		if(indi2pheno_weight.find(pid) != indi2pheno_weight.end() && indi2covar.find(pid) != indi2covar.end()) 
			indi_keep_num ++;
	}
	ifs.clear();
	ifs.close();
	if(indi_keep_num == 0) {
		std::cerr<<"Individual set retained for analysis is empty."<<std::endl;
		exit(1);
	}
	// read marker list file
	std::string marker_file = binary_genotype_file_prefix + ".bim";
	ifs.open(marker_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<marker_file<<std::endl;
		exit(1);
	}
	marker_num = 0;
	std::map<std::string, bool > group_keep;
	int ibuf=0;
	std::string str_buf;
	std::string chr_buf;
	std::string snp_buf;
	while(ifs) {
		ifs>>chr_buf;  // chromosome
		if(ifs.eof()) break;
		ifs>>snp_buf;  // variant name
		ifs>>str_buf;  // genetic dist
		ifs>>ibuf;     // bp position
		ifs>>str_buf;  // allele 1
		ifs>>str_buf;  // allele 2
		
		if(marker2group_weight.find(snp_buf) != marker2group_weight.end()) {
			group_keep[marker2group_weight[snp_buf].first] = true;
			marker_num ++;
		}
	}
	ifs.clear();
	ifs.close();
	group_keep_num = group_keep.size();
	if(marker_num == 0) {
		std::cerr<<"Marker set retained for analysis is empty."<<std::endl;
		exit(1);
	}
}

// some code are adopted from GCTA with modifications
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
	Ref<VectorXf> rvec)
{
	std::ifstream ifs;
	unsigned i=0;
	// read individual list file
	std::string indi_file = binary_genotype_file_prefix + ".fam";
	ifs.open(indi_file);
	std::cout<<"\nReading PLINK FAM file from ["+indi_file+"]."<<std::endl;
	std::vector<bool> bindi;
	indi_keep.clear();
	std::stringstream ss;
	std::string line;
	std::string pid;
	while( std::getline(ifs, line) ){
	// line content: Family_ID Within-family_ID Within-family_ID_of_father Within-family_ID_of_mother Sex Case-control_phenotype
		if(line.size() < 3) continue;
		ss.str(line);
		ss >> pid >> pid;
		if(indi2pheno_weight.find(pid) != indi2pheno_weight.end() && indi2covar.find(pid) != indi2covar.end()) {
			indi_keep.push_back(pid);
			bindi.push_back(true);
		}
		else bindi.push_back(false);
	}
	ifs.clear();
	ifs.close();
	if(indi_keep.empty()) {
		std::cerr<<"Individual set is empty."<<std::endl;
		exit(1);
	}
	std::cout<<"Number of individuals retained in analysis = "<<indi_keep.size()<<std::endl;
	for(i=0; i<indi_keep.size(); ++i) {
		yvec[i] = indi2pheno_weight[indi_keep[i]].first;
		rvec[i] = indi2pheno_weight[indi_keep[i]].second;
		xmat.row(i) = indi2covar[indi_keep[i]];
	}
	// read marker list file
	std::string marker_file = binary_genotype_file_prefix + ".bim";
	ifs.open(marker_file);
	int ibuf=0;
	std::string str_buf;
	std::string chr_buf;
	std::string snp_buf;
	std::string a1_buf, a2_buf;
	std::cout<<"\nReading PLINK BIM file from ["+marker_file+"]."<<std::endl;
	marker2pos.clear();
	std::vector<std::string> marker_set_in_kmat;
	std::map<std::string, unsigned> marker2index;
	std::vector<bool> bmarker;
	i = 0;
	while(ifs) {
		ifs>>chr_buf;  // chromosome
		if(ifs.eof()) break;
		ifs>>snp_buf;  // variant name
		ifs>>str_buf;  // genetic dist
		ifs>>ibuf;     // bp position
		ifs>>a1_buf;  // allele 1
		ifs>>a2_buf;  // allele 2
		
		if(marker2group_weight.find(snp_buf) != marker2group_weight.end()) {
			marker2pos[snp_buf] = std::make_pair(chr_buf, ibuf);
			marker2alleles[snp_buf] = std::make_pair(a1_buf, a2_buf);
			marker_set_in_kmat.push_back(snp_buf);
			marker2index[snp_buf] = i++;
			bmarker.push_back(true);
		} else bmarker.push_back(false);
	}
	ifs.clear();
	ifs.close();
	if(marker_set_in_kmat.empty()) {
		std::cerr<<"Marker set for analysis is empty."<<std::endl;
		exit(1);
	}
	std::cout<<"Size of marker set before MAF filtering = "<<marker_set_in_kmat.size()<<std::endl;
	//read binary genotype file
	std::string binary_genotype_file = binary_genotype_file_prefix + ".bed";
	ifs.open(binary_genotype_file, std::ifstream::binary);
	std::cout<<"\nReading PLINK BED file from ["+binary_genotype_file+"] in SNP-major format..."<<std::endl;
	unsigned snp_num = bmarker.size();
	unsigned indi_num = bindi.size();
	unsigned j=0, k=0;
	char ch[1];
	std::bitset<8> b;
	for(i=0; i<3; i++) ifs.read(ch,1); // skip the first three bytes
	// Read genotype in SNP-major mode
	// 00  Homozygous for first allele in .bim file
	// 11  Homozygous for second allele in .bim file
	// 01  Missing genotype
	// 10  Heterozygous
	unsigned snp_indx=0, indi_indx=0;
	for(j=0, snp_indx=0; j<snp_num; j++){
		if(!bmarker[j]){
			for(i=0; i<indi_num; i+=4) ifs.read(ch,1);
			continue;
		}
		for(i=0, indi_indx=0; i<indi_num;){
			ifs.read(ch,1);
			if(!ifs) throw("\nError: problem with the BED file. Has the FAM/BIM file been changed?\n");
			b=ch[0];
			k=0;
			while(k < 7 && i < indi_num){
				if(!bindi[i]) k+=2;
				else{
					// kmat[snp_indx][indi_indx++] = b[k++];
					// kmat[snp_indx][indi_indx++] = b[k++];
					if(b[k++]) indi_indx++;
					else kmat[snp_indx][indi_indx++] = false;
					if(b[k++]) indi_indx++;
					else kmat[snp_indx][indi_indx++] = false;
				}
				i++;
			}
		}
		snp_indx++;
		if(snp_indx==marker_set_in_kmat.size()) break;
	}
	ifs.clear();
	ifs.close();
	
	// calculate minor allele frequency
	std::cout<<"Calculating minor allele frequency..."<<std::endl;
	VectorXf maf(marker_set_in_kmat.size());
	cal_maf_from_genobits(kmat, maf);
	std::cout<<(maf.array() <= min_maf).count()<<" markers of minor allele frequency <= "<<min_maf<<" are removed."<<std::endl;
	VectorXf hwep(marker_set_in_kmat.size());
	cal_hwep_from_genobits(kmat, hwe_midp, hwep);
	std::cout<<(hwep.array() < min_hwep).count()<<" markers fail in HWE tests: P<"<<min_hwep<<std::endl;
	
	std::map<std::string, std::vector<std::string> > group2marker;
	marker2maf.clear();
	marker2hwep.clear();
	for(i=0; i<marker_set_in_kmat.size(); ++i) {
		if(maf[i] > min_maf && hwep[i] >= min_hwep) {
			std::string this_marker = marker_set_in_kmat[i];
			group2marker[marker2group_weight[this_marker].first].push_back(this_marker);
			marker2maf[this_marker] = maf[i];
			marker2hwep[this_marker] = hwep[i];
		}
	}
	if(marker2maf.empty()) {
		std::cerr<<"Marker set is empty."<<std::endl;
		exit(1);
	}
	// update marker weights using minor allele frequency
	if(beta_pdf_params.first != 1 && beta_pdf_params.second != 1) {
		update_weight_by_maf(beta_pdf_params, marker2maf, marker2group_weight);
	}

	// Set parameters for each marker group kept, including 1) the index of the first marker, 
	// 2) the number of markers, 3) shape and 4) scale.
	// Get the markers and groups to be used
	// Set weights of effect variance for each group of variants
	i = 0;
	j = 0;
	k = 0;
	marker_keep.clear();
	group_keep.clear();
	gdiag.clear();
	VectorXf marker_weight;
	for(auto it=group2marker.begin(); it!=group2marker.end(); ++it) {
		marker_keep.insert(marker_keep.end(), it->second.begin(), it->second.end());
		group_keep.push_back(it->first);
		group_sindex[i] = j;
		group_ncols[i] = it->second.size();
		
		marker_weight.resize(it->second.size());
		for(k=0; k<it->second.size(); ++k) {
			marker_weight[k] = marker2group_weight[it->second[k]].second;
		}
		gdiag.push_back(marker_weight);
		
		i ++;
		j += it->second.size();
	}
	std::cout<<"Number of markers retained = "<<marker_keep.size()<<std::endl;
	std::cout<<"Number of groups retained = "<<group_keep.size()<<std::endl;
	// Sort the columns of K matrix according to the markers retained
	for(i=0; i<marker_keep.size(); ++i) {
		j = marker2index[marker_keep[i]];
		if(j == i) continue;
		kmat[i].swap(kmat[j]);
		
		marker_set_in_kmat[j] = marker_set_in_kmat[i];
		marker2index[marker_set_in_kmat[j]] = j;
	}
	kmat.resize(marker_keep.size());
}

