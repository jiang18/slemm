#include <bitset>
#include "slemm.h"

void computeTotalGeneticValueSSGP(
	const std::string marker_estimate_file,
	const std::string binary_genotype_file_prefix,
	const std::string output_file)
{
	std::ifstream ifs;
	std::string line;
	unsigned i=0;
	// marker effects
	ifs.open(marker_estimate_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<marker_estimate_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading SNP ESTIMATE file from ["+marker_estimate_file+"]."<<std::endl;
	std::getline(ifs, line);
	std::map<std::string, double> marker2eff;
	std::map<std::string, std::string> marker2allele;
	while (std::getline(ifs, line)) {
		std::vector<std::string> cols = StrFunc::split(line, ',');
		if(cols.size() < 10 || cols[0].empty() || cols[9].empty()) continue;
		marker2eff[cols[0]] = std::stod(cols[9]);
		marker2allele[cols[0]] = cols[3];
	}
	ifs.close();
	// read individual list file
	std::string indi_file = binary_genotype_file_prefix + ".indi";
	ifs.open(indi_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<indi_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading SSGP INDI file from ["+indi_file+"]."<<std::endl;
	std::vector<std::string> indi;
	while (std::getline(ifs, line)) indi.push_back(line);
	ifs.close();
	// read marker list file
	std::string marker_file = binary_genotype_file_prefix + ".mrk";
	ifs.open(marker_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<marker_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading SSGP MRK file from ["+marker_file+"]."<<std::endl;
	std::vector<std::string> marker;
	std::vector<std::string> allele;
	while (std::getline(ifs, line)) {
		std::vector<std::string> cols = StrFunc::split(line, ',');
		if(cols.size() < 5) throw("\nError: not enough columns in MRK file.\n");
		marker.push_back(cols[0]);
		allele.push_back(cols[3]);
	}
	ifs.close();
	//read binary genotype file
	std::string binary_genotype_file = binary_genotype_file_prefix + ".bin";
	ifs.open(binary_genotype_file, std::ifstream::binary);
	std::cout<<"Reading SSGP BIN file from ["+binary_genotype_file+"]."<<std::endl;
	unsigned data_n,data_m;
	ifs.read((char*)(&data_n),sizeof(unsigned));
	ifs.read((char*)(&data_m),sizeof(unsigned));
	//to check data_n == indi_num or data_m != marker_num
	if(data_n != indi.size() || data_m != marker.size()) {
		std::cout<<"# of individuals in .indi != that in .bin or # of markers in .mrk != that in .bin"<<std::endl;
		std::cout<<"Has one (or more) of the three files been modified?"<<std::endl;
		exit(1);
	}
	// g = Z * alpha
	VectorXd alpha(marker.size());
	for(i=0; i<marker.size(); ++i) {
		if(marker2eff.find(marker[i]) != marker2eff.end()) {
			if(allele[i] == marker2allele[marker[i]]) alpha[i] = marker2eff[marker[i]];
			else alpha[i] = -marker2eff[marker[i]];
		
		} else alpha[i] = 0;
	}
	std::ofstream ofs;
	ofs.open(output_file);
	Matrix<int8_t, Dynamic, 1> indi_geno(data_m);
	for(i=0; i<data_n; ++i)
	{
		ifs.read((char*)indi_geno.data(), data_m*sizeof(signed char));
		double tgv = indi_geno.cast<double>().dot(alpha);
		ofs<<indi[i]<<','<<tgv<<std::endl;
	}
	ifs.close();
	
}

void computeTotalGeneticValuePLINK(
	const std::string marker_estimate_file,
	const std::string binary_genotype_file_prefix,
	const std::string output_file)
{
	std::ifstream ifs;
	std::string line;
	// marker effects
	ifs.open(marker_estimate_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<marker_estimate_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading SNP ESTIMATE file from ["+marker_estimate_file+"]."<<std::endl;
	std::getline(ifs, line);
	std::map<std::string, double> marker2eff;
	std::map<std::string, std::string> marker2allele;
	while (std::getline(ifs, line)) {
		std::vector<std::string> cols = StrFunc::split(line, ',');
		if(cols.size() < 10 || cols[0].empty() || cols[9].empty()) continue;
		marker2eff[cols[0]] = std::stod(cols[9]);
		marker2allele[cols[0]] = cols[3];
	}
	ifs.close();
	// read individual list file
	std::string indi_file = binary_genotype_file_prefix + ".fam";
	ifs.open(indi_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<indi_file<<std::endl;
		exit(1);
	}
	std::cout<<"Reading PLINK FAM file from ["+indi_file+"]."<<std::endl;
	std::vector<std::string> indi;
	std::stringstream ss;
	std::string pid;
	while( std::getline(ifs, line) ){
	// line content: Family_ID Within-family_ID Within-family_ID_of_father Within-family_ID_of_mother Sex Case-control_phenotype
		if(line.size() < 3) continue;
		ss.str(line);
		ss >> pid >> pid;
		indi.push_back(pid);
	}
	ifs.clear();
	ifs.close();

	// read marker list file
	std::string marker_file = binary_genotype_file_prefix + ".bim";
	ifs.open(marker_file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<marker_file<<std::endl;
		exit(1);
	}
	std::vector<double> alpha;
	std::vector<bool> bmarker;
	std::string str_buf;
	std::string snp_buf;
	std::string a1_buf;
	std::cout<<"Reading PLINK BIM file from ["+marker_file+"]."<<std::endl;
	while(ifs) {
		ifs>>str_buf;  // chromosome
		if(ifs.eof()) break;
		ifs>>snp_buf;  // variant name
		ifs>>str_buf;  // genetic dist
		ifs>>str_buf;     // bp position
		ifs>>a1_buf;  // allele 1
		ifs>>str_buf;  // allele 2
		
		if(marker2eff.find(snp_buf) != marker2eff.end()) {
			double this_effect = marker2eff[snp_buf];
			if(a1_buf != marker2allele[snp_buf]) this_effect *= -1;
			alpha.push_back(this_effect);
			bmarker.push_back(true);
		} else {
			alpha.push_back(0.0);
			bmarker.push_back(false);
		}
	}
	ifs.clear();
	ifs.close();

	unsigned snp_num = bmarker.size();
	unsigned indi_num = indi.size();
	unsigned i=0, j=0, k=0;
	//get the index of last true SNP
	unsigned last_true_snp = -1;
	for(i=snp_num-1; i>=0; i--)
		if(bmarker[i] == true) {last_true_snp = i; break;}
	if(last_true_snp == -1) throw("\nNo SNPs found both in estimate file and in genotype file.\n");

	//read binary genotype file
	std::string binary_genotype_file = binary_genotype_file_prefix + ".bed";
	ifs.open(binary_genotype_file, std::ifstream::binary);
	std::cout<<"Reading PLINK BED file from ["+binary_genotype_file+"] in SNP-major format..."<<std::endl;
	// g = Z * alpha
	Matrix<int8_t, Dynamic, 1> zvec(indi_num);
	VectorXd tgv;
	tgv.setZero(indi_num);

	char ch[1];
	std::bitset<8> b;
	for(i=0; i<3; i++) ifs.read(ch,1); // skip the first three bytes
	// Read genotype in SNP-major mode
	// 00  Homozygous for first allele in .bim file
	// 11  Homozygous for second allele in .bim file
	// 01  Missing genotype
	// 10  Heterozygous
	for(j=0; j<snp_num; j++){
		if(!bmarker[j]){
			for(i=0; i<indi_num; i+=4) ifs.read(ch,1);
			continue;
		}
		for(i=0; i<indi_num;){
			ifs.read(ch,1);
			if(!ifs) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
			b=ch[0];
			k=0;
			while(k < 7 && i < indi_num){
				zvec(i) = b[k++] + b[k++];
				i++;
		   }
		}
		tgv += zvec.cast<double>() * alpha[j];
		if(j==last_true_snp) break;
	}
	ifs.clear();
	ifs.close();
	
	std::ofstream ofs;
	ofs.open(output_file);
	for(i=0; i<indi_num; ++i)
		ofs<<indi[i]<<','<<tgv[i]<<std::endl;
	ifs.close();
	
}

void cal_dgv(
	const std::string marker_estimate_file,
	const std::string binary_genotype_file,
	const std::string output_file)
{
	if(file_check(binary_genotype_file + ".fam") && file_check(binary_genotype_file + ".bim") && file_check(binary_genotype_file + ".bed"))
		computeTotalGeneticValuePLINK(marker_estimate_file, binary_genotype_file, output_file);
	else throw("\nError: bed/bim/fam fileset NOT found.\n");
}
