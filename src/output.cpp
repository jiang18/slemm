#include "slemm.h"
#include "chisq.h"

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
	const float& vg,
	const float& ve,
	const float& lrt) 
{
	unsigned i=0, j=0, k=0;
	std::ofstream ofs;
	// residuals
	ofs.open(output_file_prefix + ".reml.py.txt");
	ofs<<"IID Py";
	for(i=0; i<covar_names.size(); ++i) 
		ofs<<' '<<covar_names[i];
	ofs<<" R\n";
	for(i=0; i<indi_keep.size(); ++i) {
		ofs<<indi_keep[i]<<' '<<Py(i)<<' ';
		ofs<<xmat.row(i)<<' '<<rvec(i)<<"\n";
	}
	ofs.close();
	
	// snp effects
	k=0;
	ofs.open(output_file_prefix + ".reml.snp.csv");
	ofs<<"SNP,Chr,Pos,Allele1,Allele2,MAF,HWE_Pval,Group,Weight,Effect"<<std::endl;
	for(i=0; i<group_keep.size(); ++i) {
		for(j=0; j<ncols[i]; ++j) {
			std::string this_marker = marker_keep[k];
			ofs<<this_marker<<','<<marker2pos[this_marker].first<<','<<marker2pos[this_marker].second<<',';
			ofs<<marker2alleles[this_marker].first<<','<<marker2alleles[this_marker].second<<','<<marker2maf[this_marker]<<',';
			ofs<<marker2hwep[this_marker]<<','<<group_keep[i]<<','<<1<<',';
			ofs<<snp_blup(k)<<std::endl;
			++k;
		}
	}
	ofs.close();
	
	// blue of fixed effects
	ofs.open(output_file_prefix + ".reml.blue.csv");
	ofs<<"covar,blue\n";
	for(i=0; i<covar_names.size(); ++i) 
		ofs<<covar_names[i]<<','<<blue(i)<<"\n";
	ofs.close();
	
	// variance component estimates
	ofs.open(output_file_prefix + ".reml.vc.csv");
	ofs<<"pop_size,"<<indi_keep.size()<<std::endl;
	ofs<<"num_markers,"<<ncols.sum()<<std::endl;
	ofs<<"Vg,"<<vg<<std::endl;
	ofs<<"Ve,"<<ve<<std::endl;
	ofs<<"h2,"<<vg/(vg+ve)<<std::endl;
	ofs<<"LRT,"<<lrt<<std::endl;
	ofs<<"pval,"<<getOneDfChisqPval(lrt)<<std::endl;
	ofs.close();
}

void write_gstat_into_files(
	const std::string output_file_prefix,
	const int pop_size,
	const int window_size,
	const float vg,
	const std::map<std::string, float>& marker2maf,
	const std::map<std::string, float>& marker2hwep,
	std::map<std::string, std::pair<std::string, int> >& marker2pos,
	const std::vector<std::string>& marker_keep,
	const std::vector<gstat>& lmm_out ) 
{
	std::ofstream ofs;
	
	ofs.open(output_file_prefix + ".gstat.txt");
	ofs<<"variant\tChr\tPos\tMAC\tMAF\tHWEP\txPy\txPx\txBx\txSx\tWindowSize\tChiSq\tPval"<<std::endl;
	for(int j=0; j<lmm_out.size(); ++j) {
		if(lmm_out[j].vi < 0) {
			ofs<<"fake"+std::to_string(lmm_out[j].vi)<<"\tNA\tNA\tNA\tNA\tNA\t";
			ofs<<lmm_out[j].xPy<<"\t"<<lmm_out[j].xPx<<"\tNA\t"<<lmm_out[j].xSx<<"\t";
			
			float chisq = (lmm_out[j].xPy * lmm_out[j].xPy)/(vg*lmm_out[j].xPx);
			ofs<<"NA\t"<<chisq<<"\t"<<getOneDfChisqPval(chisq)<<std::endl;
		} else {
			std::string varid = marker_keep[lmm_out[j].vi];
			float maf = marker2maf.at(varid);
			
			ofs<<varid<<"\t"<<marker2pos[varid].first<<"\t"<<marker2pos[varid].second<<"\t";
			ofs<<int(maf*2*pop_size+0.5)<<"\t"<<maf<<"\t"<<marker2hwep.at(varid)<<"\t";
			ofs<<lmm_out[j].xPy<<"\t"<<lmm_out[j].xPx<<"\t"<<lmm_out[j].xBx<<"\t"<<lmm_out[j].xSx<<"\t";
			
			float chisq = (lmm_out[j].xPy * lmm_out[j].xPy)/(vg*lmm_out[j].xPx);
			ofs<<window_size<<"\t"<<chisq<<"\t"<<getOneDfChisqPval(chisq)<<std::endl;
		}
	}
	ofs.close();
}
