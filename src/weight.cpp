#include "slemm.h"

void cal_maf_from_genobits(const std::vector<std::vector<bool> >& geno, Ref<VectorXf> maf)
{
	const int n = geno[0].size();
	#pragma omp parallel for
	for(unsigned i=0; i<geno.size(); ++i) {
		int nz_bits = 0;
		for(unsigned j=0; j<n; ++j) {
			if(geno[i][j]) nz_bits ++;
		}
		maf(i) = nz_bits/float(n);
	}
	maf = ((maf.array() > 0.5).cast<float>() - maf.array()).abs();
}

float beta_pdf (const float a, const float b, const float rval )
{
	float value;

	if ( rval <= 0.0 || 1.0 <= rval )
	{
		value = 0.0;
	}
	else
	{
		value = pow ( rval, a - 1.0 ) * pow ( 1.0 - rval, b - 1.0 );
	}
	return value;
}

void update_weight_by_maf (
	const std::pair<float, float>& beta_pdf_params,
	const std::map<std::string, float>& marker2maf,
	std::map<std::string, std::pair<std::string, float> >& marker2group_weight )
{
	std::cout<<"\nWeighting markers for heritability using Beta PDF:\n";
	printf("Weight ~ X^(%.3f-1)*(1-X)^(%.3f-1), where X=MAF.\n\n", beta_pdf_params.first, beta_pdf_params.second);
	for(auto it=marker2maf.begin(); it!=marker2maf.end(); ++it)
		marker2group_weight[it->first].second = beta_pdf(beta_pdf_params.first, beta_pdf_params.second, it->second);
}

