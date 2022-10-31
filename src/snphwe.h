#ifndef SNPHWE_H_
#define SNPHWE_H_
#include <cstdint>

double SNPHWE2(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp);
int32_t SNPHWE_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);
int32_t SNPHWE_midp_t(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double thresh);
double hwep_chisq (int obs_hets, int obs_hom1, int obs_hom2);

#endif
