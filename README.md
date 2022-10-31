## SLEMM (Stochastic-Lanczos-Expedited Mixed Models)
SLEMM is a software tool for large-scale genomic predictions and genome-wide association studies.

## Author and contact
[Jicai Jiang](https://cals.ncsu.edu/animal-science/people/jicai-jiang)

## Compilation
#### Requirements
1. x86-64 Linux
2. [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) 2019 or above
3. Eigen 3.3.9
#### Linux
1. Kernel version >= 3.10.0 (not tested on older versions)
2. Intel C++ >= 19.0 (not tested on older versions)

## Executable
https://github.com/jiang18/slemm/releases/tag/20220822

---

## Dependencies of slemm_gamma.py and slemm_gwa.py
- Linux packages: python3, python3-devel, gcc, and gcc-c++
- Python packages: cython, numpy, scipy, and [pgenlib](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)

It is usually easy to fullfill the dependencies if your machine has python3 and pip. If it is difficult for you to pip the python3 packages, please try the ones I have compiled on CentOS Linux 7 with GCC 4.8.5 (https://github.com/jiang18/slemm#executable). 

## Error weight
The argument ```--error_weight_name``` can be added in --reml/--wls/--lmm to weight individual error variance. To enable this, add one additional column in the phenotype CSV file and specify the header name by --error_weight_name. Values are typically 1/r^2-1, where r^2 is the reliability of deregressed EBVs. 
