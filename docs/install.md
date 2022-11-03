## Installation
SLEMM consists of one C++ program (`slemm`) and two Python programs (`slemm_gamma.py` and `slemm_gwa.py`).
- The C++ program is used for genomic predictions and GWAS step 1.
- The two Python programs are used for GWAS step 2.

### Compilation of `slemm`
1. Requirements
    - Linux x86_64
    - [Intel C++](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) >= 19.0 (not tested on older versions)
    - [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)
    - Eigen >= 3.3.7
2. Change the Eigen path in `Makefile`. 
3. On the command line type `make` while in the directory containing the `Makefile`.
4. This should produce the executable called `slemm`.

### Use of `slemm_gamma.py` and `slemm_gwa.py`
1. Dependencies
    - Linux packages: python3, python3-devel, gcc, and gcc-c++
    - Python packages: cython, numpy, scipy, and [pgenlib](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)
2. `python3 slemm_gwa.py --help`

### Binary executables
- Binary executables are available in the [Releases](https://github.com/jiang18/slemm/releases).
- `slemm` is statically compiled with the Intel Math Kernel Library on CentOS Linux 7.
- `slemm_gamma` and `slemm_gwa` are statically compiled on CentOS Linux 7.
