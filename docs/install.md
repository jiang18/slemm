SLEMM consists of one C++ program (**slemm**) and two Python programs (**slemm_gamma.py** and **slemm_gwa.py**).

## Installation
### Compilation of slemm
1. Requirements
    - Linux x86_64
    - Intel C++ >= 19.0 (not tested on older versions)
    - [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)
    - Eigen >= 3.3.7
2. Change the Eigen path in `Makefile`. 
3. On the command line type `make` while in the directory containing the Makefile.
4. This should produce the executable called `slemm`.

### Dependencies of slemm_gamma.py and slemm_gwa.py
- Linux packages: python3, python3-devel, gcc, and gcc-c++
- Python packages: cython, numpy, scipy, and [pgenlib](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)

## Pre-compiled binaries
