## Installation
SLEMM consists of one C++ program (`slemm`) and two Python programs (`slemm_gamma.py` and `slemm_gwa.py`).
- The C++ program is used for genomic predictions and GWAS step 1.
- The two Python programs are used for GWAS step 2.

### Binary executables
- Binary executables are available in the [Releases](https://github.com/jiang18/slemm/releases/latest).
- `slemm` in `*-x86_64-linux-mkl.zip` is statically compiled with Intel Math Kernel Library 2020 on CentOS Linux 7.
- `slemm` in `*-x86_64-wsl-ubuntu-onemkl.zip` is statically compiled with [Intel oneAPI Base and HPC Toolkits 2023.1](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) on Ubuntu 22.
- Standalone executables of `slemm_gamma` and `slemm_gwa` are compiled on CentOS Linux 7.
> **Note**
> `slemm` should be executable. Type `chmod +x slemm*` on the command line if needed.
> .

### Compilation of `slemm`
1. Requirements
    - Linux x86_64
    - [Intel C++](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) >= 19.0 (not tested on older versions)
    - [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)
    - Eigen >= 3.3.7
2. Change the Eigen path in `Makefile`. 
3. On the command line type `make` while in the directory containing the `Makefile`.
4. This should produce the executable named `slemm`.

### Use of `slemm_gamma.py` and `slemm_gwa.py`
1. Dependencies
    - Linux packages: python3, python3-devel, gcc, and gcc-c++
    - Python packages: cython, numpy, scipy, and [pgenlib](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)
2. `python3 slemm_gwa.py --help`

## :computer::desktop_computer:For Windows
Windows Subsystem for Linux (WSL) can be used to run SLEMM on a Windows machine.
1. [Install Linux on Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/install)
    - Ubuntu will be installed by default, and other Linux distributions are available. 
2. Run [pre-compiled SLEMM binaries](https://github.com/jiang18/slemm/releases/latest) on a WSL Linux distribution as on a Linux machine.
    - Use the SLEMM binary that is statically compiled on WSL Ubuntu 22 with Intel oneMKL and Eigen 3.3.9.
