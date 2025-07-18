## Installation
SLEMM consists of one C++ program (`slemm`) and two Python programs (`slemm_gamma.py` and `slemm_gwa.py`).
- The C++ program is used for genomic predictions and GWAS step 1.
- The two Python programs are used for GWAS step 2.

### Binary executables
- Binary executables are available in the [Releases](https://github.com/jiang18/slemm/releases/latest).
- `*-x86_64-linux.zip` contains three standalone executables: `slemm`, `slemm_gamma`, and `slemm_gwa`.
- `slemm` is statically compiled with [Intel oneAPI Base & HPC Toolkit 2023.1](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) and Eigen 3.3.9 on Ubuntu 22.
- `slemm_gamma` and `slemm_gwa` are statically compiled on CentOS Linux 7.

> [!NOTE]
> `slemm` should be executable. Type `chmod +x slemm*` on the command line if needed.

### Compilation of `slemm`
1. Requirements
    - Linux x86_64
    - [Intel oneAPI Base & HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)
    - Eigen >= 3.3.7
2. Change the Eigen path in `Makefile`. 
3. On the command line type `make` while in the directory containing the `Makefile`.
4. This should produce the executable named `slemm`.

### Use of [`slemm_gamma.py`](../bin/slemm_gamma.py) and [`slemm_gwa.py`](../bin/slemm_gwa.py)
1. Dependencies
    - Linux packages: python3, python3-devel, gcc, and gcc-c++
    - Python packages: cython, numpy, scipy, and [pgenlib](https://github.com/chrchang/plink-ng/tree/master/2.0/Python)
2. `python3 slemm_gwa.py --help`

## :computer::desktop_computer:For Windows
Windows Subsystem for Linux (WSL) can be used to run SLEMM on a Windows machine.
1. [Install Linux on Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/install)
    - Ubuntu will be installed by default, and other Linux distributions are available. 
2. Run [pre-compiled SLEMM binaries](https://github.com/jiang18/slemm/releases/latest) on a WSL Linux distribution as on a Linux machine.
