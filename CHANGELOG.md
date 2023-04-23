# Change Log

## Jan 31, 2023 (v0.89.5)
- Changed over from ctime to chrono for timing. 
- Fixed a bug that caused a segmentation fault when using the binary executables on recent Linux versions.

## Oct 30, 2022 (v0.89.3)
- Auto-disabled --iter_weighting if --window_size is larger than the total number of SNPs.

## Aug 22, 2022
- Fixed a minor issue in genotype reading.

## Jul 31, 2022
- Improved the speed of association tests in slemm_gamma and slemm_gwa. 
- Added a header line in output of slemm_gamma and slemm_gwa.

## Jul 02, 2022
- Improved the setting of --num_random_probes in --reml. 

## Jun 10, 2022
- Improved the computation for small samples. 

## May 17, 2022
- Added the empirical BLUE of fixed effects.
- Added two options (--sig_chisq and --indep_r2) to specify how to identify significant SNPs and fit them as fixed effects in the 2nd REML of --iter_weighting.
- Changed the option of --iter_weighting to a flag. 

## Apr 09, 2022
- Improved the Lanczos recursion.

## Mar 16, 2022
- Removed the flag of --reweight.
- Added an option (--iter_weighting) for specifying the number of iterations in iterative SNP weighting. 

## Mar 10, 2022
- SLEMM replaces [SSGP](https://github.com/jiang18/ssgp).

## Jan 22, 2022
- Added --error_weight_name for --lmm/--reml/--wls to weight individual error variances.

## Dec 13, 2021
- Added a flag (--lrt) to enable/disable the computation of logLL in --reml. Without --lrt, --reml time is reduced by ~50%.

## Dec 09, 2021
- Added a weighted least squares routine, --wls.

## Dec 05, 2021
- Reduced memory usage to 30-40% at a cost of 20-50% more time.
- Removed the variational Bayes routine and cleaned the code.

## Nov 15, 2021
- Improved genomic prediction.

## Oct 25, 2021
- Improved genomic prediction.
- Fixed a bug in --reml.

## Sep 06, 2021
- Fixed a bug in reading CSV files.
- Modified QR decomposition in --reml/--lmm.

## Jul 04, 2021
- Initial release
