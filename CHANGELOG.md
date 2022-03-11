# Change Log

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
