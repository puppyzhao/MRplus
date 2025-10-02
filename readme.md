#MRplus
MRplus is an R package designed to enhance Mendelian Randomization (MR) analysis by providing intuitive visualization tools. It helps researchers evaluate the sensitivity of MR results to different instrumental variable (IV) selection thresholds, LD pruning parameters, and GWAS dataset choices. The package is particularly useful for assessing robustness and potential publication bias in previously published MR studies.

## Features
Visualize the impact of different IV selection thresholds (e.g., 5e-8, 5e-7, 5e-6, 1e-5) on MR results.
Explore the effect of LD parameters (r² and kb) on the number of IVs.
Compare MR outcomes across multiple thresholds.
Display multi-exposure and multi-outcome p-value distributions using bubble plots
Facilitate quick retrospective checks for robustness and potential bias.

## Installation
You can install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("puppyzhao/MRplus")

Dependencies and Setup
MRplus requires access to the OpenGWAS database via the ieuopengwas package.
To use the package, please follow these steps:
Register for an account at OpenGWAS.
Obtain your personal token after registration.
In R, set your token using:
Sys.setenv(OPENGWAS_JWT = "...your token...")
Without this setup, MRplus will not function properly.

Usage Example
library(MRplus)
# Usage example1:
mr_iv_count_across_ld_and_pval_threshold("ieu-b-4872")
# Usage example2:
mr_result_across_pval_thresholds("ieu-b-4872", "ukb-b-16890")
# Usage example3:
mr_result_matrix_across_datasets(
  exposure_ids = c("ieu-b-4872", "bbj-a-80"),
  outcome_ids = c("ukb-b-16890", "ebi-a-GCST90013972"),
  p_values = c(1e-5, 5e-6, 5e-7, 5e-8)
)

License
This package is released under the MIT License.

