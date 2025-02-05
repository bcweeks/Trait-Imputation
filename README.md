# Trait-Imputation
This is the pipeline used for imputing missing values in a large-scale trait dataset assembled using Skelevision.

Using phylogenetic models, we estimate values for traits that were not confidently measured by Skelevision. This approach is outlined in detail in:

Weeks, B.C., Z. Zhou, C.M. Probst, J.S. Berv, B. O'Brien, B.W. Benz, H.R. Skeen, M. Ziebell, L. Bodt, and D.F. Fouhey. 2024. Skeletal trait measurements for thousands of bird species. bioRxiv. doi: https://doi.org/10.1101/2024.12.19.629481

The trait imputation is done in R; users can implement the approach for the Weeks et al. (2024) dataset using the Trait Imputation_2-5-25.R script. The script requires functions found in the data_cleaning_functions.R script. 
