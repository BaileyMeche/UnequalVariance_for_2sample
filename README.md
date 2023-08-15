# UnequalVariance_for_two_sample

## Partially Overlapping Data Power Simulation

This R notebook runs simulations to evaluate the statistical power of mean difference statistical tests for mixed paired and two sample data. 

The code allows specifying the sample sizes for the paired data (n), unpaired x data (n1), and unpaired y data (n2). It also allows specifying the effect size (delta), distribution shape (gamma), and correlation (rho) to use in simulating the data.

The key steps are:

1. Set input parameters to generate overlapping samples data 
2. Calculate test statistics and p-values for:
   - Bhoj et al. Zb and Tb
   - Magel and Fu M
   - Dubnicka R, Rw
   - Lin & Stivers Zls
   - Derrick et al. Tnew1, Tnew2
   - The paper's proposed T_adj  
   - Paired t-test
   - Signed rank test
3. Repeat simulation multiple times
4. Calculate power as the proportion of p-values below alpha across simulations
5. Output power results for different parameter combinations to Excel

The code is structured into sections:

- Preliminaries: Install and load R packages
- Generate data functions 
- Statistics functions to calculate each test
- Initialization of simulation parameters
- Outer loop over parameters
- Execute simulation and export results

The code allows easily modifying the parameters and tests evaluated. It outputs a tidy summary of power results across different scenarios that can be further analyzed in an Excel file. 

*Warning*: If running this code on a Jupyter Notebook, the final executing loop may need to be split over multiple code blocks to avoid storage issues. 
