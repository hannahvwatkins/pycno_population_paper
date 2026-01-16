[![DOI](https://zenodo.org/badge/751104985.svg)](https://doi.org/10.5281/zenodo.18274840)
## *Pycnopodia helianthoides* population decline manuscript

This repository contains the code, model outputs and figures for the time series and spatial analyses for *Pycnopodia helianthoides* in Canada for our upcoming mauscript. It is closely related to our recent [COSEWIC repository](https://github.com/hannahvwatkins/pycno_population_analysis), with updated figures and supplemental details.

*Code*

All `.R` files with the prefix `cleaning_` contain code to process the raw data from each source for use in the analyses. The `custom_functions.R` script contains custom functions designed for use with these analyses. These include those created specficially for handling ordinal response output by Dan Greenberg, and functions for validating Stan models in R by myself and Helen Yan - this must be sourced at the beginning of each R script to be able to access these functions. All code is heavily annotated to allow for greater readibility, but if you have any questions or concerns, please feel free to get in touch.

*Figures*

The figures folder contains all of the figures in the main text and appendix.


Session Info
```
R version 4.3.3 (2024-02-29)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Vancouver
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] sf_1.0-15          janitor_2.2.0      loo_2.7.0          DHARMa_0.4.6      
 [5] rstan_2.32.6       StanHeaders_2.32.6 cmdstanr_0.7.1     lubridate_1.9.3   
 [9] forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2       
[13] readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1     
[17] tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] gtable_0.3.4         tensorA_0.36.2.1     xfun_0.42            QuickJSR_1.1.3      
 [5] processx_3.8.4       inline_0.3.19        lattice_0.22-5       tzdb_0.4.0          
 [9] vctrs_0.6.5          tools_4.3.3          ps_1.7.6             generics_0.1.3      
[13] stats4_4.3.3         curl_5.2.1           parallel_4.3.3       proxy_0.4-27        
[17] fansi_1.0.6          pkgconfig_2.0.3      KernSmooth_2.23-22   Matrix_1.6-5        
[21] checkmate_2.3.1      distributional_0.4.0 RcppParallel_5.1.7   lifecycle_1.0.4     
[25] compiler_4.3.3       munsell_0.5.0        codetools_0.2-19     snakecase_0.11.1    
[29] class_7.3-22         nloptr_2.0.3         pillar_1.9.0         MASS_7.3-60.0.1     
[33] classInt_0.4-10      boot_1.3-29          abind_1.4-5          nlme_3.1-164        
[37] posterior_1.5.0      tidyselect_1.2.1     stringi_1.8.3        splines_4.3.3       
[41] grid_4.3.3           colorspace_2.1-0     cli_3.6.2            magrittr_2.0.3      
[45] pkgbuild_1.4.3       utf8_1.2.4           e1071_1.7-14         withr_3.0.0         
[49] scales_1.3.0         backports_1.4.1      timechange_0.3.0     matrixStats_1.2.0   
[53] lme4_1.1-35.1        gridExtra_2.3        hms_1.1.3            knitr_1.45          
[57] V8_4.4.2             rlang_1.1.3          Rcpp_1.0.12          DBI_1.2.2           
[61] glue_1.7.0           minqa_1.2.6          rstudioapi_0.15.0    jsonlite_1.8.8      
[65] R6_2.5.1             units_0.8-5        
```

