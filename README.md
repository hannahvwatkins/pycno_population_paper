## *Pycnopodia helianthoides* population decline manuscript

This repository contains the code, model outputs and figures for the time series and spatial analyses for *Pycnopodia helianthoides* in Canada for our upcoming mauscript. It is closely related to our recent [COSEWIC repository](https://github.com/hannahvwatkins/pycno_population_analysis), with updated figures and supplemental details.

*Code*

All `.R` files with the prefix `cleaning_` contain code to process the raw data from each source for use in the analyses. The `custom_functions.R` script contains custom functions designed for use with these analyses. These include those created specficially for handling ordinal response output by Dan Greenberg, and functions for validating Stan models in R by myself and Helen Yan - this must be sourced at the beginning of each R script to be able to access these functions. All code is heavily annotated to allow for greater readibility, but if you have any questions or concerns, please feel free to get in touch.

*Figures*

The figures folder contains all of the figures in the main text and appendix.


Session Info
```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8   
[3] LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C                   
[5] LC_TIME=English_Canada.utf8    

time zone: America/Vancouver
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.1.2     DHARMa_0.4.6        rstan_2.26.22       StanHeaders_2.26.27
 [5] cmdstanr_0.6.0      lubridate_1.9.2     forcats_1.0.0       stringr_1.5.0      
 [9] dplyr_1.1.2         purrr_1.0.1         readr_2.1.4         tidyr_1.3.0        
[13] tibble_3.2.1        ggplot2_3.4.2       tidyverse_2.0.0    

loaded via a namespace (and not attached):
 [1] gridExtra_2.3        inline_0.3.19        rlang_1.1.1          magrittr_2.0.3      
 [5] matrixStats_1.0.0    compiler_4.3.1       mgcv_1.8-42          loo_2.6.0           
 [9] systemfonts_1.0.4    callr_3.7.3          vctrs_0.6.3          pkgconfig_2.0.3     
[13] crayon_1.5.2         fastmap_1.1.1        backports_1.4.1      ellipsis_0.3.2      
[17] labeling_0.4.2       utf8_1.2.3           promises_1.2.0.1     rmarkdown_2.23      
[21] markdown_1.7         tzdb_0.4.0           ps_1.7.5             nloptr_2.0.3        
[25] ragg_1.2.5           bit_4.0.5            xfun_0.39            jsonlite_1.8.7      
[29] later_1.3.1          parallel_4.3.1       prettyunits_1.1.1    R6_2.5.1            
[33] gap.datasets_0.0.5   stringi_1.7.12       qgam_1.3.4           boot_1.3-28.1       
[37] Rcpp_1.0.11          iterators_1.0.14     knitr_1.43           httpuv_1.6.11       
[41] Matrix_1.6-0         splines_4.3.1        timechange_0.2.0     tidyselect_1.2.0    
[45] rstudioapi_0.15.0    abind_1.4-5          doParallel_1.0.17    ggtext_0.1.2        
[49] codetools_0.2-19     curl_5.0.1           processx_3.8.2       pkgbuild_1.4.2      
[53] lattice_0.21-8       plyr_1.8.8           shiny_1.7.4.1        withr_2.5.0         
[57] posterior_1.4.1      evaluate_0.21        RcppParallel_5.1.7   xml2_1.3.5          
[61] pillar_1.9.0         gap_1.5-1            tensorA_0.36.2       checkmate_2.2.0     
[65] foreach_1.5.2        stats4_4.3.1         distributional_0.3.2 generics_0.1.3      
[69] vroom_1.6.3          hms_1.1.3            commonmark_1.9.0     munsell_0.5.0       
[73] scales_1.2.1         minqa_1.2.5          xtable_1.8-4         glue_1.6.2          
[77] tools_4.3.1          lme4_1.1-34          grid_4.3.1           colorspace_2.1-0    
[81] nlme_3.1-162         cli_3.6.1            textshaping_0.3.6    fansi_1.0.4         
[85] V8_4.3.2             gtable_0.3.3         digest_0.6.33        farver_2.1.1        
[89] htmltools_0.5.5      lifecycle_1.0.3      mime_0.12            gridtext_0.1.5      
[93] bit64_4.0.5          MASS_7.3-60           
```

