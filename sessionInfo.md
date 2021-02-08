sessionInfo
================

The R/RStudio and R packages version used to generate all codes and
documentation

    ## Warning: package 'survival' was built under R version 3.6.3

    ## Loading required package: Hmisc

    ## Loading required package: lattice

    ## Loading required package: Formula

    ## Loading required package: ggplot2

    ## 
    ## Attaching package: 'Hmisc'

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

    ## Loading required package: SparseM

    ## 
    ## Attaching package: 'SparseM'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

    ## Loading required package: gsubfn

    ## Loading required package: proto

    ## Loading required package: RSQLite

    ## Loading required package: prodlim

    ## Loading required package: data.table

    ## riskRegression version 2019.01.29

    ## 
    ## Attaching package: 'riskRegression'

    ## The following objects are masked from 'package:pec':
    ## 
    ##     ipcw, selectCox

    ## Warning: package 'survAUC' was built under R version 3.6.3

    ## Warning: package 'timeROC' was built under R version 3.6.3

    ## Warning: package 'table1' was built under R version 3.6.3

    ## 
    ## Attaching package: 'table1'

    ## The following objects are masked from 'package:Hmisc':
    ## 
    ##     label, label<-, units

    ## The following objects are masked from 'package:base':
    ## 
    ##     units, units<-

    ## Warning: package 'kableExtra' was built under R version 3.6.3

    ## 
    ## Attaching package: 'boot'

    ## The following object is masked from 'package:lattice':
    ## 
    ##     melanoma

    ## The following object is masked from 'package:survival':
    ## 
    ##     aml

    ## -- Attaching packages ------------------------------------------------------------------------------------------------------------------------------ tidyverse 1.2.1 --

    ## v tibble  2.1.3     v purrr   0.3.2
    ## v tidyr   1.0.0     v dplyr   0.8.3
    ## v readr   1.3.1     v stringr 1.4.0
    ## v tibble  2.1.3     v forcats 0.4.0

    ## -- Conflicts --------------------------------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::between()    masks data.table::between()
    ## x dplyr::filter()     masks stats::filter()
    ## x dplyr::first()      masks data.table::first()
    ## x dplyr::group_rows() masks kableExtra::group_rows()
    ## x dplyr::lag()        masks stats::lag()
    ## x dplyr::last()       masks data.table::last()
    ## x dplyr::src()        masks Hmisc::src()
    ## x dplyr::summarize()  masks Hmisc::summarize()
    ## x purrr::transpose()  masks data.table::transpose()

    ## Warning: package 'rsample' was built under R version 3.6.3

    ## Warning: package 'webshot' was built under R version 3.6.3

    ## It seems that the version of `phantomjs` installed is greater than or equal to the requested version.To install the requested version or downgrade to another version, use `force = TRUE`.

``` r
sI<-sessionInfo()
print(sI)
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19041)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] splines   stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] webshot_0.5.2             rsample_0.0.6            
    ##  [3] forcats_0.4.0             stringr_1.4.0            
    ##  [5] dplyr_0.8.3               purrr_0.3.2              
    ##  [7] readr_1.3.1               tidyr_1.0.0              
    ##  [9] tibble_2.1.3              tidyverse_1.2.1          
    ## [11] boot_1.3-22               kableExtra_1.1.0         
    ## [13] table1_1.2                knitr_1.25               
    ## [15] plotrix_3.7-6             timeROC_0.4              
    ## [17] survivalROC_1.0.3         survAUC_1.0-5            
    ## [19] riskRegression_2019.01.29 data.table_1.12.4        
    ## [21] pec_2018.07.26            prodlim_2018.04.18       
    ## [23] sqldf_0.4-11              RSQLite_2.1.2            
    ## [25] gsubfn_0.7                proto_1.0.0              
    ## [27] rms_5.1-3.1               SparseM_1.77             
    ## [29] Hmisc_4.2-0               ggplot2_3.2.1            
    ## [31] Formula_1.2-3             lattice_0.20-38          
    ## [33] survival_3.2-7           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] TH.data_1.0-10      colorspace_1.4-1    htmlTable_1.13.2   
    ##  [4] base64enc_0.1-3     rstudioapi_0.10     listenv_0.8.0      
    ##  [7] furrr_0.1.0         MatrixModels_0.4-1  bit64_0.9-7        
    ## [10] mvtnorm_1.0-11      lubridate_1.7.4     xml2_1.2.2         
    ## [13] codetools_0.2-16    jsonlite_1.6        broom_0.5.2        
    ## [16] cluster_2.1.0       compiler_3.6.1      httr_1.4.1         
    ## [19] backports_1.1.5     assertthat_0.2.1    Matrix_1.2-17      
    ## [22] lazyeval_0.2.2      cli_1.1.0           acepack_1.4.1      
    ## [25] htmltools_0.4.0     quantreg_5.51       tools_3.6.1        
    ## [28] gtable_0.3.0        glue_1.3.1          Rcpp_1.0.2         
    ## [31] cellranger_1.1.0    vctrs_0.3.1         nlme_3.1-140       
    ## [34] iterators_1.0.12    xfun_0.20           ps_1.3.0           
    ## [37] globals_0.12.5      rvest_0.3.4         lifecycle_0.1.0    
    ## [40] future_1.17.0       polspline_1.1.16    MASS_7.3-51.4      
    ## [43] zoo_1.8-6           scales_1.0.0        hms_0.5.1          
    ## [46] parallel_3.6.1      sandwich_2.5-1      RColorBrewer_1.1-2 
    ## [49] yaml_2.2.0          memoise_1.1.0       gridExtra_2.3      
    ## [52] rpart_4.1-15        latticeExtra_0.6-28 stringi_1.4.3      
    ## [55] foreach_1.4.7       checkmate_1.9.4     lava_1.6.6         
    ## [58] chron_2.3-54        rlang_0.4.6         pkgconfig_2.0.3    
    ## [61] evaluate_0.14       htmlwidgets_1.5.3   cmprsk_2.2-9       
    ## [64] processx_3.4.1      bit_1.1-14          tidyselect_0.2.5   
    ## [67] magrittr_1.5        R6_2.4.0            generics_0.0.2     
    ## [70] multcomp_1.4-10     DBI_1.0.0           pillar_1.4.2       
    ## [73] haven_2.3.1         foreign_0.8-71      withr_2.1.2        
    ## [76] abind_1.4-5         nnet_7.3-12         modelr_0.1.5       
    ## [79] crayon_1.3.4        rmarkdown_2.6       timereg_1.9.4      
    ## [82] grid_3.6.1          readxl_1.3.1        callr_3.3.2        
    ## [85] blob_1.2.0          digest_0.6.21       numDeriv_2016.8-1.1
    ## [88] munsell_0.5.0       viridisLite_0.3.0   tcltk_3.6.1
