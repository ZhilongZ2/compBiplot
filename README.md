
<!-- README.md is generated from README.Rmd. Please edit that file -->

# compBiplot

<!-- badges: start -->
<!-- badges: end -->

The goal of compBiplot is to provide tools for constructing log-ratio
biplots for compositional datasets. The package implements centered and
isometric log-ratio transformations (clr and ilr), and generates
covariance- or correlation-based biplots via log-ratio PCA (LR-PCA).

## Installation

You can install the development version of **compBiplot** from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ZhilongZ2/compBiplot")
library(compBiplot)
```

## Example

Here are some basic examples which show you how to use **compBiplot**:

``` r
library(compBiplot)
## basic example code

# The package includes the cleaned Ecuador–Mexico blood panel:
data("ECU.MF")
head(ECU.MF)
#>   age neutrophils lymphocytes monocytes eosinophils basophils redbloodcells
#> 1  33        5400        3870       880         230        50           5.5
#> 2  35        3320        4790       610         200        30           5.8
#> 3  70        3280        1410       271          10        20           5.1
#> 4  53        5220        2380       670         130        70           5.2
#> 5  30        2620        1350       310          60        10           4.8
#> 6  27        2520        2170       500         270        20           5.3
#>    mcv  mch mchc rdwP hemoglobin hematocritP platelets  mpv Condition
#> 1 87.8 29.1 33.1 12.7       15.9        48.0       497  9.8         0
#> 2 89.3 31.0 34.7 12.6       18.0        51.9       342  8.9         0
#> 3 88.5 29.6 33.5 14.0       15.0        44.8       264  9.0         1
#> 4 92.3 31.9 34.6 12.4       16.5        47.7       394  9.4         0
#> 5 87.3 28.3 32.4 13.8       13.6        42.0       178 11.1         1
#> 6 86.0 28.2 32.7 14.7       14.9        45.0       402  9.7         0

# ?ECU.MF
# See ‘data-raw/ECU_MF.R’ for full preprocessing code.
```

**compBiplot** also includes two common log-ratio transformations that
deal with compositional data. They are implemented to be able to handle
zero entries:

``` r
# Extract the five compositional parts

X <- as.matrix(ECU.MF[, c("neutrophils", "lymphocytes",
"monocytes", "eosinophils", "basophils")])

clr_X <- clr_transform(X, zero.method = "cmultRepl")
head(clr_X)
#>   neutrophils lymphocytes  monocytes eosinophils basophils
#> 1    1.997117    1.663972 0.18288426  -1.1589583 -2.685015
#> 2    1.768724    2.135290 0.07446307  -1.0406785 -2.937799
#> 3    2.846123    2.001869 0.35264285  -2.9468909 -2.253744
#> 4    2.168570    1.383173 0.11559454  -1.5241487 -2.143188
#> 5    2.428471    1.765402 0.29411401  -1.3481137 -3.139873
#> 6    1.767362    1.617831 0.14995631  -0.4662298 -3.068920

ilr_X <- ilr_transform(X, zero.method = "cmultRepl")
head(ilr_X)
#>         ilr1      ilr2      ilr3      ilr4
#> 1 -0.2355687 -1.345309 -2.113347 -3.001938
#> 2  0.2592010 -1.533008 -2.049741 -3.284559
#> 3 -0.5969775 -1.691253 -4.053376 -2.519762
#> 4 -0.5553595 -1.355610 -2.378620 -2.396157
#> 5 -0.4688611 -1.471998 -2.463071 -3.510485
#> 6 -0.1057349 -1.259560 -1.424277 -3.431156
```

Remaining *TODOS* for compBiplot:  
\* Add R functions: corrbiplot, covbiplot  
\* Add tests and examples  
\* Add vignettes  
\* Refine documents
