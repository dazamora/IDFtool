




<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.com/dazamora/IDFtool.svg?branch=master)](https://travis-ci.com/dazamora/IDFtool) [![codecov](https://codecov.io/gh/dazamora/IDFtool/branch/master/graph/badge.svg)](https://codecov.io/gh/dazamora/IDFtool)

IDFtool
=======

IDFtool computes intensity-duration-frequency curves per specific time duration and different return periods. The intensity-duration-frequency curves are used in hydrology to express in a synthetic way, fixed a return period (T) and a duration (d) of a rainfall event. IDFtool included an uncertainty analysis in PDFs and IDF curves, by bootstrap method.

Instalation
-----------

Currently, you can install the version under development from [Github](https://github.com/dazamora/IDFtool), using these commands:

``` r
install.packages("devtools")
devtools::install_github("dazamora/IDFtool")
```

Example
-------

Meteorology station in the Farfan Airport in Tulua, Colombia.

``` r
library(IDFtool)
data(inten)
Test.idftool <- IDFCurve(Data = inten, Station='2610516', Duration = FALSE,
                         Periods = FALSE, Type = "gumbel", M.fit = "lmoments",
                         Plot = 1234, Strategy = 1, logaxe = "", CI = FALSE, 
                         CIpdf = TRUE, iter = 50, goodtest = FALSE,
                         Resolution = 300, SAVE = FALSE, name = TRUE)
#> [1] "Just compute a strategy"
```

<img src="README-unnamed-chunk-2-1.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-2.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-3.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-4.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-5.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-6.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-7.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-8.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-9.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-10.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-11.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-12.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-13.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-14.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-15.png" width="100%" style="display: block; margin: auto;" /><img src="README-unnamed-chunk-2-16.png" width="100%" style="display: block; margin: auto;" />
