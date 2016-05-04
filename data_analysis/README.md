
Virome evaluation
=============================


Data from Knowles et al. 2016 is re-evaluated using R (version 3.2.3) and packages MASS (version 7.3-44) and resample (version 0.4).

## Data

To perform the analysis the data used by **Extended Figure 4a** in Knowles et al. needs to be aquired from [doi:10.1038/nature17193](http://dx.doi.org/10.1038/nature17193) and placed in this folder named as:
 * *nature17193-sf4.csv*  

## Analysis

 * *DataReanalysis.R*   &nbsp; - &nbsp;   code to recreate Extended Figure 4a; and to calculate 95% confidence intervals for several correlation and robust regression techniques.

To perform the analysis (once the data has been aquired and placed in this folder) start R in this folder and run:

```r
	source("DataReanalysis.R")
```



###References


Knowles B. et al. 2016. Lytic to temperate switching of viral communities. *Nature*. 531: 466–470. [doi:10.1038/nature17193](http://dx.doi.org/10.1038/nature17193)

R Core Team. 2015. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.  URL https://www.R-project.org/.

Venables W. N., Ripley B. D. 2002. Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0

Hesterberg T. 2015. resample: Resampling Functions. R package version 0.4. https://CRAN.R-project.org/package=resample

