TOC R package repository
======

The official release of the TOC R package can be found at [CRAN](http://cran.r-project.org/package=TOC). To install from CRAN enter the following command line in R (on Windows OS, run R as Administrator): 

```{r}
install.packages("TOC")
```

To install the development version of the TOC R package from this account on Github, first install the devtools R package using the following command:

```{r}
install.packages("devtools")
```

Then install the TOC package from Github using the devtools package:

```{r}
devtools::install_github("amsantac/TOC")
```

Once installed, the TOC package can be loaded. The package includes sample data to run the examples:

```{r}
library(TOC)
index <- raster(system.file("external/Prob_Map2.rst", package="TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package="TOC"))
mask <- raster(system.file("external/MASK3.rst", package="TOC"))
```

The TOC and ROC curves can be generated and plotted using the following instructions:

```{r}
tocd <- TOC(index, boolean, mask, nthres=100, progress=TRUE)
tocd
plot(tocd)

rocd <- ROC(index, boolean, mask, nthres=100, progress=TRUE)
rocd
plot(rocd)
```

For help on the functions implemented in the TOC package see the [reference manual](/TOC-manual.pdf) or run the following command lines in R:

```{r}
?'TOC-package'
?TOC
```
An interactive web application for TOC developed with [Shiny](http://shiny.rstudio.com) can be found [here](https://amsantac.shinyapps.io/TOCapp). The [Shiny TOC app](https://amsantac.shinyapps.io/TOCapp) allows to construct the TOC and ROC curves for spatial and non-spatial data. 