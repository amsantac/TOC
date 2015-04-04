TOC R package repository
======

To install the TOC R package from this account on Github, first install the devtools R package using the following command line in R (on Windows OS, run R as Administrator):

```{r}
install.packages("devtools")
```

Then load the devtools package and install the TOC R package from Github:

```{r}
library(devtools)
install_github("amsantac/TOC")
```

Now the TOC package can be loaded, as well as sample data included in the package:

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

For help on the functions implemented in the TOC package see the [PDF reference manual](https://github.com/amsantac/TOC/blob/master/TOC-manual.pdf?raw=true) or run the following command lines in R:

```{r}
?'TOC-package'
?TOC
```
