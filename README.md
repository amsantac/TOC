TOC R package repository
======

To install the TOC R package from this account on Github, first run the following command line in R in order to install the devtools R package:

```{r}
install.packages("devtools")
```

Then load the devtools package and install the TOC R package from Github:

```{r}
library(devtools)
install_github("amsantac/TOC")
```

Now the TOC package can be loaded, as well as sample data from inside the package:

```{r}
library(TOC)
index <- raster(system.file("external/Prob_Map2.rst", package="TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package="TOC"))
mask <- raster(system.file("external/MASK3.rst", package="TOC"))
```

The TOC and ROC curves can be generated and plotted:

```{r}
tocd <- TOC(index, boolean, mask, NAval=0, nthres=100, uncertainty=TRUE, progress=TRUE)
tocd
plot.TOC(tocd)

rocd <- ROC(index, boolean, mask, NAval=0, nthres=100, uncertainty=TRUE, progress=TRUE)
rocd
plot.ROC(rocd)
```

For help on the functions implemented in the TOC package see:

```{r}
?TOC
```
