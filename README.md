[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![minimal R version](https://img.shields.io/badge/R%3E%3D-2.1.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kotzeb0912)](https://cran.r-project.org/) [![packageversion](https://img.shields.io/badge/Package%20version-0.0.0.9000-orange.svg?style=flat-square)](https://github.com/michaeldumelle/DumelleEtAl2021STLMM)
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--04--26-blue.svg)](https://github.com/michaeldumelle/DumelleEtAl2021STLMM)

## DumelleEtAl2021STLMM

### A Supplementary R Package to to "A Linear Mixed Model Formulation for Spatio-Temporal Random Processes with Computational Advances for the Product, Sum, and Product-Sum Covariance Functions"

##### Michael Dumelle<sup>1</sup>, Jay M. Ver Hoef<sup>2</sup>, Claudio Fuentes<sup>1</sup>, Alix Gitelman<sup>1</sup>

##### <sup>1</sup>Department of Statistics, Oregon State University, Corvallis, OR, USA
##### <sup>2</sup>NOAA Fisheries (NMFS) Alaska Fisheries Science Center, Marine Mammal Laboratory, Seattle, WA, USA

### Abstract
To properly characterize a spatio-temporal random process, it is necessary to understand the process' dependence structure. It is common to describe this dependence using a single random error having a complicated covariance. Instead of using the single random error approach, we describe spatio-temporal random processes using linear mixed models having several random errors; each random error describes a specific quality of the covariance. This linear mixed model formulation is general, intuitive, and contains many commonly used covariance functions as special cases. We focus on using the linear mixed model formulation to express three covariance functions: product (separable), sum (linear), and product-sum. We discuss benefits and drawbacks of each covariance function and propose novel algorithms using Stegle eigendecompositions, a recursive application of the Sherman-Morrison-Woodbury formula, and Helmert-Wolf blocking to efficiently invert their covariance matrices, even when every spatial location is not observed at every time point. Via a simulation study and an analysis of temperature data in Oregon, USA, we assess computational and model performance of these covariance functions when estimated using restricted maximum likelihood (likelihood-based) and Cressie's weighted least squares (semivariogram-based). We end by offering guidelines for choosing among combinations of the covariance functions and estimation methods based on properties of observed data and the desired balance between computational efficiency and model performance. 


### Package Overview

This supplementary R package contains all files used in creation of this document. 

### Installation

The easiest way to install this R package is to run
```
install.packages("devtools")
library(devtools)
devtools::install_github("michaeldumelle/DumelleEtAl2021STLMM@main")
```

After installation, the associated package files can be located on your machine via the `system.file()` function.

### Preprint

The preprint is available at
```
system.file("preprint/main.pdf", package = "DumelleEtAl2021STLMM")
```

or alternatively can be downloaded [here](https://github.com/michaeldumelle/DumelleEtAl2021STLMM/blob/main/inst/preprint/main.pdf).

Supplementary material for the preprint is available at 

```
system.file("preprint/supplementary.pdf", package = "DumelleEtAl2021STLMM")
```

or alternatively can be downloaded [here](https://github.com/michaeldumelle/DumelleEtAl2021STLMM/blob/main/inst/preprint/supplementary.pdf).

### Images

All images can be found at
```
system.file("images", package = "DumelleEtAl2021STLMM")
```

All R scripts used to create the images can be found at
```
system.file("scripts/images", package = "DumelleEtAl2021STLMM")
```

These files are named corresponding to the figure numbers in the preprint.

### Inversion and Empirical Semivariogram Computational Benchmarks

All R scripts used to study inversion and empirical semivariogram computational benchmarks (Section 4.3) can be found at
```
system.file("scripts/inverses", package = "DumelleEtAl2021STLMM")
```

All output from these R scripts can be found at
```
system.file("output/inverses", package = "DumelleEtAl2021STLMM")
```

### Simulation Study 

All R scripts used in the simulation study (Section 5) can be found at
```
system.file("scripts/simulations", package = "DumelleEtAl2021STLMM")
```

All output from these R scripts can be found at
```
system.file("output/simulations", package = "DumelleEtAl2021STLMM")
```

### Data Analysis 

The data can be loaded after installing the R package and running
```
library(DumelleEtAl2021STLMM)
data("or_data")
```

All R scripts used in the data analysis (Section 6) can be found at
```
system.file("scripts/dataanalysis", package = "DumelleEtAl2021STLMM")
```

All output from these R scripts can be found at
```
system.file("output/dataanalysis", package = "DumelleEtAl2021STLMM")
```

### Thank you for visiting -- we hope you enjoyed the work!
