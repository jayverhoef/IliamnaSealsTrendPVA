
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.1.1-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kotzeb0912)](https://cran.r-project.org/package=kotzeb0912) [![packageversion](https://img.shields.io/badge/Package%20version-1.0-orange.svg?style=flat-square)](commits/master)

[![Last-changedate](https://img.shields.io/badge/last%20change-2018--02--13-yellowgreen.svg)](/commits/master)

[![DOI](https://zenodo.org/badge/90825285.svg)](https://zenodo.org/badge/latestdoi/90825285)

# IliamnaSealsTrendPVA
## An R package in support of publication, "A Bayesian Analysis of Abundance, Trend and Population Viability for Harbor Seals in Iliamna Lake, Alaska." 

#### Jay M. Ver Hoef<sup>a</sup>

#### <sup>a</sup>NOAA Fisheries (NMFS) Alaska Fisheries Science Center 

As a scientific work, and in keeping with common scientific practicies, I kindly request that you cite my research project and applicable publications if you use my work(s) or data in your publications or presentations. Additionally, I strongly encourage and welcome collaboration to promote use of these data in the proper context and scope.  The publication is currently submitted:

#### Peter L. Boveng, Jay M. Ver Hoef, David E. Withrow and Josh M. London. A Bayesian Analysis of Abundance, Trend and Population Viability for Harbor Seals in Iliamna Lake, Alaska. In press *Risk Analysis*.


Executive Summary
-----------------

 Harbor seals in Iliamna Lake, Alaska, are a small, isolated population, and one of only two freshwater populations of harbor seals in the world, yet little is known about their abundance or risk for extinction. Bayesian hierarchical models were used to estimate abundance and trend of this population. Observational models were developed from aerial survey and harvest data, and they included effects for time-of-year and time-of-day on survey counts. Underlying models of abundance and trend were based on a Leslie matrix model that used prior information on vital rates from the literature. We developed three scenarios for variability in the priors and used them as part of a sensitivity analysis. The models were fitted using Markov chain Monte Carlo methods. The population production rate implied by the vital rate estimates was about 5\% per year, very similar to the average annual harvest rate. After a period of growth in the 1980s, the population appears to be relatively stable at around 400 individuals. A population viability analysis assessing the risk of quasi-extinction, defined as any reduction to 50 animals or below in the next 100 years, ranged from 1\% to 3\%, depending on the prior scenario. Although this is moderately low risk, it does not include genetic or catastrophic environmental events, which may have occurred to the population in the past, so our results should be applied cautiously. 

Installation
------------

Installation of this R data package is done through the `devtools::install_github()` function or by downloading the [source package from the latest release](https://github.com/jayverhoef/KrigLinCaution).

```
library("devtools")
install_github("jayverhoef/IliamnaSealsTrendPVA")
```

Run R Scripts
-------------

The knitr document used to create the manuscript can be found here on your computer system:

```
system.file("doc/IliamnaSealsTrendPVA.Rnw", package = "IliamnaSealsTrendPVA")
```

which contains all of the R code embedded in the Latex manuscript.  Stripping out the R code with the "purl" command yields a pure R script, which can be found here:


```
system.file("doc/IliamnaSealsTrendPVA.R", package = "IliamnaSealsTrendPVA")
```

To run the whole script from within R use:

```
# set working directory to /doc path in package
setwd(system.file("doc", package = "IliamnaSealsTrendPVA"))
# make a list of files and directories
file.list <- list.files()
# copy files to R temporary directory (or change to a permanent one of your choice)
file.copy(file.list, tempdir(), recursive = TRUE)
# set working directory to R temporary directory (or whereever you copied stuff)
setwd(tempdir())
# Run the R script
source("IliamnaSealsTrendPVA.R")
```

A pure Latex document can be found here:

```
system.file("doc/IliamnaSealsTrendPVA.tex", package = "IliamnaSealsTrendPVA")
```

-------------
##### Disclaimer

<sub>This repository is a scientific product and is not official communication of the Alaska Fisheries Science Center, the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. AFSC-MML has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.</sub>
