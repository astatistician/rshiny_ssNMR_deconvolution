# Automated ssNMR Spectral Deconvolution Shiny App

This repository contains an R Shiny application designed for the automated deconvolution of solid-state NMR (ssNMR) mixture spectra based on two template/reference spectra. The app provides an intuitive interface for spectral decomposition analysis and supports interactive plots for detailed visualization.

## Publication
This Shiny application is described in the following manuscript:

> Prostko et al., **R Shiny App for the Automated Deconvolution of NMR Spectra to Quantify the Solid-State Forms of Pharmaceutical Mixtures**, *Metabolites*, 2022, 12(12), 1248.  
> DOI: [10.3390/metabo12121248](https://doi.org/10.3390/metabo12121248)


## Dependencies
To run this application, the following R packages must be installed beforehand:

- **PepsNMR** (available from Bioconductor)
- **plotly**
- **shiny**
- **shinyBS**
- **shinyFeedback**
- **shinyhelper**
- **stringi**
- **tidyverse**

### Installation
Use the following R code to install the required packages:

```r
# Install CRAN packages
cran_packages <- c("plotly", "tidyverse", "stringi", "shiny", "shinyhelper", "shinyFeedback", "shinyBS")
install.packages(setdiff(cran_packages, installed.packages()[,"Package"]))

# Install PepsNMR from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!"PepsNMR" %in% installed.packages()[,"Package"]) {
    BiocManager::install("PepsNMR")
}
```

## Usage
Open up either **server.R** or **ui.R** and click *Run*. Alternatively, go to **[https://valkenborg-lab.shinyapps.io/ssNMRdeconvolution/](https://valkenborg-lab.shinyapps.io/ssNMRdeconvolution/)** and use the application on-line.
For more information on how to use the app, refer to the user manual located here [metabolites-1975964-supplementary](./metabolites-1975964-supplementary.pdf).