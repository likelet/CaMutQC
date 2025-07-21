<img src="https://github.com/likelet/CaMutQC/blob/WX/vignettes/CaMutQC_logo.png" align="left" height="120" /></a>


## &emsp; Integrative quality control of cancer somatic mutations with  <br />&emsp; [Ca]()ncer [Mut]()ation [Q]()uality [C]()ontrol

<br />

## Introduction

The quality control of cancer somatic mutations has a great significance in the research of cancer genomics. It helps reduce the noise, systematic bias, as well as the false positive rate in dataset, and to pick out candidate mutations related to tumorigenesis. However, existing tools were not designed especially for cancer somatic mutations, and parameters and filtering standards used are different among them. Besides, there is no such a platform, database or tool stores or summarizes the standards for filtration of cancer somatic mutations applied in previous studies, which increases the time and energy spent on related research.  

Therefore, we present this R package, CaMutQC, for the comprehensive filtration and selection of cancer somatic mutations for tumor-normal paired samples. CaMutQC is able to filter false positive mutations generated due to technical issues, as well as to select candidate cancer mutations through a series of well-structured functions by labeling mutations with various flags.

Also, a detailed and vivid filter report will be generated after completing a whole filtration or selection section.

## Installation
1. Through Bioconductor (**recommended**) to get the most stable version:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# the latest version is on the devel branch, but might not be the most stable version
# BiocManager::install(version='devel')
# the most stable version is on the release branch (by default)
BiocManager::install("CaMutQC")
```

2. Through Github to get the latest version: 
```R
#  install the latest version from GitHub
if(!require(devtools)) install.packages("devtools")
devtools::install_github("likelet/CaMutQC")
```

## Usage
This is a simple workflow of CaMutQC.   
<div  align="left">   

<img src="https://github.com/likelet/CaMutQC/blob/WX/vignettes/CaMutQC-workflow.png" height="500" width="700" alt = "CaMutQC framework"/>

</div>

A detailed manual can be found [here](https://likelet.github.io/CaMutQC/).

## Shiny APP
We developed a Shiny application for CaMutQC to enhance its accessibility. Users can use the following code to launch Shiny app build with the package.
```R
pkg_required_shiny <- c('shiny','shinyjs', 'shinyFiles', 'DT')

# Install required packages
checkPackage <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package ", pkg, " needed for shiny app. Please install it."), call. = FALSE)
  }
}
invisible(lapply(pkg_required_shiny, checkPackage))
# run shiny app from shiny package
shiny::runApp(system.file("shiny", package = "CaMutQC"))
```

## Author

This software was mainly developed by:

* Xin Wang, sylviawang555@gmail.com

## Supervised by 

* [Jian Ren](renjian@sysucc.org.cn) from School of Life Sciences in Sun Yat-sen University 
* [Qi Zhao](zhaoqi@sysucc.org.cn) from Bioinformatic Center of Sun Yat-sen University Cancer Center 

## Maintainer
[Xin Wang](sylviawang555@gmail.com)  

## Citation
Xin Wang, Tengjia Jiang, Ao Shen, Yaru Chen, Yanqing Zhou, Jie Liu, Shuhan Zhao, Shifu Chen, Jian Ren, and Qi Zhao. CaMutQC: An R package for integrative quality control and filtration of cancer somatic mutations. Computational and Structural Biotechnology Journal. 2025;27:3147–3154. https://doi.org/10.1016/j.csbj.2025.07.011

