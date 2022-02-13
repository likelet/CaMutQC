# [Ca]()ncer [Mut]()ation [Q]()uality [C]()ontrol

The quality control of cancer somatic mutations has a great significance in the research of cancer genomics. It helps reduce the noise, systematic bias, as well as the false positive rate in dataset, and to pick out candidate mutations related to tumorigenesis. However, existing tools were not designed especially for cancer somatic mutations, and parameters and filtering standards used are different among them. Besides, there is no such a platform, database or tool stores or summarizes the standards for filtration of cancer somatic mutations applied in previous studies, which increases the time and energy spent on related research.  

Therefore, we present this R package, CaMutQC, for the comprehensive filtration and selection of cancer somatic mutations. CaMutQC is able to filter false positive mutations generated due to technical issues, as well as to select candidate cancer mutations through a series of well-structured functions by labeling mutations with various flags.

Also, a detailed and vivid filter report will be offered after completing a whole filtration or selection section.

## Installation

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

A detailed manual can be found [here](https://seqworld.com/CaMutQC/).

## Author

This software was mainly developed by:

* Xin Wang, sylviawang555@gmail.com

## Supervised by 

* [Jian Ren](renjian@sysucc.org.cn) from School of Life Sciences in Sun Yat-sen University 
* [Qi Zhao](zhaoqi@sysucc.org.cn) from Bioinformatic Center of Sun Yat-sen University Cancer Center 

## Maintainer
[Xin Wang](sylviawang555@gmail.com), Sun Yat-sen university  <br/>
