# Visualize quality control data from methylation array dataset with methview.qc <img src="img/methview.qc_hexsticker.png" align="right" height="140" />  

![GitHub repo size](https://img.shields.io/github/repo-size/YoannPa/methview.qc)
![GitHub issues](https://img.shields.io/github/issues-raw/YoannPa/methview.qc)
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/YoannPa/methview.qc)  

_Methylation arrays have become a standard in clinics for epigenetic tumor profiling. While there is now a broad choice of tools to process and analyze methylation array data (RnBeads, Minfi, Sesame), none of them propose neat and complete quality control plots._  
_Here we introduce the R package **Methview.qc** which allows you to generate quality control plots from methylation array dataset._  
_**Methview.qc** utilizes RnBSet to store methylation datasets from both HumanMethylation450 and MethylationEPIC arrays. Its internal use of data.table boosts the memory-efficient processing of data, and the ggplot2 grammar provides R users with neat, yet customizable, publication-ready graphics to visualize quality control data._  

**Author: PAGEAUD Y.<sup>1</sup>**  
**Contributors: RATHGEBER A.<sup>1</sup>**  
**1-** [**DKFZ - Division of Applied Bioinformatics, Germany.**](https://www.dkfz.de/en/applied-bioinformatics/index.php)  
**How to cite:** _Pageaud Y. et al., Visualize quality control data from methylation array dataset with methview.qc_  

![GitHub R package version](https://img.shields.io/github/r-package/v/YoannPa/methview.qc?label=Package%20version&logo=RStudio&logoColor=white&style=for-the-badge)  
<img src="https://img.shields.io/static/v1?label=compatibility&message=4.3.0&color=blue&logo=R&logoColor=white&style=for-the-badge" />  
![GitHub last commit](https://img.shields.io/github/last-commit/YoannPa/methview.qc?logo=git&style=for-the-badge)  
![GitHub](https://img.shields.io/github/license/YoannPa/methview.qc?color=brightgreen&style=for-the-badge) 

## Content
Currently the package methview.qc contains **23 functions**:

* `devscore.fluo()` - Computes a deviation score between samples fluorescence and an internal HM450K or MethylationEPIC reference.  
* `devscore.heatmap()` - Plots QC deviation heatmaps based on samples fluorescence deviation score.  
* `get_expected_intensity()` - Provides expected intensity for a given methylation array probe ID.  
* `get_IDATs_runinfo()` - Retrieves runinfo data from a methylation array sample's IDAT files.  
* `get_platform()` - Detects platform used to generate data in the RnBSet.
* `gp_density_map()` - Displays the distribution of genotyping probes values in a cohort.
* `load_metharray_QC_meta()` - Loads methylation array QC metadata as a data.table.  
* `mergeQC_intensities_and_meta()` - Merges red and green channels intensities with QC probes metadata.  
* `plot_all_qc()` - Draws and saves all quality control plots available in methview.qc.  
* `plot_array_QCprobe()` - Plots fluorescence intensities barplots for a single QC methylation array probe.  
* `plot_array_QCtarget()` - Plots samples fluorescence intensities distribution for QC probes of a specific target type.  
* `plot_asso_all_annot()` - Draws association test results between all annotations from an RnBSet.  
* `plot_asso_annot_PC()` - Draws association test results between annotations from an RnBSet and PCs from a prcomp object.  
* `plot_asso_annot_QC()` - Draws association test results between annotations from an RnBSet and QC probes intensities.
* `RnB2PCA()`	- Computes a PCA from an RnBSet on a subset of selected probes.
* `rnb_add_runinfo()`	- Adds IDATs runinfo to an RnBSet pheno table.
* `rnb_biplot()` - Draws a customizable PCA biplots on samples methylation array QC data.  
* `rnb_crossbiplot()` - Draws a customizable PCA cross biplot on samples methylation array QC data.  
* `rnb_test_asso_all_annot()` - Tests associations between all annotations in a RnBSet.  
* `rnb_test_asso_annot_PC()` - Tests associations between annotations from an RnBSet and PCs from a prcomp object.  
* `rnb_test_asso_annot_QC()` - Tests associations between annotations from an RnBSet and QC probes intensities.  
* `snp_heatmap()` - Draws a heatmap from methylation array genotyping probes.  
* `target.biplot()` - Draws a PCA biplots on methylation array quality control targets.  

## Prerequisites
### Install Bioconductor dependencies
In R do:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(pkgs = c("RnBeads", "RnBeads.hg19", "minfiData", "minfiDataEPIC", "IlluminaDataTestFiles"))
```
### Install CRAN dependencies
In R do:
```R
install.packages(c('ggplot2', 'data.table', 'parallel', 'RColorBrewer', 'grDevices'))
```
### Install the [BiocompR](https://github.com/YoannPa/BiocompR) package
In R do:
```R
devtools::install_github("YoannPa/BiocompR")
```
## Installing Methview.qc
In R do:
```R
devtools::install_github("YoannPa/methview.qc")
```
## Problems ? / I need help !
For any questions related to bugs please check the section "**Known Issues**" available below.  
If the issue you experience is not adressed in the known issues please check if an existing issue adresses your point [here](https://github.com/YoannPa/methview.qc/issues/). If not, create a [new issue here](https://github.com/YoannPa/methview.qc/issues/new).

### Known Issues
**❎  Error: C++14 standard requested but CXX14 is not defined**  
In CentOS Linux release 7.4 some users experienced this error happening during the installation of methview.qc. The error message is raised during installation of the sparseMatrixStats package.  
To fix this issue proceed as following:
1. Go to your personnal Linux directory:  
```bash
cd ~/
```
2. Create directory **.R/** and the file **.R/Makevars**  
```bash
mkdir .R
nano .R/Makevars
```
3. In the opened file .R/Makevars paste the following line:  
```bash
CXX14 = g++ -std=c++1y -Wno-unused-variable -Wno-unused-function -fPIC
```
4. Save with **Ctrl + O** and close the file with **Ctrl + X**.   
5. Restart your R session and try again to install methview.qc:  
```R
devtools::install()
```
Many thanks to [**@Lena-Vo**](https://github.com/Lena-Vo) who provided this solution.  

**❎  Error: \`row.names<-.data.frame\`(\`*tmp*\`, value = value) : invalid 'row.names' length**  
```R
`row.names<-.data.frame`(`*tmp*`, value = value) : 
  invalid 'row.names' length
```
This error can arise from Bioconductor packages incompatibilities. Incriminated packages are dependencies of methview.qc.  
In order to fix this issue it is recommended to run the following command:  
```R
BiocManager::valid()
```
Running this command allows you to list Bioconductor packages that are either "out-of-date" or "too new" regarding your current installation.  
Among logs returned, you will be offered to "create a valid installation with
" a dedicated Biocmanager::install command.  
Using the command usually fix packages incompatibilities, and solve the issue.  
If you still have this error after reinstalling these packages, feel free to
create an issue.  

## Development & Feature requests
If you wish to contribute to the development of this package, or if you would like me to add a new feature that would be useful in Methview.qc, please write me at [y.pageaud@dkfz.de](y.pageaud@dkfz.de).

## References
1. [_Assenov, Y. et al. Comprehensive analysis of DNA methylation data with RnBeads. **Nature Methods 11, 1138–1140** (2014)._](https://www.nature.com/articles/nmeth.3115)  
2. [_Pageaud Y. et al., BiocompR - Advanced visualizations for data comparison._](https://github.com/YoannPa/BiocompR)  

## Licence
The repository methview.qc is currently under the GPL-3.0 licence.  
