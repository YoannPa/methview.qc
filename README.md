# Visualize quality control data from methylation array dataset with methview.qc <img src="img/methview.qc_hexsticker.png" align="right" height="140" />  
_**methview.qc** allows you to generate quality control plots from your methylation array dataset._  

**Author: PAGEAUD Y.<sup>1</sup>**  
**Contributors: RATHGEBER A.<sup>1</sup>**  
**1-** [**DKFZ - Division of Applied Bioinformatics, Germany.**](https://www.dkfz.de/en/applied-bioinformatics/index.php)  

**Version: 0.0.33 (Beta)**  
**R Compatibility: Version 4.0.5**  
**Last Update: 24/08/2022**  
**How to cite:** _Pageaud Y. et al., Visualize quality control data from methylation array dataset with methview.qc_  

## Content
Currently the package methview.qc contains **13 functions**:

* `cohort.gp.density()` - Displays the distribution of genotyping probes values in a cohort.
* `devscore.fluo()` - Computes a deviation score between samples fluorescence and an internal HM450K reference.  
* `devscore.heatmap` - Plots QC deviation heatmaps based on samples fluorescence deviation score.  
* `get_expected_intensity()` - Provides expected intensity for a given methylation array probe ID.  
* `get_platform()` - Detects platform used to generate data in the RnBSet.  
* `load_metharray_QC_meta()` - Loads methylation array QC metadata as a data.table.  
* `mergeQC_intensities_and_meta()` - Merges red and green channels intensities with QC probes metadata.  
* `plot_all_qc()` - Draws and saves all quality control plots available in methview.qc  
* `plot_array_QCprobe()` - Plots fluorescence intensities barplots for a single QC methylation array probe.  
* `plot_array_QCtarget()` - Plots samples fluorescence intensities distribution for QC probes of a specific target type.  
* `sampleQC.biplot()` - Draws a customizable PCA biplots on samples methylation array QC data.  
* `snp_heatmap()` - Draws a heatmap from methylation array genotyping probes.  
* `target.biplot()` - Draws a PCA biplots on methylation array quality control targets.  

## Prerequisites
### Install Bioconductor dependencies
In R do:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RnBeads")
```
### Install CRAN dependencies
```R
inst.pkgs = c('ggplot2', 'data.table', 'parallel', 'RColorBrewer', 'grDevices')
install.packages(inst.pkgs)
```

### Install the [BiocompR](https://github.com/YoannPa/BiocompR) package
The git repository of BiocompR is currently private.  
If you want to get access to it, please send a message at [**y.pageaud@dkfz.de**](y.pageaud@dkfz.de) with your Github username. Once you have accepted the invitation, the link in this section's title link should work.
1. In the Git repository click on "Clone or Download".
2. Copy the HTTPS link.
3. Open a terminal and type:
```bash
git clone https://github.com/YoannPa/BiocompR.git
```
4. Open the folder BiocompR/ and open the "BiocompR.Rproj" file in RStudio.
5. In the RStudio console, type:
```R
devtools::install()
```
If no error is displayed BiocompR should be installed now.  

## Installing
1. In the Git repository click on "Clone or Download".
2. Copy the HTTPS link.
3. Open a terminal and type:
```bash
git clone https://github.com/YoannPa/methview.qc.git
```
4. Open the folder methview.qc and open the "methview.qc.Rproj" file in RStudio.
5. In the RStudio console, type:
```R
devtools::install()
```

## Problems ? / I need help !
For any questions **Not related to bugs or development** please check the section "**Known Issues**" available below. If the issue you experience is not adressed in the known issues you can write me at [y.pageaud@dkfz.de](y.pageaud@dkfz.de).  

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

## Technical questions / Development / Feature request
If you encounters issues or if a feature you would expect is not available in a methview.qc function, please check if an existing issue adresses your point [here](https://github.com/YoannPa/methview.qc/issues/). If not, create a [new issue here](https://github.com/YoannPa/methview.qc/issues/new).  

## References
1. [_Assenov, Y. et al. Comprehensive analysis of DNA methylation data with RnBeads. **Nature Methods 11, 1138–1140** (2014)._](https://www.nature.com/articles/nmeth.3115)  
2. _Pageaud Y. et al., BiocompR - Advanced visualizations for data comparison._  

## Licence
The repository methview.qc is currently under the GPL-3.0 licence.  


