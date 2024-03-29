% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_fun.R
\name{mergeQC_intensities_and_meta}
\alias{mergeQC_intensities_and_meta}
\title{Merges red and green channels intensities with QC probes metadata.}
\usage{
mergeQC_intensities_and_meta(RnBSet, DT.QC.meta)
}
\arguments{
\item{RnBSet}{An \code{RnBSet} basic object for storing methylation
array DNA methylation and experimental quality information
(Bisulfite data not supported).
\itemize{
 \item{For more information about RnBSet object read
 \link[RnBeads]{RnBSet-class}.}
 \item{To create an RnBSet object run
 \link[RnBeads]{rnb.execute.import}.}
 \item{For additionnal options to import methylation array
 data in the RnBSet see options available in
 \link[RnBeads]{rnb.options}.}
}}

\item{DT.QC.meta}{A \code{data.table} with methylation array quality control
metadata obtained using \link{load_metharray_QC_meta}.}
}
\value{
A \code{data.table} list matching QC metadata with green channel and
        red channel intensities.
}
\description{
Merges red and green channels intensities with QC probes metadata.
}
\examples{
#Create an RnBSet for MethylationEPIC data
require(Biobase)
idat.dir <- system.file("extdata", package = "minfiDataEPIC")
sample.annotation <- system.file(
    "extdata", "SampleSheet.csv", package = "minfiDataEPIC")
data.source <- c(idat.dir, sample.annotation)
rnb.set <- RnBeads::rnb.execute.import(
    data.source = data.source, data.type = "idat.dir")
# Specify the column containing samples identifiers (sentrix barcodes).
RnBeads::rnb.options(identifiers.column = "barcode")
#Create the data.table with quality control metadata
dt.meta <- load_metharray_QC_meta(array.meta = "controlsEPIC")
# Merge red and green channels intensities with QC metadata
dt.mrg <- mergeQC_intensities_and_meta(
   RnBSet = rnb.set, DT.QC.meta = dt.meta)
dt.mrg
}
\author{
Yoann Pageaud.
}
