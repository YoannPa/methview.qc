% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_fun.R
\name{prep_annot_asso_fromRnB}
\alias{prep_annot_asso_fromRnB}
\title{Prepares annotations to be tested for associations.}
\usage{
prep_annot_asso_fromRnB(RnBSet, verbose = FALSE)
}
\arguments{
\item{RnBSet}{An \code{RnBSet} basic object for storing methylation array
DNA methylation and experimental quality information
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

\item{verbose}{A \code{logical} to display information about the step-by-step
processing of the data if TRUE (Default: verbose = FALSE).}
}
\value{
A \code{list} containing updated annotations, the number of
        annotations available, and the RnBSet annotation table it contains.
}
\description{
Prepares annotations to be tested for associations.
}
\examples{
# Create an RnBSet for Human Methylation 450K data
require(Biobase)
idat.dir <- system.file("extdata", package = "minfiData")
sample.annotation <- system.file(
    "extdata", "SampleSheet.csv", package = "minfiData")
data.source <- c(idat.dir, sample.annotation)
rnb.set <- RnBeads::rnb.execute.import(
    data.source = data.source, data.type = "idat.dir")
# Prepare annotations before association tests
methview.qc:::prep_annot_asso_fromRnB(RnBSet = rnb.set)
}
\author{
Yoann Pageaud.
}
\keyword{internal}
