
#' Detects platform used to generate data in the RnBSet.
#' 
#' @param RnBSet An \code{RnBSet} basic object for storing methylation array DNA
#'               methylation and experimental quality information (Bisulfite
#'               data not supported).
#'               \itemize{
#'                \item{For more information about RnBSet object read
#'                \link[RnBeads]{RnBSet-class}.}
#'                \item{To create an RnBSet object run
#'                \link[RnBeads]{rnb.execute.import}.}
#'                \item{For additionnal options to import methylation array data
#'                in the RnBSet see options available in
#'                \link[RnBeads]{rnb.options}.}
#'               }
#' @return A \code{character} string:
#'         \itemize{
#'          \item{"HM450K" if the RnBSet contains Human Methylation 450K data}
#'          \item{"MethylationEPIC" if the RnBSet contains MethylationEPIC data}
#'         }
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' library(RnBeads)
#' idat.dir <- "~/data/my_idat_dir/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' #Check platform used to generate the dataset.
#' get.platform(my.rnbset)

get.platform <- function(RnBSet){
  #Get the 65 or 59 genotyping (rs) probes
  rs.probes <- rownames(RnBSet@sites)[
    grepl(pattern = "rs", x = rownames(RnBSet@sites))]
  
  if(length(rs.probes) == 59){
    array.type <- "MethylationEPIC"
  } else if(length(rs.probes) == 65){
    array.type <- "HM450K"
  } else {
    stop(paste(
      "RnBSet platform not supported.",
      "Supported platforms are HM450K and MethylationEPIC.",
      "Please contact developper to request support for your methylation data.")
    )
  }
  return(array.type)
}

#' Loads methylation array QC metadata as a data.table.
#' @param array.meta A \code{character} string specifying the array type to load
#'                   quality control metadata from
#'                   (Supported: array.meta = c("controls450", "controlsEPIC")).
#' @return A \code{data.table} with HM450K or MethylationEPIC quality control
#'         metadata depending on the QC metadata requested.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Load quality control metadata for Human Methylation 450K array
#' DT.QC.meta <- load.metharray.QC.meta(array.meta = "controls450")
#' @references Assenov Y. et al., Comprehensive analysis of DNA methylation data
#'             with RnBeads.

load.metharray.QC.meta <- function(array.meta){
  # Get methylation array QC metadata as a data.table
  QC.meta <- RnBeads::rnb.get.annotation(array.meta)
  DT.QC.meta <- data.table::as.data.table(QC.meta)
  DT.QC.meta[, ID := as.character(ID)]
  # Rename target levels to lowercase
  bit::setattr(DT.QC.meta$Target, "levels", vapply(
    X = levels(DT.QC.meta$Target), USE.NAMES = FALSE,
    FUN.VALUE = character(length = 1), FUN = methview.qc::smart.tolower))
  return(DT.QC.meta)
}

#' Merges red and green channels intensities with QC probes metadata.
#'
#' @param RnBSet     An \code{RnBSet} basic object for storing methylation
#'                   array DNA methylation and experimental quality information
#'                   (Bisulfite data not supported).
#'                   \itemize{
#'                    \item{For more information about RnBSet object read
#'                    \link[RnBeads]{RnBSet-class}.}
#'                    \item{To create an RnBSet object run
#'                    \link[RnBeads]{rnb.execute.import}.}
#'                    \item{For additionnal options to import methylation array
#'                    data in the RnBSet see options available in
#'                    \link[RnBeads]{rnb.options}.}
#'                   }
#' @param DT.QC.meta A \code{data.table} with methylation array quality control
#'                   metadata obtained using \link{load.metharray.QC.meta}.
#' @return A \code{data.table} list matching QC metadata with green channel and
#'         red channel intensities.
#' @author Yoann Pageaud.
#' @export merge.QC.intensities.and.meta
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' library(RnBeads)
#' idat.dir <- "~/data/MethylationEPIC/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' rnb.options(identifiers.column = "barcode")
#' #Create the data.table with quality control metadata
#' dt.meta <- load.metharray.QC.meta(array.meta = "controlsEPIC")
#' # Merge red and green channels intensities with QC metadata
#' dt.mrg <- merge.QC.intensities.and.meta(RnBSet = rnb.set, DT.QC.meta = dt.meta)

merge.QC.intensities.and.meta <- function(RnBSet, DT.QC.meta){
  #Get sample IDs
  column.names <- RnBSet@pheno$ID
  #Get QC data
  qc.data <- RnBeads::qc(RnBSet) #Cy3 is Green; Cy5 is Red.
  #Merge Red and Green intensities matrices with QC probes metadata
  QC.data <- lapply(X = names(qc.data), FUN = function(i){
    colnames(qc.data[[i]]) <- column.names
    DT.QC <- data.table::as.data.table(qc.data[[i]], keep.rownames = TRUE)
    DT.QC <- merge(
      x = DT.QC.meta, y = DT.QC, by.x = "ID", by.y = "rn", all = TRUE)
    data.table::setnames(x = DT.QC, old = "ID", new = "QC.probe.IDs")
    DT.QC
  })
  names(QC.data) <- c("Cy3 - Electric Lime Green", "Cy5 - Dark Red")
  #Return QC.data
  return(QC.data)
}


#' Computes ratio of Cy3 & Cy5 intensities with associated color shades.
#'
#' @param DT.probe.ratio A \code{data.table} with 3 columns:
#' \itemize{
#'  \item{column 1 "Samples" must contain all sample IDs.}
#'  \item{column 2 "Cy3 intensity" must contain Cy3 intensity values.}
#'  \item{column 3 "Cy5 intensity" must contain Cy5 intensity values.}
#' }
#' @return A \code{data.table} with 7 columns:
#' \itemize{
#'  \item{column 1 "Samples" contains all sample IDs.}
#'  \item{column 2 "Cy3 intensity" contains Cy3 intensity values.}
#'  \item{column 3 "Cy5 intensity" contains Cy5 intensity values.}
#'  \item{column 4 "Intensity ratio" contains Cy3/Cy5 or Cy5/Cy3 ratio values.}
#'  \item{column 5 "Ratio.type" specify if the ratio is Cy3/Cy5 or Cy5/Cy3.}
#'  \item{column 6 "color.ratio" contains colors associated to ratio values.}
#'  \item{column 7 "round.ratio" contains the rounded ratio values.}
#' }
#' @author Yoann Pageaud.
#' @export
#' @keywords internal

compute.intensity.ratio <- function(DT.probe.ratio){
  #Calculate ratio Cy5/Cy3
  DT.probe.ratio[, c("Intensity ratio", "Ratio.type") := .(
    `Cy5 intensity`/`Cy3 intensity`, "Cy5/Cy3")]
  #Compute ratio color
  DT.probe.ratio[, color.ratio := unlist(lapply(
    X = DT.probe.ratio$`Intensity ratio`, FUN = function(j){
      if(!is.na(j)){
        grDevices::colorRampPalette(
          colors = c("#ff0000", "yellow", "#00ff00"))(round(j))[2]
      } else { "white" }
    }))]
  #If some ratio < 1 try reverse ratio Cy3/Cy5
  if(nrow(DT.probe.ratio[`Intensity ratio` < 1]) > 0){
    DT.probe.ratio[
      `Intensity ratio` < 1, c("Intensity ratio", "Ratio.type") := .(
        `Cy3 intensity`/`Cy5 intensity`, "Cy3/Cy5")]
    DT.probe.ratio[Ratio.type == "Cy3/Cy5", color.ratio := unlist(lapply(
      X = DT.probe.ratio[Ratio.type == "Cy3/Cy5"]$`Intensity ratio`,
      FUN = function(j){
        rev(grDevices::colorRampPalette(
          colors = c("#ff0000", "yellow", "#00ff00"))(round(j)))[2]
      }))]
  }
  #Calculate round ratio
  DT.probe.ratio[, round.ratio := round(`Intensity ratio`)]
  #If some ratio close to 1
  if(nrow(DT.probe.ratio[round.ratio == 1]) > 0){
    DT.probe.ratio[round.ratio == 1, color.ratio := grDevices::colorRampPalette(
      colors = c("#ff0000", "yellow", "#00ff00"))(3)[2]]
  }
  #If some ratio close to 2
  if(nrow(DT.probe.ratio[round.ratio == 2]) > 0){
    DT.probe.ratio[round.ratio == 2 & Ratio.type == "Cy5/Cy3", color.ratio :=
                     grDevices::colorRampPalette(
                       colors = c("#ff0000", "yellow", "#00ff00"))(4)[2]]
    DT.probe.ratio[round.ratio == 2 & Ratio.type == "Cy3/Cy5",
                   color.ratio := rev(grDevices::colorRampPalette(
                     colors = c("#ff0000", "yellow", "#00ff00"))(4))[2]]
  }
  #If some ratio close to 3
  if(nrow(DT.probe.ratio[round.ratio == 3]) > 0){
    DT.probe.ratio[
      round.ratio == 3 & Ratio.type == "Cy5/Cy3", color.ratio :=
        grDevices::colorRampPalette(
          colors = c("#ff0000", "yellow", "#00ff00"))(4)[2]]
    DT.probe.ratio[
      round.ratio == 3 & Ratio.type == "Cy3/Cy5", color.ratio :=
        rev(grDevices::colorRampPalette(
          colors = c("#ff0000", "yellow", "#00ff00"))(4))[2]]
  }
  return(DT.probe.ratio)
}


#' Provides expected intensity for a given methylation array probe ID.
#'
#' @param DT.QC.meta    A \code{data.table} with methylation array quality
#'                      control metadata obtained using
#'                      \link{load.metharray.QC.meta}.
#' @param probe.id      A \code{character} string matching a methylation array
#'                      probe ID listed in 'DT.QC.meta'.
#' @param channel.names A \code{character} vector specifying the names of the 2
#'                      color channels to use for creating the data.table of
#'                      expected intensity (Default:
#'                      channel.names = c("Cy3 - Electric Lime Green",
#'                      "Cy5 - Dark Red")).
#' @return A \code{data.table} with a 'Channel' column and an
#'         'Expected intensity' column providing the probe intensity expected
#'         for green and red channels.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create the data.table with quality control metadata from HM450K
#' dt.meta <- load.metharray.QC.meta(array.meta = "controls450")
#' #Check expected intensity for QC probe "27630314"
#' get.expected.intensity(DT.QC.meta = dt.meta, probe.id = "27630314")

get.expected.intensity <- function(
  DT.QC.meta, probe.id,
  channel.names = c("Cy3 - Electric Lime Green", "Cy5 - Dark Red")){
  #Create DT.expected.intensity
  if(DT.QC.meta[ID == probe.id]$`Evaluate Green` == "+" &
     DT.QC.meta[ID == probe.id]$`Evaluate Red` == "-" &
     DT.QC.meta[ID == probe.id]$`Expected Intensity` == "High"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("High", "Background"))
  } else if(
    DT.QC.meta[ID == probe.id]$`Evaluate Green` == "+" &
    DT.QC.meta[ID == probe.id]$`Evaluate Red` == "-" &
    DT.QC.meta[ID == probe.id]$`Expected Intensity` == "Medium"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("Medium", "Background"))
  } else if(
    DT.QC.meta[ID == probe.id]$`Evaluate Green` == "+" &
    DT.QC.meta[ID == probe.id]$`Evaluate Red` == "-" &
    DT.QC.meta[ID == probe.id]$`Expected Intensity` == "Low"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("Low", "Background"))
  } else if(
    DT.QC.meta[ID == probe.id]$`Evaluate Green` == "-" &
    DT.QC.meta[ID == probe.id]$`Evaluate Red` == "+" &
    DT.QC.meta[ID == probe.id]$`Expected Intensity` == "High"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("Background", "High"))
  } else if(
    DT.QC.meta[ID == probe.id]$`Evaluate Green` == "-" &
    DT.QC.meta[ID == probe.id]$`Evaluate Red` == "+" &
    DT.QC.meta[ID == probe.id]$`Expected Intensity` == "Background"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("Background", "Background"))
  } else if(
    DT.QC.meta[ID == probe.id]$`Evaluate Green` == "+" &
    DT.QC.meta[ID == probe.id]$`Evaluate Red` == "-" &
    DT.QC.meta[ID == probe.id]$`Expected Intensity` == "Background"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("Background", "Background"))
  } else if(
    DT.QC.meta[ID == probe.id]$`Evaluate Green` == "+" &
    DT.QC.meta[ID == probe.id]$`Evaluate Red` == "+" &
    DT.QC.meta[ID == probe.id]$`Expected Intensity` == "Background"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("Background", "Background"))
  } else if(
    DT.QC.meta[ID == probe.id]$`Evaluate Green` == "-" &
    DT.QC.meta[ID == probe.id]$`Evaluate Red` == "-" &
    DT.QC.meta[ID == probe.id]$`Expected Intensity` == "High"){
    DT.expected.intensity <- data.table::data.table(
      "Channel" = channel.names,
      "Expected intensity" = c("Background", "Background"))
  } else { stop("Don't know how to handle this probe!") }
  #TODO: Maybe the definition of "Background" should be thought again ?
  return(DT.expected.intensity)
}
