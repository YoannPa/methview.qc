
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
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiDataEPIC")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiDataEPIC")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' #Check platform used to generate the dataset.
#' get_platform(rnb.set)

get_platform <- function(RnBSet){
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
#' DT.QC.meta <- load_metharray_QC_meta(array.meta = "controls450")
#' @references Assenov Y. et al., Comprehensive analysis of DNA methylation data
#'             with RnBeads.

load_metharray_QC_meta <- function(array.meta){
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
#'                   metadata obtained using \link{load_metharray_QC_meta}.
#' @return A \code{data.table} list matching QC metadata with green channel and
#'         red channel intensities.
#' @author Yoann Pageaud.
#' @export mergeQC_intensities_and_meta
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiDataEPIC")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiDataEPIC")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' RnBeads::rnb.options(identifiers.column = "barcode")
#' #Create the data.table with quality control metadata
#' dt.meta <- load_metharray_QC_meta(array.meta = "controlsEPIC")
#' # Merge red and green channels intensities with QC metadata
#' dt.mrg <- mergeQC_intensities_and_meta(
#'    RnBSet = rnb.set, DT.QC.meta = dt.meta)

mergeQC_intensities_and_meta <- function(RnBSet, DT.QC.meta){
  # Get sample IDs
  if(is.null(rnb.options()$identifiers.column)){
    column.names <- RnBSet@pheno[, 1]
  } else { column.names <- RnBSet@pheno[, rnb.options()$identifiers.column] }
  # Check that all identifiers provided are unique
  if(any(duplicated(column.names))){
    stop(paste(
      "Identifiers column of the RnBSet contains duplicated elements.",
      "Please specify a column with unique identifiers using",
      "RnBeads::rnb.options(identifiers.column = 'barcode')"))
  }
  #Get QC data
  qc.data <- RnBeads::qc(RnBSet) #Cy3 is Green; Cy5 is Red.
  
  #Check if some probes have missing intensities data in all samples and rise
  # warning if any.
  qc.new <- lapply(X = names(qc.data), FUN = function(i){
    missing.probes <- rownames(qc.data[[i]])[
      apply(X = is.na(qc.data[[i]]), MARGIN = 1, FUN = all)]
    if(length(missing.probes) > 0){
      warning(paste(i, "fluorescence intensity missing for probe IDs:",
                    paste0(paste(missing.probes, collapse = ", "), ".")),
              " Probes removed from final data.table.")
      #Remove empty rows from qc.data[[i]]
      qc.data[[i]][!rownames(qc.data[[i]]) %in% missing.probes,]
    } else { qc.data[[i]] }
  })
  names(qc.new) <- names(qc.data)
  qc.data <- qc.new
  #Merge Red and Green intensities matrices with QC probes metadata
  QC.data <- lapply(X = names(qc.data), FUN = function(i){
    colnames(qc.data[[i]]) <- column.names
    DT.QC <- data.table::as.data.table(qc.data[[i]], keep.rownames = TRUE)
    DT.QC <- merge(
      x = DT.QC.meta, y = DT.QC, by.x = "ID", by.y = "rn", all.y = TRUE)
    data.table::setnames(x = DT.QC, old = "ID", new = "QC.probe.IDs")
    DT.QC
  })
  names(QC.data) <- c("Cy3 - Green", "Cy5 - Red")
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
#' @keywords internal

compute_intensity_ratio <- function(DT.probe.ratio){
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
#'                      \link{load_metharray_QC_meta}.
#' @param probe.id      A \code{character} string matching a methylation array
#'                      probe ID listed in 'DT.QC.meta'.
#' @param channel.names A \code{character} vector specifying the names of the 2
#'                      color channels to use for creating the data.table of
#'                      expected intensity (Default:
#'                      channel.names = c("Cy3 - Green", "Cy5 - Red")).
#' @return A \code{data.table} with a 'Channel' column and an
#'         'Expected intensity' column providing the probe intensity expected
#'         for green and red channels.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create the data.table with quality control metadata from HM450K
#' dt.meta <- load_metharray_QC_meta(array.meta = "controls450")
#' #Check expected intensity for QC probe "27630314"
#' get_expected_intensity(DT.QC.meta = dt.meta, probe.id = "27630314")

get_expected_intensity <- function(
  DT.QC.meta, probe.id,
  channel.names = c("Cy3 - Green", "Cy5 - Red")){
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


#' Updates a target intensities metadata.
#' 
#' @param QC.data    A \code{data.table} list matching QC metadata with green
#'                   channel and red channel intensities, obtained with the
#'                   function \link{mergeQC_intensities_and_meta}.
#' @param DT.QC.meta A \code{data.table} with methylation array quality control
#'                   metadata obtained with the function
#'                   \link{load_metharray_QC_meta}.
#' @param target     A \code{character} specifying the name of the target step
#'                   for which the metadata should be updated.
#' @param ncores     An \code{integer} specifying the number of cores or threads
#'                   to be used for parallel processing.      
#' @return A \code{data.table} containing the updated metadata for a given
#'         target step of a given quality control dataset.
#' @author Yoann Pageaud.
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiDataEPIC")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiDataEPIC")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' RnBeads::rnb.options(identifiers.column = "barcode")
#' #Create the data.table with quality control metadata
#' dt.meta <- load_metharray_QC_meta(array.meta = "controlsEPIC")
#' # Merge red and green channels intensities with QC metadata
#' dt.mrg <- mergeQC_intensities_and_meta(
#'    RnBSet = rnb.set, DT.QC.meta = dt.meta)
#' update_target_meta(QC.data = QC.data, target = "Hybridization")
#' @keywords internal

update_target_meta <- function(QC.data, DT.QC.meta, target, ncores = 1){
  #Melt Green & Red QC data
  DT.target.Cy3 <- data.table::melt.data.table(
    data = QC.data$`Cy3 - Green`[Target == target],
    measure.vars = colnames(QC.data$`Cy3 - Green`)[-c(1:10)],
    variable.name = "Samples", value.name = "Cy3 intensity")
  DT.target.Cy5 <- data.table::melt.data.table(
    data = QC.data$`Cy5 - Red`[Target == target],
    measure.vars = colnames(QC.data$`Cy5 - Red`)[-c(1:10)],
    variable.name = "Samples", value.name = "Cy5 intensity")
  #Rbind data.tables
  ls.dt.target <- list(DT.target.Cy3, DT.target.Cy5)
  names(ls.dt.target) <- names(QC.data)
  DT.target <- data.table::rbindlist(
    l = ls.dt.target, idcol = "Cyanine", use.names = FALSE)
  #Check expected intensities for each probes
  ls.exp.intens <- parallel::mclapply(
    X = unique(DT.target$QC.probe.IDs), mc.cores = ncores, FUN = function(i){
      methview.qc::get_expected_intensity(
        DT.QC.meta = DT.QC.meta, probe.id = i, channel.names = names(QC.data))
    })
  names(ls.exp.intens) <- unique(DT.target$QC.probe.IDs)
  DT.exp.intens <- data.table::rbindlist(l = ls.exp.intens, idcol = "Probe.ID")
  #Modify DT.target with expected intensities
  invisible(lapply(X = seq(nrow(DT.exp.intens)), FUN = function(i){
    DT.target[QC.probe.IDs == DT.exp.intens[i,]$Probe.ID &
                Cyanine == DT.exp.intens[i,]$Channel,
              `Expected Intensity` := DT.exp.intens[i,]$`Expected intensity`]
  }))
  #Change order of levels in expected intensity
  DT.target[, `Expected Intensity` := factor(
    `Expected Intensity`, levels = levels(`Expected Intensity`)[c(2, 4, 3, 1)])]
  
  return(DT.target)
}

#' Computes a PCA from an RnBSet on quality control probes intensities
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
#' @return A \code{list} containing a prcomp object with all results from the
#'         PCA, and a data.table with all the RnBSet data.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create an RnBSet for MethylationEPIC data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiDataEPIC")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiDataEPIC")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' # Compute PCA on quality control probes intensities from the RnBSet
#' comp_RnBqc2PCA(RnBSet = rnb.set)     

comp_RnBqc2PCA <- function(RnBSet){
  if(methview.qc::get_platform(RnBSet = RnBSet) == "MethylationEPIC"){
    DT.QC.meta <- methview.qc::load_metharray_QC_meta(
      array.meta = "controlsEPIC")
  } else if(methview.qc::get_platform(RnBSet = RnBSet) == "HM450K"){
    DT.QC.meta <- methview.qc::load_metharray_QC_meta(
      array.meta = "controls450")
  }
  #Merge Red and Green intensities matrices with QC probes metadata
  QC.data <- methview.qc::mergeQC_intensities_and_meta(
    RnBSet = RnBSet, DT.QC.meta = DT.QC.meta)
  QC.data <- data.table::rbindlist(
    l = QC.data, use.names = TRUE, idcol = "Channel")
  QC.data[, Target := as.factor(Target)]
  melt.QC.dt <- data.table::melt(
    QC.data, id.vars = colnames(QC.data)[1:11], variable.name = "Samples")
  t.QC.dt <- data.table::dcast(
    melt.QC.dt, formula = Samples ~ Channel + Description)
  RnBSet@pheno[, 1] <- as.factor(RnBSet@pheno[, 1])
  t.QC.dt <- data.table::merge.data.table(
    x = RnBSet@pheno, y = t.QC.dt,by.x = colnames(RnBSet@pheno)[1],
    by.y = "Samples", all.y = TRUE)
  pca_t.res <- prcomp(t.QC.dt[, -c(1:ncol(RnBSet@pheno)), ], scale. = FALSE)
  #Return the results of the PCA
  ls_res <- list("prcomp" = pca_t.res, "data" = t.QC.dt)
  return(ls_res)
}

#' Computes a deviation score between samples fluorescence and an internal
#' HM450K reference.
#'
#' @param RnBSet  An \code{RnBSet} basic object for storing HM450K DNA
#'                methylation and experimental quality information (Bisulfite
#'                data are not supported).
#'                \itemize{
#'                 \item{For more information about RnBSet object read
#'                 \link[RnBeads]{RnBSet-class}.}
#'                 \item{To create an RnBSet object run
#'                 \link[RnBeads]{rnb.execute.import}.}
#'                 \item{For additionnal options to import methylation array
#'                 data in the RnBSet see options available in
#'                 \link[RnBeads]{rnb.options}.}
#'                }
#' @param samples A \code{character} vector specifying the samples to include
#'                for the deviation score calculation. You can catch the sample
#'                IDs you wish to evaluate running \code{RnBSet@pheno[,1]}
#' @param target  A \code{character} specifying the name of the target step to
#'                evaluate fluorescence from.
#' @param ncores  An \code{integer} specifying the number of cores or threads to
#'                be used for parallel processing.
#' @return A \code{data.table} with probes data and samples computed deviation
#'         scores.
#' @details The deviation score is calculated as following:
#'          \enumerate{
#'           \item{the square root of samples intensities fluorescence is
#'           calculated for both red & green channels.}
#'           \item{the square root of the reference intensities fluorescence is
#'           calculated the same way.}
#'           \item{Then the ratio of the samples values out of the reference
#'           values is calculated.}
#'           \item{The previous result is multiplied by 100 to get a percentage.
#'           }
#'           \item{Then 100 is substracted to the result in order to obtain the
#'           delta between this percentage and the reference (which counts as
#'           100\% signal intensity). The resulting values are deviation scores
#'           of fluorescence computed for each probes and each samples. It can
#'           be extracted from the resulting data.table in the column
#'           'percent.diff.sqrt'.}
#'          }
#'          Formula: \eqn{(sqrt(sample\_signal)/sqrt(ref\_signal)) * 100 - 100}
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create an RnBSet for Human Methylation 450K data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiData")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiData")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' # Compute deviation score on selected samples for hybridization QC probes
#' devscore.fluo(
#'     RnBSet = rnb.set, samples = c("GroupA_3","GroupA_2","GroupB_2"),
#'     target = "Hybridization")

devscore.fluo <- function(RnBSet, samples, target, ncores = 1){
  platform <- methview.qc::get_platform(RnBSet = RnBSet)
  #Check it is HM450K
  if(platform == "HM450K"){
    #Load quality control metadata for Human Methylation 450K array
    DT.QC.meta <- methview.qc::load_metharray_QC_meta(array.meta = "controls450")
  } else if(platform == "MethylationEPIC"){
    #Load quality control metadata for Methylation EPIC array
    DT.QC.meta <- methview.qc::load_metharray_QC_meta(array.meta = "controlsEPIC")
  } else { stop("devscore.fluo() only supports HM450K and EPIC data.") }
  QC.data <- methview.qc::mergeQC_intensities_and_meta(
    RnBSet = RnBSet, DT.QC.meta = DT.QC.meta)
  #Keep only samples requested
  QC.data <- lapply(X = QC.data, FUN = function(i){
    keep <- c(colnames(i)[1:10], samples)
    i[, ..keep]
  })
  #Update target metadata
  DT.target <- methview.qc:::update_target_meta(
    QC.data = QC.data, DT.QC.meta = DT.QC.meta, target = target,
    ncores = ncores)
  #Create unique combination cyanine & ID
  DT.target[, Cyanine.probe.ID := paste(Cyanine, QC.probe.IDs, sep = ".")]
  if(platform == "HM450K"){
    hm450ref <- methview.qc:::HM450K_sysdata[Target == target]
    hm450ref[, Cyanine.probe.ID := paste(Channel, QC.probe.IDs, sep = ".")]
    #Map HM450K reference fluorescence to DT.target
    DT.target <- data.table::merge.data.table(
      x = DT.target,
      y = hm450ref[, c("Cyanine.probe.ID", "PCAWG.avg.intensity"), ],
      by = "Cyanine.probe.ID", all.x = TRUE)
    #Compute percentage difference of the square root of fluorescence intensities
    DT.target[, percent.diff.sqrt := (
      (sqrt(`Cy3 intensity`)/sqrt(PCAWG.avg.intensity))*100) - 100]
  } else if(platform == "MethylationEPIC"){
    epicref <- methview.qc:::EPIC_sysdata[Target == target]
    epicref[, Cyanine.probe.ID := paste(Channel, QC.probe.IDs, sep = ".")]
    #Map HM450K reference fluorescence to DT.target
    DT.target <- data.table::merge.data.table(
      x = DT.target,
      y = epicref[, c("Cyanine.probe.ID", "GSE197678.avg.intensity"), ],
      by = "Cyanine.probe.ID", all.x = TRUE)
    #Compute percentage difference of the square root of fluorescence intensities
    DT.target[, percent.diff.sqrt := (
      (sqrt(`Cy3 intensity`)/sqrt(GSE197678.avg.intensity))*100) - 100]
  }
  #Make percentage difference table
  DT.target <- DT.target[, c(2:13, 16), ]
  return(DT.target)
}

#' Retrieves runinfo data from a methylation array sample's IDAT files.
#' 
#' @param sentrix_barcode A \code{character} string to specify the sample from
#'                        which you wish to retrieve the runinfo (from both Red
#'                        & Green IDAT files).
#' @param IDATs_dir       A \code{character} specifying the directory where IDAT
#'                        files are stored that will be search recursively for
#'                        the files matching the sentrix barcode. 
#' @param data_format     A \code{character} string to specify how much
#'                        information you want from runinfo:
#'                        \itemize{
#'                         \item{'full': to get all details available from
#'                         the sample's IDAT files runinfo.}
#'                         \item{'short': to get a one line summary of runinfo
#'                         data from the sample's IDAT files (Scan year,
#'                         Decoding tool version & Scanning tool version).}
#'                         \item{'both': the default value, to get both the full
#'                         details and the one line summary of runinfo data from
#'                         the sample's IDAT files in 2 separate data.tables.}
#'                        }
#' @return A \code{data.table} or a \code{list} of data.tables with the sample's
#'         runinfo data.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Set IDAT directory
#' idat.dir <- system.file("extdata", "idat", package = "IlluminaDataTestFiles")
#' # Get runinfo data from a IDAT file name from within the directory
#' get_IDATs_runinfo(sentrix_barcode = "4019585376_B_Red", IDATs_dir = idat.dir)
#' @references ML Smith, KA Baggerly, H Bengtsson, ME Ritchie & KD Hansen.
#'             illuminaio: An open source IDAT parsing tool for Illumina
#'             microarrays, F1000Research, (2013) 2:264.

get_IDATs_runinfo <- function(sentrix_barcode, IDATs_dir, data_format = "both"){
  files <- list.files(IDATs_dir, recursive = TRUE)
  IDATs <- files[files %like% sentrix_barcode]
  if(length(IDATs) != 2){
    stop("There should be 2 IDAT files per sample ID.")
  } else {
    ls_content <- lapply(
      X = file.path(IDATs_dir, IDATs), FUN = illuminaio::readIDAT)
    if(all.equal(
      target = ls_content[[1]]$RunInfo,
      current = ls_content[[2]]$RunInfo)){
      dt_scan_info <- as.data.table(ls_content[[1]]$RunInfo)
      if(nrow(dt_scan_info) != 0){
        dt_scan_info[, RunTime := as.POSIXct(
          RunTime, format = "%m/%d/%Y %H:%M:%S")]
        dt_scan_info[, Tool := paste(BlockCode, CodeVersion)]
        if(data_format == 'short' | data_format == 'both'){
          scan_year <- format(x = dt_scan_info$RunTime, "%Y")
          if(length(unique(scan_year)) != 1){
            scan_year <- unique(format(
              x = dt_scan_info[BlockType == "Scan"]$RunTime,
              "%Y"))
            if(length(scan_year) != 1){
              stop("All scans did not happened the same year.")
            }
          } else { scan_year <- unique(scan_year) }
          
          decoding <- dt_scan_info[BlockType == "Decoding"]$Tool
          if(length(unique(decoding)) != 1){
            stop("Decoding has been done with different tools or versions.")
          } else { decoding <- unique(decoding) }
          
          scan_tool <- dt_scan_info[BlockType == "Scan"]$Tool
          if(length(unique(scan_tool)) != 1){
            stop("Scan has been done with different tools or versions.")
          } else { scan_tool <- unique(scan_tool) }
          dt_scan_simple <- data.table(
            "Scan_year" = scan_year, "Decoding" = decoding,
            "Scan" = scan_tool)
        }
        dt_scan_info <- dt_scan_info[, -c("Tool"), ]
        
      } else {
        warning(
          "The sample's IDAT files do not contain any runinfo data.")
        class(dt_scan_info$RunTime) <- class(as.POSIXct(NA))
        dt_scan_info <- rbind(
          dt_scan_info, data.table(as.POSIXct(NA), NA, NA, NA, NA),
          use.names = FALSE)
        if(data_format == 'short' | data_format == 'both'){
          dt_scan_simple <- data.table(
            "Scan_year" = NA, "Decoding" = NA, "Scan" = NA)
        }
      }
    } else { stop("Green and Red IDATs runinfos do not match.") }
  }
  if(data_format == 'short'){
    return(dt_scan_simple)
  } else if(data_format == 'full'){
    return(dt_scan_info)
  } else if(data_format == 'both'){
    ls_both <- list(
      "short_runinfo" = dt_scan_simple, "full_runinfo" = dt_scan_info)
    return(ls_both)
  }
}

#' Prepares annotations to be tested for associations.
#'
#' @param RnBSet  An \code{RnBSet} basic object for storing methylation array
#'                DNA methylation and experimental quality information
#'                (Bisulfite data not supported).
#'                \itemize{
#'                 \item{For more information about RnBSet object read
#'                 \link[RnBeads]{RnBSet-class}.}
#'                 \item{To create an RnBSet object run
#'                 \link[RnBeads]{rnb.execute.import}.}
#'                 \item{For additionnal options to import methylation array
#'                 data in the RnBSet see options available in
#'                 \link[RnBeads]{rnb.options}.}
#'                }
#' @param verbose A \code{logical} to display information about the step-by-step
#'                processing of the data if TRUE (Default: verbose = FALSE).
#' @return A \code{list} containing updated annotations, the number of
#'         annotations available, and the RnBSet annotation table it contains.
#' @author Yoann Pageaud.
#' @examples
#' # Create an RnBSet for Human Methylation 450K data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiData")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiData")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' # Prepare annotations before association tests
#' methview.qc:::prepare_annot_asso(RnBSet = rnb.set)
#' @keywords internal

prepare_annot_asso <- function(RnBSet, verbose = FALSE){
  annot.table <- pheno(RnBSet)
  annots <- lapply(X = annot.table, FUN = function(x) {
    # Exclude columns containing the same value for all rows
    # (except when missing)
    if(length(unique(na.omit(x))) < 2){ return(NULL) }
    # Convert as factor character variables
    if(is.character(x)){ x <- as.factor(x) }
    # Convert as factor TRUE/FALSE variables
    if(is.logical(x)){ x <- as.factor(x) }
    if(is.factor(x)){ if(anyDuplicated(na.omit(x)) == 0) { return(NULL) } }
    x
  })
  annots <- annots[!vapply(
    X = annots, FUN = is.null, FUN.VALUE = logical(length = 1L))]
  n.annot <- length(annots)
  if(n.annot == 0){
    stop("No suitable annotation found for association tests.")}
  if(verbose){
    cat(c("Testing the following annotations for associations:\n",
          paste(names(annots), collapse = ", "), ".\n"))
  }
  res <- list(
    "annotations" = annots, "n.annot" = n.annot,
    "annot.table" = annot.table)
  return(res)
}

#' Tests association of an annotation with another one or with a PC.
#'
#' @param x           A \code{vector} of numerics or factors containing data
#'                    from one annotation.
#' @param y           A \code{vector} of numerics or factors containing data
#'                    from a PCA principal component or another annotation.
#' @param perm.matrix A \code{matrix} containing permutations of the order of
#'                    rows, as integers, based on annotations length. The
#'                    permutation matrix can optionally be used to test the
#'                    significance of a correlation value between 2 annotations
#'                    or an annotation and a PCA principal component.
#' @return A \code{data.table} containing all results from an association test
#'         between x and y.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Most basic statistic test between 2 vectors of integers
#' test.annots(x = 1:15, y = 24:10)

test.annots <- function (x, y, perm.matrix = NULL){
  # Set nominal parameter
  nominal <- TRUE
  #If annotation x or y is a date then convert it into integers
  if(class(x) == "Date"){ x <- as.integer(x) }
  if(class(y) == "Date"){ y <- as.integer(y) }
  inds <- which(!(is.na(x) | is.na(y)))
  if (length(inds) < 2) {
    nominal <- FALSE
    warning("Not enough common values.")
  }
  x <- x[inds]
  if (is.factor(x)) {
    x <- as.factor(as.character(x))
    if (nlevels(x) < 2) {
      nominal <- FALSE
      warning("Not enough categories.")
    }
  }
  y <- y[inds]
  if (is.factor(y)) {
    y <- as.factor(as.character(y))
    if (nlevels(y) < 2) {
      nominal <- FALSE
      warning("Not enough categories.")
    }
  }
  if(nominal){ # If a stat. test can be run
    # Retrieve test p-value
    get.p <- function(expr) {
      tryCatch(suppressWarnings(expr$p.value), error = function(er){
        as.double(NA)
      })
    }
    if(is.factor(x)){
      if(is.factor(y)) {
        # If both vectors are factors do Non-parametric Fisher test with
        # 50 000 replicates for Monte Carlo test.
        simulate <- (nlevels(x) > 2 || nlevels(y) > 2)
        dt.res <- data.table::data.table(
          test = "Fisher", correlation = NA, pvalue = get.p(
            fisher.test(
              x, y, conf.int = FALSE, simulate.p.value = simulate,
              B = 50000)))
      } else if (nlevels(x) == 2) {
        # If Y is a numeric vector and X has only 2 possible values do
        # Non-parametric Wilcoxon–Mann–Whitney test
        # aka Wilcoxon rank-sum test
        values <- tapply(y, x, identity)
        dt.res <- data.table::data.table(
          test = "Wilcoxon-Mann-Whitney", correlation = NA,
          pvalue = get.p(wilcox.test(
            values[[1]], values[[2]], alternative = "two.sided"))
        )
      } else {
        # If X has multiple levels & Y is a numeric vector do
        # Non-parametric Kruskal-Wallis test
        dt.res <- data.table::data.table(
          test = "Kruskal-Wallis", correlation = NA, pvalue = get.p(
            kruskal.test(y, x)))
      }
    } else if(is.factor(y)) {
      if (nlevels(y) == 2) {
        # If X is a numeric vector and Y has only 2 possible values do
        # Non-parametric Wilcoxon–Mann–Whitney test
        # aka Wilcoxon rank-sum test
        values <- tapply(x, y, identity)
        dt.res <- data.table::data.table(
          test = "Wilcoxon-Mann-Whitney", correlation = NA,
          pvalue = get.p(wilcox.test(
            values[[1]], values[[2]], alternative = "two.sided")))
      } else {
        # If X is a numeric vector and Y has multiple levels do
        # Non-parametric Kruskal-Wallis test
        dt.res <- data.table::data.table(
          test = "Kruskal-Wallis", correlation = NA, pvalue = get.p(
            kruskal.test(x, y)))
      }
    } else {
      # If both X & Y are numeric vectors do a correlation test
      N <- length(inds)
      if(is.null(perm.matrix)){
        cor_res <- cor(x, y) # Calculate Pearson correlation
        # Get T statistic
        t_stat <- (cor_res*sqrt(N-2))/(sqrt(1-cor_res^2))
        dt.res <- data.table::data.table(
          test = "Pearson Corr.", correlation = cor_res,
          pvalue = 2*pt(-abs(t_stat), N-2))
      } else {
        values <- apply(X = perm.matrix, MARGIN = 2, FUN = function(i){
          cor(x[i[i <= N]], y)
        })
        abs_val <- abs(values)
        # Correlation p-value here is the proportion of times where the
        # 1st correlation value is inferior or equal to calculated
        # correlations in all permutations. 10 000 permutations to make
        # sure that 1 result out of 10 000 will be significant enough if
        # correlated.
        dt.res <- data.table::data.table(
          test = "Pearson Corr.", correlation = values[1],
          pvalue = mean(abs_val[1] <= abs_val))
      }
    }
  } else { # Empty data.table result when no test possible
    dt.res <- data.table::data.table(
      test = NA, correlation = NA, pvalue = NA)
  }
  return(dt.res)
}

#' Tests associations between annotations from an RnBSet and PCs from a prcomp
#' object.
#'
#' @param RnBSet     A \code{RnBSet} basic object for storing methylation array
#'                   data and experimental quality information (Bisulfite data
#'                   not supported).
#'                   \itemize{
#'                    \item{For more information about RnBSet object read
#'                    \link[RnBeads]{RnBSet-class}.}
#'                    \item{To create an RnBSet object run
#'                    \link[RnBeads]{rnb.execute.import}.}
#'                    \item{For additionnal options to import methylation array
#'                    data in the RnBSet see options available in
#'                    \link[RnBeads]{rnb.options}.}
#'                   }
#' @param prcomp.res A PCA result of classes \code{prcomp} or
#'                   \code{irlba_prcomp} resulting from stats::prcomp() or
#'                   irlba::prcomp_irlba().
#' @param perm.count An \code{integer} specifying the number of permutations to
#'                   realize on a vector, for the permutations matrix
#'                   initialization, to be used for calculating the significance
#'                   of a correlation test (Default: perm.count = 10000).
#' @param max.PCs    An \code{integer} specifying the maximum number of
#'                   principal components to consider for association tests with
#'                   annotations (Default: max.PCs = 8).
#' @param verbose    A \code{logical} to display information about the
#'                   step-by-step processing of the data if TRUE
#'                   (Default: verbose = FALSE).
#' @return A \code{data.table} containing all association test results.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create an RnBSet for Human Methylation 450K data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiData")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiData")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' # Compute PCA on QC probes from the RnBSet     
#' pca_res <- comp_RnBqc2PCA(RnBSet = rnb.set)
#' # Compute association test between annotations and PCs from QC probes
#' res <- rnb_test_asso_annot_PC(RnBSet = rnb.set, prcomp.res = pca_res$prcomp)

rnb_test_asso_annot_PC <- function(
  RnBSet, prcomp.res, perm.count = 10000, max.PCs = 8, verbose = FALSE){
  # Prepare the RnBset annotations
  prep_res <- prepare_annot_asso(RnBSet = RnBSet, verbose = verbose)
  annots <- prep_res$annotations
  n.annot <- prep_res$n.annot
  annot.table <- prep_res$annot.table
  
  if(perm.count != 0 && (
    (!is.null(prcomp.res)) || sum(!vapply(
      X = annots, FUN = is.factor,
      FUN.VALUE = logical(length = 1L))) >= 2)) {
    # Create the random permutation matrix
    perm.matrix <- mapply(
      FUN = sample, rep(nrow(annot.table), times = perm.count))
    perm.matrix[, 1] <- 1:nrow(perm.matrix)
  } else {
    warning("Cannot initialize the permutations matrix.")
    perm.matrix <- NULL
  }
  
  pc.association.count <- max.PCs
  # Get PCs % variance explained
  PCA_metrics <- BiocompR::prepare_pca_data(
    prcomp.res = prcomp.res, dt.annot = annot.table, PCs = 1:max.PCs,
    scale = 1)
  dpoints <- prcomp.res$x
  if (!is.null(dpoints)) {
    if (ncol(dpoints) > pc.association.count) {
      # Reduce principal components coordinates to the maximum number of
      # dimensions wanted
      dpoints <- dpoints[, 1:pc.association.count]
    }
    # Test all annotations against all PCs
    ls_allres <- lapply(X = seq(n.annot), FUN = function(i){
      if(verbose){ cat("Testing association of", names(annots)[i], "&\n") }
      ls_tres <- lapply(X = seq(ncol(dpoints)), FUN = function(j){
        if(verbose){ cat("\t", colnames(dpoints)[j], "\n") }
        t.result <- test.annots(
          x = annots[[i]], y = dpoints[, j],
          perm.matrix = perm.matrix)
        t.result[, c("annotation", "PC", "var.explained") := .(
          names(annots)[i], colnames(dpoints)[j],
          (PCA_metrics$var.explained*100)[j])]
        t.result
      })
      data.table::rbindlist(l = ls_tres)
    })
    # Rbind all results
    dt_allres <- data.table::rbindlist(l = ls_allres)
    dt_allres[, log_trans_pval := -log10(pvalue)]
    rm(ls_allres)
    # Convert PC as factor to keep the right order
    dt_allres[, PC := as.factor(x = PC)]
    dt_allres[, PC := factor(
      x = PC, levels = paste0("PC", seq(pc.association.count)))]
  }
  rm(dpoints)
  # Return association test results
  return(dt_allres)
}

#' Tests associations between annotations from an RnBSet and QC probes
#' intensities.
#'
#' @param RnBSet       A \code{RnBSet} basic object for storing methylation
#'                     array data and experimental quality information
#'                     (Bisulfite data not supported).
#'                     \itemize{
#'                      \item{For more information about RnBSet object read
#'                      \link[RnBeads]{RnBSet-class}.}
#'                      \item{To create an RnBSet object run
#'                      \link[RnBeads]{rnb.execute.import}.}
#'                      \item{For additionnal options to import methylation
#'                      array data in the RnBSet see options available in
#'                      \link[RnBeads]{rnb.options}.}
#'                     }
#' @param perm.count   An \code{integer} specifying the number of permutations
#'                     to realize on a vector, for the permutations matrix
#'                     initialization, to be used for calculating the
#'                     significance of a correlation test
#'                     (Default: perm.count = 10000).
#' @param max.QCprobes An \code{integer} specifying how many QC probes should be
#'                     kept as the top QC probes associated with RnBSet
#'                     annotations (Default: max.QCprobes = 50).
#' @param verbose      A \code{logical} to display information about the
#'                     step-by-step processing of the data if TRUE
#'                     (Default: verbose = FALSE).
#' @param ncores       An \code{integer} to specify the number of cores/threads
#'                     to be used to parallel-compute association tests between
#'                     annotations and QC probes intensities.
#' @return A \code{data.table} containing all association test results.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create an RnBSet for Human Methylation 450K data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiData")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiData")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' # Compute association tests between annotations and QC probes intensities
#' # (if ncores = 1 execution may take several minutes...)
#' res <- rnb_test_asso_annot_QC(RnBSet = rnb.set)

rnb_test_asso_annot_QC <- function(
  RnBSet, perm.count = 10000, max.QCprobes = 50, verbose = FALSE, ncores = 1){
  # Prepare the RnBset annotations
  prep_res <- methview.qc:::prepare_annot_asso(
    RnBSet = RnBSet, verbose = verbose)
  annots <- prep_res$annotations
  n.annot <- prep_res$n.annot
  annot.table <- prep_res$annot.table
  #Create the data.table with quality control metadata
  if(methview.qc::get_platform(RnBSet = RnBSet) == "HM450K"){
    dt.meta <- load_metharray_QC_meta(array.meta = "controls450")
  } else if(methview.qc::get_platform(RnBSet = RnBSet) == "MethylationEPIC"){
    dt.meta <- load_metharray_QC_meta(array.meta = "controlsEPIC")
  } else{ stop("Methylation array platform unknown.") }
  # Merge red and green channels intensities with QC metadata
  QC.data <- mergeQC_intensities_and_meta(
    RnBSet = RnBSet, DT.QC.meta = dt.meta)
  DTQC <- rbindlist(l = QC.data, idcol = "Channel")
  DTQC[Channel == "Cy3 - Green", Channel := "Green"]
  DTQC[Channel == "Cy5 - Red", Channel := "Red"]
  DTQC[, Probe_name := paste(Description, QC.probe.IDs, Channel, sep = "_")]
  DTQC <- data.table::as.data.table(
    x = t(as.matrix(DTQC[, -c(1:11), ], rownames = "Probe_name")),
    keep.rownames = "ID")
  if(perm.count != 0 && sum(!vapply(
    X = annots, FUN = is.factor, FUN.VALUE = logical(length = 1L))) >= 2) {
    # Create the random permutation matrix
    perm.matrix <- mapply(
      FUN = sample, rep(nrow(annot.table), times = perm.count))
    perm.matrix[, 1] <- 1:nrow(perm.matrix)
  } else {
    warning("Cannot initialize the permutations matrix.")
    perm.matrix <- NULL
  }
  if(all.equal(target = annot.table[[1]], current = DTQC$ID)){
    # Test all annotations against QC probes intensities
    ls_allres <- lapply(X = seq(n.annot), FUN = function(i){
      if(verbose){ cat(
        "Testing association of", names(annots)[i],
        "against all QC probes intensities...") }
      ls_tres <- mclapply(
        X = seq(2, ncol(DTQC)), mc.cores = ncores, FUN = function(j){
          t.result <- test.annots(
            x = annots[[i]], y = DTQC[[j]],
            perm.matrix = perm.matrix)
          t.result[, c("annotation", "QC_probe") := .(
            names(annots)[i], colnames(DTQC)[j])]
          t.result
        })
      if(verbose){ cat("Done.\n") }
      data.table::rbindlist(l = ls_tres)
    })
    # Rbind all results
    dt_allres <- data.table::rbindlist(l = ls_allres)
    dt_allres[, log_trans_pval := -log10(pvalue)]
    rm(ls_allres)
  } else { stop("ID columns in tables don't have the same order.") }
  # Get top significant QC probes
  dt_allres[, min_pval := min(pvalue, na.rm = TRUE), by = QC_probe]
  top_probes <- utils::head(
    x = unique(dt_allres, by = "QC_probe")[order(min_pval)],
    n = max.QCprobes)$QC_probe
  dt_allres <- dt_allres[QC_probe %in% top_probes]
  # Convert PC as factor to keep the right order
  dt_allres[, QC_probe := as.factor(x = QC_probe)]
  dt_allres[, QC_probe := factor(x = QC_probe, levels = top_probes)]
  return(dt_allres)
}

#' Tests associations between all annotations in a RnBSet.
#'
#' @param RnBSet       A \code{RnBSet} basic object for storing methylation
#'                     array data and experimental quality information
#'                     (Bisulfite data not supported).
#'                     \itemize{
#'                      \item{For more information about RnBSet object read
#'                      \link[RnBeads]{RnBSet-class}.}
#'                      \item{To create an RnBSet object run
#'                      \link[RnBeads]{rnb.execute.import}.}
#'                      \item{For additionnal options to import methylation
#'                      array data in the RnBSet see options available in
#'                      \link[RnBeads]{rnb.options}.}
#'                     }
#' @param perm.count   An \code{integer} specifying the number of permutations
#'                     to realize on a vector, for the permutations matrix
#'                     initialization, to be used for calculating the
#'                     significance of a correlation test
#'                     (Default: perm.count = 10000).
#' @param verbose      A \code{logical} to display information about the
#'                     step-by-step processing of the data if TRUE
#'                     (Default: verbose = FALSE).
#' @return A \code{data.table} containing all association test results.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Create an RnBSet for Human Methylation 450K data
#' require(Biobase)
#' idat.dir <- system.file("extdata", package = "minfiData")
#' sample.annotation <- system.file(
#'     "extdata", "SampleSheet.csv", package = "minfiData")
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- RnBeads::rnb.execute.import(
#'     data.source = data.source, data.type = "idat.dir")
#' # Compute association tests between annotations and QC probes intensities
#' res <- rnb_test_asso_all_annot(RnBSet = rnb.set)

rnb_test_asso_all_annot <- function(
  RnBSet, perm.count = 10000, verbose = FALSE){
  
  # Prepare the RnBset annotations
  prep_res <- prepare_annot_asso(RnBSet = RnBSet, verbose = verbose)
  annots <- prep_res$annotations
  n.annot <- prep_res$n.annot
  annot.table <- prep_res$annot.table
  
  if(perm.count != 0 && sum(!vapply(
    X = annots, FUN = is.factor, FUN.VALUE = logical(length = 1L))) >= 2) {
    # Create the random permutation matrix
    perm.matrix <- mapply(
      FUN = sample, rep(nrow(annot.table), times = perm.count))
    perm.matrix[, 1] <- 1:nrow(perm.matrix)
  } else {
    warning("Cannot initialize the permutations matrix.")
    perm.matrix <- NULL
  }
  
  if (n.annot > 1) {
    # Create matrix of tests combinations
    test_matrix <- utils::combn(x = names(annots), m = 2)
    # Test association between all annotations available
    ls_annotres <- apply(X = test_matrix, MARGIN = 2, FUN = function(i){
      if(verbose){ cat("Testing association of", i[1], "&", i[2], "\n") }
      t.result <- test.annots(
        x = annots[[i[1]]], y = annots[[i[2]]],
        perm.matrix = perm.matrix)
      t.result[, c("annotation1", "annotation2") := .(i[1], i[2])]
    })
    dt_annotres <- data.table::rbindlist(l = ls_annotres)
    #Duplicate results for the full table
    dt_annotres_bis <- data.table::copy(dt_annotres)
    data.table::setnames(
      x = dt_annotres_bis, old = "annotation1", new = "annotation2_new")
    data.table::setnames(
      x = dt_annotres_bis, old = "annotation2", new = "annotation1")
    data.table::setnames(
      x = dt_annotres_bis, old = "annotation2_new", new = "annotation2")
    # Rbind all results
    dt_annotres <- rbind(dt_annotres, dt_annotres_bis, use.names = TRUE)
    rm(dt_annotres_bis, perm.matrix)
    dt_annotres[, log_trans_pval := -log10(pvalue)]
  }
  return(dt_annotres)
}
