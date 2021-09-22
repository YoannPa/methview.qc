
#' Draws a heatmap from methylation array genotyping probes.
#'
#' @param RnBSet             A \code{RnBSet} basic object for storing
#'                           methylation array data and experimental quality
#'                           information (Bisulfite data not supported).
#'                           \itemize{
#'                            \item{For more information about RnBSet object
#'                            read \link[RnBeads]{RnBSet-class}.}
#'                            \item{To create an RnBSet object run
#'                            \link[RnBeads]{rnb.execute.import}.}
#'                            \item{For additionnal options to import
#'                            methylation array data in the RnBSet see options
#'                            available in \link[RnBeads]{rnb.options}.}
#'                           }
#' @param dist.method        A \code{character} vector to specify the distance
#'                           methods to be used on the matrix rows and columns.
#'                           \itemize{
#'                            \item{If the vector is of length 1: the given
#'                                  method will be applied on rows and columns
#'                                  of the matrix.}
#'                            \item{If the vector is of length 2: the first
#'                                  method will be applied on matrix rows, and
#'                                  the second method will be applied on matrix
#'                                  columns.}
#'                           }
#'                           (Default: dist.method = 'manhattan';
#'                           Supported: dist.method = c('euclidean', 'maximum',
#'                           'manhattan', 'canberra', 'binary', 'minkowski',
#'                           'none')).
#' @param annot.grps         A \code{list} of vectors of groups to which
#'                           variables belongs for the annotation sidebars.
#'                           Vectors' lengths have to match the number of
#'                           variables.
#' @param annot.pal          A \code{vector} or a list of vectors containing
#'                           colors as characters for the annotation sidebars.
#'                           The length of vectors has to match the number of
#'                           levels of vectors listed in 'annot.grps'.
#'                           \itemize{
#'                            \item{If annot.pal is a list: its length must
#'                                  match the length of the list provided to
#'                                  'annot.grps'.}
#'                            \item{If annot.pal is a vector: make sure that the
#'                                  levels content of annotations listed in
#'                                  'annot.grps' is the same, and that no
#'                                  annotation contains less or more levels than
#'                                  another one in 'annot.grps'.}
#'                           }
#' @param heatmap.pal        A \code{character} vector to be used as a color
#'                           palette for the heatmap.
#' @param x.lab              A \code{character} to specify X-axis title
#'                           (Default: x.lab = 'Samples').
#' @param plot.title         A \code{character} to be used as title for the
#'                           plot (Default: plot.title = NULL ; ).
#' @param htmp.text.x        An \code{element_text} object to setup X axis text 
#'                           of the heatmap (Default: htmp.text.x = 
#'                           element_text(size = 10, angle = -45, hjust = 0,
#'                           vjust = 0.5, color = "black")).
#' @param htmp.text.y.right  An \code{element_text} object to setup right Y axis
#'                           text of the heatmap (Default: htmp.text.y.right = 
#'                           element_text(size = 7, color = "black")).
#' @param htmp.title.y.right An \code{element_text} object to setup right Y axis
#'                           title of the heatmap (Default: htmp.title.y.right = 
#'                           element_text(size = 11)).
#' @param anno.text.y.right  An \code{element_text} object to setup right Y axis
#'                           text of the top annotation bar (Default: 
#'                           anno.text.y.right = 
#'                           element_text(size = 7, color = "black")).
#' @param anno.ticks.y.right An \code{element_line} object to setup right Y axis
#'                           ticks on the top annotation bar (Default:
#'                           anno.ticks.y.right =
#'                           element_line(color = "black")).
#' @param lgd.text           An \code{element_text} object to setup legend
#'                           labels (Default: lgd.text =
#'                           element_text(size = 10)).
#' @param lgd.space.width    A \code{numeric} specifying the width of the legend
#'                           space (Default: lgd.space.width = 1).
#' @param show.annot         A \code{logical} to specify whether annotations
#'                           should be displayed at the top of the heatmap
#'                           (show.annot = TRUE) or not (show.annot = FALSE).
#' @param annot.size         A \code{numeric} defining the width of the
#'                           annotation bars (Default: annot.size = 1).
#' @param dend.size          A \code{numeric} defining the height of the
#'                           dendrogram made on columns
#'                           (Default: dend.col.size = 1).
#' @return A \code{grob} of the heatmap created on methylation array
#'         genotyping probes.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data.
#' library(RnBeads)
#' idat.dir <- "~/data/MethylationEPIC/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' rnb.options(identifiers.column = "barcode")
#' #Load scientific palette package ggsci.
#' library(ggsci)
#' #Plot heatmap of MethylationEPIC genotyping probes.
#' snp.htmp <- snp.heatmap(
#'   RnBSet = rnb.set,
#'   annot.grps = list("Donors" = rnb.set@pheno[, 1]),
#'   annot.pal = pal_npg("nrc", alpha = 1)(10))
#' #Save heatmap in a PDF file.
#' ggsave(
#'   filename = "heatmap.pdf", plot = snp.htmp$result.grob, device = "pdf",
#'   path = "~/")
#' @references Pageaud Y. et al., BiocompR - Advanced visualizations for data
#'             comparison.
#' @references Assenov Y. et al., Comprehensive analysis of DNA methylation data
#'             with RnBeads.

snp.heatmap <- function(
  RnBSet, dist.method = "manhattan", annot.grps = NULL, annot.pal = NULL,
  heatmap.pal = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0","#FDDBC7",
                  "#F4A582", "#D6604D", "#B2182B"),
  x.lab = "Samples", plot.title = NULL,
  htmp.text.x = element_text(
    size = 10, angle = -45, hjust = 0, vjust = 0.5, color = "black"),
  htmp.text.y.right = ggplot2::element_text(size = 7, color = "black"),
  htmp.title.y.right = ggplot2::element_text(size = 12),
  anno.text.y.right = ggplot2::element_text(size = 12, color = "black"),
  anno.ticks.y.right = ggplot2::element_line(color = "black"),
  lgd.text = ggplot2::element_text(size = 10), lgd.space.width = 1,
  show.annot = FALSE, annot.size = 1, dend.size = c(0, 2)){
  
  #Get the 65 or 59 genotyping (rs) probes
  rs.probes <- rownames(RnBSet@sites)[
    grepl(pattern = "rs", x = rownames(RnBSet@sites))]
  #Extract methylation matrix from RnBSet
  meth.mat <- RnBeads::meth(RnBSet, row.names = TRUE)
  rs.meth.mat <- meth.mat[rs.probes, ]
  #Get platform
  array.type <- methview.qc::get.platform(RnBSet = RnBSet)
  #Make plot title if none
  if(is.null(plot.title)){
    plot.title <- paste("Heatmap of", array.type, "genotyping probes")
  }
  #Set annot.grps and annot.pal if none defined
  if(is.null(annot.grps) & is.null(annot.pal)){
    annot.grps <- list(Groups = seq(ncol(rs.meth.mat)))
    annot.pal <- grDevices::rainbow(n = ncol(rs.meth.mat))
  }
  #Plot SNP CpG heatmap using genotyping probes from methylation array data
  snp.htmp <- BiocompR::gg2heatmap(
    m = rs.meth.mat, dist.method = dist.method, row.type = "genotyping probes",
    y.lab = paste(array.type, "genotyping probes"), x.lab = x.lab,
    theme_heatmap = theme(
      axis.text.y.right = htmp.text.y.right,
      axis.ticks.y.right = ggplot2::element_line(color = "black"),
      axis.text.x = htmp.text.x, axis.title.y.right = htmp.title.y.right),
    annot.grps = annot.grps, annot.pal = annot.pal, theme_annot = theme(
      axis.text.y.right = anno.text.y.right,
      axis.ticks.y.right = anno.ticks.y.right), show.annot = show.annot,
    guide_custom_bar = ggplot2::guide_colorbar(
      title = "Biallelic SNPs version", barwidth = 15, ticks.linewidth = 2,
      ticks.colour = "black", title.vjust = 0.86),
    scale_fill_grad = ggplot2::scale_fill_gradientn(
      colors = heatmap.pal, limits = c(0, 1), breaks = seq(0, 1, by = 0.5),
      labels = c("Homozygous\nV1", "Heterozygous\nV1/V2","Homozygous\nV2"),
      na.value = "black"),
    lgd.text = lgd.text, annot.size = annot.size, dend.size = dend.size,
    dendrograms = TRUE, y.axis.right = TRUE, plot.title = plot.title,
    lgd.space.width = lgd.space.width)
  return(snp.htmp)
}


#' Loads QC ggplot2 default theme.
#'
#' @return A ggplot2 \code{theme} object used as a default theme for QC
#'         plots.
#' @author Yoann Pageaud.
#' @export
#' @examples theme_qc <- methview.qc::load.metharray.QC.theme()
#' @references
#' @keywords internal

load.metharray.QC.theme <- function(){
  #Create plot theme
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 10, color = "black"),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.title.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(color = "grey"),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "grey"),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", fill = NA),
    plot.margin = grid::unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
}


#' Plots fluorescence intensities barplots for a single QC methylation array probe.
#'
#' @param array.type A \code{character} string specifying the type of array used
#'                   to generate the data (Default: array.type = "HM450K";
#'                   Supported: array.type = c("HM450K", "EPIC")).
#' @param probe.ID   A \code{character} string specifying the ID of a
#'                   methylation array quality control probe.
#' @param QC.data    A \code{data.table} list matching QC metadata with green
#'                   channel and red channel intensities, obtained with the
#'                   function \link{merge.QC.intensities.and.meta}.
#' @param DT.QC.meta A \code{data.table} with methylation array quality control
#'                   metadata obtained with the function
#'                   \link{load.metharray.QC.meta}.
#' @param cohort     A \code{character} string to specify the name of the cohort
#'                   to be displayed as part of the plot title
#'                   (Default: cohort = "RnBSet").
#' @return A \code{gtable} barplot of the QC probe fluorescence
#'         intensities.
#' @author Yoann Pageaud.
#' @export plot.array.QC.probe
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
#' #Draw probe specific QC plot for QC probe "21630339"
#' probe.plot <- plot.array.QC.probe(
#'   array.type = "EPIC", probe.ID = "21630339", QC.data = dt.mrg,
#'   DT.QC.meta = dt.meta)
#' #Save plot in a PDF file
#' ggsave(filename = "QCprobe_21630339.pdf", plot = probe.plot, device = "pdf",
#'        path = "~/")

plot.array.QC.probe <- function(
  array.type = "HM450K", probe.ID, QC.data, DT.QC.meta, cohort = "RnBSet"){
  #Melt Cy3 & Cy5 data.tables
  DT.probe.Cy3 <- data.table::melt.data.table(
    data = QC.data$`Cy3 - Electric Lime Green`[QC.probe.IDs == probe.ID],
    measure.vars = colnames(QC.data$`Cy3 - Electric Lime Green`)[-c(1:10)],
    variable.name = "Samples", value.name = "Cy3 intensity")[, c(
      "Samples", "Cy3 intensity"), ]
  DT.probe.Cy5 <- data.table::melt.data.table(data = QC.data$`Cy5 - Dark Red`[
    QC.probe.IDs == probe.ID],
    measure.vars = colnames(QC.data$`Cy5 - Dark Red`)[-c(1:10)],
    variable.name = "Samples", value.name = "Cy5 intensity")[, c(
      "Samples", "Cy5 intensity"), ]
  #Load metharray Quality Control theme
  theme_qc <- methview.qc::load.metharray.QC.theme()
  #Check if any intensity equals to 0
  if(nrow(DT.probe.Cy3[`Cy3 intensity` == 0]) +
     nrow(DT.probe.Cy5[`Cy5 intensity` == 0]) > 0){
    warning(paste("Intensity values missing for probe", probe.ID, paste(
      ":", DT.QC.meta[ID == probe.ID]$Target,
      DT.QC.meta[ID == probe.ID]$Description,
      DT.QC.meta[ID == probe.ID]$Index)))
  }
  #Check if none of the values are missing
  if(!(all(is.na(DT.probe.Cy3$`Cy3 intensity`)) &
       all(is.na(DT.probe.Cy5$`Cy5 intensity`)))){
    
    #Plot Barplot on Cy3 intensity for probe
    cy3plot <- ggplot2::ggplot(
      data = DT.probe.Cy3,
      mapping = ggplot2::aes(x = Samples, y = `Cy3 intensity`)) +
      ggplot2::geom_bar(stat = "identity", fill = "#00ff00", color = "black") +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      theme_qc
    #Plot Barplot on Cy5 intensity for probe
    cy5plot <- ggplot2::ggplot(
      data = DT.probe.Cy5,
      mapping = ggplot2::aes(x = Samples, y = `Cy5 intensity`)) +
      ggplot2::geom_bar(stat = "identity", fill = "#ff0000", color = "black") +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      theme_qc
    
    #Plot Cy5/Cy3 ratio barplot
    DT.probe.ratio <- data.table::merge.data.table(
      x = DT.probe.Cy3, y = DT.probe.Cy5, by = "Samples", all = TRUE)
    #Compute intensity ratio values.
    DT.probe.ratio <- methview.qc::compute.intensity.ratio(
      DT.probe.ratio = DT.probe.ratio)
    
    #Create DT.expected.intensity
    DT.expected.intensity <- methview.qc::get.expected.intensity(
      DT.QC.meta = DT.QC.meta, probe.id = probe.ID,
      channel.names = names(QC.data))
    
    #Barplot filled with ratio colors
    ratioplot <- ggplot2::ggplot(
      data = DT.probe.ratio, mapping = ggplot2::aes(
        x = Samples, y = `Intensity ratio`, fill = color.ratio)) +
      ggplot2::geom_bar(stat = "identity", color = "black") +
      ggplot2::scale_fill_manual(
        values = sort(unique(DT.probe.ratio$color.ratio))) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      theme_qc + ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 10, angle = -45, hjust = 0,
                                            vjust = 0.5, color = "black"),
        axis.ticks.x = ggplot2::element_line(color = "black"),
        axis.title.x = ggplot2::element_text(size = 14))
    
    #Convert ggplots in grobs
    cy3grob <- ggplot2::ggplotGrob(cy3plot)
    cy5grob <- ggplot2::ggplotGrob(cy5plot)
    ratiogrob <- ggplot2::ggplotGrob(ratioplot)
    #Resize based on widths
    ls.qc.grobs <- BiocompR::resize.grobs(ls.grobs = list(
      'cy3grob'= cy3grob, 'cy5grob' = cy5grob, 'ratiogrob' = ratiogrob),
      dimensions = 'widths', start.unit = 3, end.unit = 5)
    
    #Final plot
    qc.plot <- gridExtra::arrangeGrob(top = grid::textGrob(paste(
      cohort, "-", array.type, "Quality control intensities for",
      DT.QC.meta[ID == probe.ID]$Target, "probe",
      DT.QC.meta[ID == probe.ID]$Description, DT.QC.meta[ID == probe.ID]$Index,
      paste0("(ID = ", DT.QC.meta[ID == probe.ID]$ID, ")"))),
      grobs = list(
        grid::textGrob("Measured intensities"),
        grid::textGrob("Expected intensity"),
        ls.qc.grobs$cy3grob, grid::textGrob(label = DT.expected.intensity[
          Channel == "Cy3 - Electric Lime Green"]$`Expected intensity`),
        ls.qc.grobs$cy5grob, grid::textGrob(label = DT.expected.intensity[
          Channel == "Cy5 - Dark Red"]$`Expected intensity`),
        ls.qc.grobs$ratiogrob),
      nrow = 4, ncol = 2, heights = c(0.2, 1, 1, 2),
      widths = c(1, 8/ncol(QC.data$`Cy5 - Dark Red`[, -c(1:10), ])))
    #Return final plot
    gridExtra::grid.arrange(qc.plot)
    return(qc.plot)
  }
}

#' Plots samples fluorescence intensities distribution for QC probes of a
#' specific target type.
#'
#' @param array.type A \code{character} string specifying the type of array used
#'                   to generate the data (Default: array.type = "HM450K";
#'                   Supported: array.type = c("HM450K", "EPIC")).
#' @param target   A \code{character} string specifying the QC target type. Each
#'                 'target' matches a specific step in Illumina array
#'                 methods. Supported values:
#'                 target = c("Bisulfite Conversion I",
#'                 "Bisulfite Conversion II", "Extension", "Hybridization",
#'                 "Negative", "Non-polymorphic", "Norm A", "Norm C", "Norm G",
#'                 "Norm T", "Specificity I", "Specificity II", "Staining",
#'                 "Target Removal").
#' @param QC.data  A \code{data.table} list matching QC metadata with green
#'                 channel and red channel intensities, obtained with the
#'                 function \link{merge.QC.intensities.and.meta}.
#' @param DT.QC.meta A \code{data.table} with methylation array quality control metadata
#'                   obtained with the function \link{load.metharray.QC.meta}.
#' @param cohort   A \code{character} string to specify the name of the cohort
#'                 to be displayed as part of the plot title
#'                 (Default: cohort = "RnBSet").
#' @param ncores   An \code{integer} to specify the number of cores/threads to
#'                 be used to parallel-compute probes intensities.
#' @return A \code{gg} plot. The resulting plot depends of 2 parameters:
#'         \itemize{
#'          \item{The target type:
#'           \itemize{
#'            \item{If target = "Negative", the result is is a list of length 4,
#'                  containing 4 plots depicting fluorescence intensities of
#'                  negative probes in the red and green channels.}
#'            \item{If the target is of any other type, only 1 plot will be
#'                  returned. The resulting plot may still vary from one target
#'                  type to another.}
#'           }
#'          }
#'          \item{The number of samples considered:
#'           \itemize{
#'            \item{If 'QC.data' contains less than 5 samples, probes
#'                  fluorescence intensities will be displayed as dots, 1 dot
#'                  per sample.}
#'            \item{If 'QC.data' contains between 5 and 29 samples, probes
#'                  fluorescence intensities will be displayed as boxplots.}
#'            \item{If 'QC.data' contains 30 samples or more, probes
#'                  fluorescence intensities will be displayed as violin plots.}
#'           }
#'          }
#'         }
#' @author Yoann Pageaud.
#' @export plot.array.QC.target
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
#' #Draw target specific QC plot for "Staining" QC probes
#' target.plot <- plot.array.QC.target(
#'   array.type = "EPIC", target = "Staining", QC.data = dt.mrg,
#'   DT.QC.meta = dt.meta, ncores = 2)
#' #Save plot in a PDF file
#' ggsave(filename = "QC_staining.pdf", plot = target.plot, device = "pdf",
#'        path = "~/")

plot.array.QC.target <- function(
  array.type = "HM450K", target, QC.data, DT.QC.meta, cohort = "RnBSet",
  ncores = 1){
  #Load metharray Quality Control theme
  theme_qc <- methview.qc::load.metharray.QC.theme()
  #Update target metadata
  DT.target <- update.target.meta(
    QC.data = QC.data, target = target, ncores = ncores)
  # #Create QC boxplots for all probes target types
  # DT.target.Cy3 <- data.table::melt.data.table(
  #   data = QC.data$`Cy3 - Electric Lime Green`[Target == target],
  #   measure.vars = colnames(QC.data$`Cy3 - Electric Lime Green`)[-c(1:10)],
  #   variable.name = "Samples", value.name = "Cy3 intensity")
  # DT.target.Cy5 <- data.table::melt.data.table(
  #   data = QC.data$`Cy5 - Dark Red`[Target == target],
  #   measure.vars = colnames(QC.data$`Cy5 - Dark Red`)[-c(1:10)],
  #   variable.name = "Samples", value.name = "Cy5 intensity")
  # #Rbind data.tables
  # ls.dt.target <- list(DT.target.Cy3, DT.target.Cy5)
  # names(ls.dt.target) <- names(QC.data)
  # DT.target <- data.table::rbindlist(
  #   l = ls.dt.target, idcol = "Cyanine", use.names = FALSE)
  # #Check expected intensities for each probes
  # ls.exp.intens <- mclapply(
  #   X = unique(DT.target$QC.probe.IDs), mc.cores = ncores, FUN = function(i){
  #     methview.qc::get.expected.intensity(
  #       DT.QC.meta = DT.QC.meta, probe.id = i, channel.names = names(QC.data))
  #   })
  # names(ls.exp.intens) <- unique(DT.target$QC.probe.IDs)
  # DT.exp.intens <- data.table::rbindlist(l = ls.exp.intens, idcol = "Probe.ID")
  # #Modify DT.target with expected intensities
  # invisible(lapply(X = seq(nrow(DT.exp.intens)), FUN = function(i){
  #   DT.target[QC.probe.IDs == DT.exp.intens[i,]$Probe.ID &
  #               Cyanine == DT.exp.intens[i,]$Channel,
  #             `Expected Intensity` := DT.exp.intens[i,]$`Expected intensity`]
  # }))
  # #Change order of levels in expected intensity
  # DT.target[, `Expected Intensity` := factor(
  #   `Expected Intensity`, levels = levels(`Expected Intensity`)[c(2, 4, 3, 1)])]
  
  #Plot intensities for probes by target type
  if(target %in% c(
    "Bisulfite Conversion I", "Bisulfite Conversion II", "Extension",
    "Hybridization", "Non-polymorphic", "Specificity I", "Specificity II",
    "Staining", "Target Removal")){
    
    #Make X-Axis labels
    probe.labels <- paste(
      paste(unique(DT.target, by = "QC.probe.IDs")$Description,
            unique(DT.target, by = "QC.probe.IDs")$Index, sep = "."),
      "\n(ID = ", unique(DT.target, by = "QC.probe.IDs")$QC.probe.IDs, ")",
      sep = "")
    names(probe.labels) <- unique(DT.target$QC.probe.IDs)
    
    #Plot target
    target.plot <- ggplot2::ggplot(
      data = DT.target, mapping = ggplot2::aes(
        x = Cyanine, y = `Cy3 intensity`, fill = Cyanine)) +
      theme_qc + ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 13),
        legend.text = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 13),
        strip.text.x = ggplot2::element_text(
          size = 10, angle = 90, color = "black"),
        strip.text.y = ggplot2::element_text(
          size = 13, angle = 0, color = "black"),
        strip.background = ggplot2::element_rect(
          color = "black", size = 0.5, fill = "white"),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.spacing.x = grid::unit(0, "cm")) +
      ggplot2::labs(
        x = paste(array.type, "probe IDs"), y = "Fluorescence intensity",
        fill = "Channels", title = paste(
          cohort, "-", array.type, "Quality control intensities for", target,
          "probes")) +
      ggplot2::facet_grid(`Expected Intensity` ~ QC.probe.IDs, scales = "free",
                          labeller = ggplot2::labeller(.cols = probe.labels)) +
      ggplot2::scale_y_continuous(
        limits = c(0, max(DT.target$`Cy3 intensity`))) +
      ggplot2::scale_x_discrete(
        breaks = unique(DT.target$QC.probe.IDs), labels = probe.labels) +
      ggplot2::scale_fill_manual(values = c("#00ff00", "#ff0000"),
                                 labels = c("Green channel", "Red channel"))
    
    #If more than or equal to 30 samples draw boxplots with violins
    if(nrow(unique(DT.target, by = "Samples")) >= 30){
      target.plot <- target.plot + ggplot2::geom_violin() +
        ggplot2::geom_boxplot(fill = "white", width = 0.1)
    } else if(nrow(unique(DT.target, by = "Samples")) >= 5){
      #If more than or equal to 5 samples draw boxplots alone
      target.plot <- target.plot + ggplot2::geom_boxplot()
    } else {
      #If less than 5 samples draw dotplots
      target.plot <- target.plot + ggplot2::geom_point(
        shape = 21, size = 3, mapping = ggplot2::aes(fill = DT.target$Cyanine))
    }
    
    #Set angle for X strip labels
    if(array.type == "HM450K"){
      if(target %in% c("Bisulfite Conversion II", "Extension", "Hybridization",
                       "Non-polymorphic", "Specificity II", "Staining",
                       "Target Removal")){
        target.plot <- target.plot + ggplot2::theme(
          strip.text.x = ggplot2::element_text(
            size = 10, angle = 0, color = "black"))
        width.target.plt <- 25
      } else if(target %in% c("Bisulfite Conversion I", "Specificity I")){
        target.plot <- target.plot + ggplot2::theme(
          strip.text.x = ggplot2::element_text(
            size = 10, angle = 90, color = "black"))
        width.target.plt <- 25
      } else { stop("Probe type not supported.") }
    } else if(array.type == "EPIC"){
      if(target %in% c("Bisulfite Conversion II", "Extension", "Hybridization",
                       "Specificity II", "Staining", "Target Removal")){
        target.plot <- target.plot + ggplot2::theme(
          strip.text.x = ggplot2::element_text(
            size = 10, angle = 0, color = "black"))
        width.target.plt <- 25
      } else if(target %in% c("Bisulfite Conversion I", "Specificity I",
                              "Non-polymorphic")){
        target.plot <- target.plot + ggplot2::theme(
          strip.text.x = ggplot2::element_text(
            size = 10, angle = 90, color = "black"))
        width.target.plt <- 25
      } else { stop("Probe type not supported.") }
    }
    
    #If more than or equal to 30 samples draw boxplots with violins
    if(nrow(unique(DT.target, by = "Samples")) >= 30){
      target.plot <- target.plot + ggplot2::geom_violin() +
        ggplot2::geom_boxplot(fill = "white", width = 0.1)
    } else if(nrow(unique(DT.target, by = "Samples")) >= 5){
      #If more than or equal to 5 samples draw boxplots alone
      target.plot <- target.plot + ggplot2::geom_boxplot()
    } else {
      #If less than 5 samples draw dotplots
      target.plot <- target.plot + ggplot2::geom_point(
        shape = 21, size = 3, mapping = ggplot2::aes(fill = DT.target$Cyanine))
    }
    #Plot and return final target plot
    target.plot
    return(target.plot)
    
  } else { #Plot Norm & Negative plots
    
    if(target == "Negative") {
      #Make Y-Axis Strip labels
      strip.labels <- c("Green channel", "Red channel")
      names(strip.labels) <- unique(DT.target$Cyanine)
      
      #Calculate ranges of negative probes
      if(array.type == "HM450K"){
        neg.target.ranges <- lapply(
          X = split(x = sort(unique(DT.target, by = "QC.probe.IDs")$Index),
                    f = rep(c(1:4), times = c(154, 154, 154, 152))),
          FUN = function(i){ c(min(i), max(i)) })
      } else if(array.type == "EPIC"){
        neg.target.ranges <- lapply(
          X = split(x = sort(unique(DT.target, by = "QC.probe.IDs")$Index),
                    f = rep(c(1:3), times = c(137, 137, 137))),
          FUN = function(i){ c(min(i), max(i)) })
      } else {
        stop(
          "Unknown 'array.type'. Supported array.type are 'HM450K' & 'EPIC'.")
      }
      
      #Negative Plot
      ls.neg.plot <- lapply(X = neg.target.ranges, FUN = function(i){
        negative.plot <- ggplot2::ggplot(
          data = DT.target[Index >= i[1] & Index <= i[2]],
          mapping = ggplot2::aes(
            x = QC.probe.IDs, y = `Cy3 intensity`, fill = Cyanine)) +
          theme_qc + ggplot2::theme(
            legend.position = "none",
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_text(size = 13),
            strip.text.y = ggplot2::element_text(
              size = 13, angle = 90, color = "black"),
            strip.background = ggplot2::element_rect(
              color = "black", size = 0.5, fill = "white"),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.spacing.x = grid::unit(0, "cm")) +
          ggplot2::labs(
            x = paste(array.type, "negative probes"),
            y = "Log-scaled fluorescence intensity", title = paste(
              cohort, "-", array.type, "Quality control intensities for",
              target, "probes", i[1], "to", i[2])) +
          ggplot2::facet_grid(
            Cyanine ~., scales = "free",
            labeller = ggplot2::labeller(.rows = strip.labels)) +
          ggplot2::scale_y_log10(limits = c(
            min(DT.target$`Cy3 intensity`[DT.target$`Cy3 intensity` > 1]),
            max(DT.target$`Cy3 intensity`))) +
          ggplot2::scale_fill_manual(values = c("#00ff00", "#ff0000"))
        
        if(nrow(unique(DT.target, by = "Samples")) >= 5){
          #If more than or equal to 5 samples draw boxplots alone
          negative.plot <- negative.plot + ggplot2::geom_boxplot()
        } else {
          #If less than 5 samples draw dotplots
          negative.plot <- negative.plot +
            ggplot2::geom_point(shape = 21, size = 1, mapping = ggplot2::aes(
              fill = DT.target[Index >= i[1] & Index <= i[2]]$Cyanine))
        }
        #Plot final negative plot
        negative.plot
      })
      #Add ranges as names to ls.neg.plot
      names(ls.neg.plot) <- unlist(lapply(X = lapply(
        X = neg.target.ranges, FUN = as.character), FUN = paste0, collapse="-"))
      #Return list of negative plots
      return(ls.neg.plot)
    } else { #Target is Norm
      
      #Make Y-Axis Strip labels
      if(target %in% c("Norm C", "Norm G")){
        strip.labels <- c("Green channel\n(High)", "Red channel\n(Background)")
        names(strip.labels) <- unique(DT.target$Cyanine)
      } else {
        strip.labels <- c("Green channel\n(Background)", "Red channel\n(High)")
        names(strip.labels) <- unique(DT.target$Cyanine)
      }
      #Make X-Axis labels
      probe.labels <- paste(
        unique(DT.target, by = "QC.probe.IDs")$Description,
        unique(DT.target, by = "QC.probe.IDs")$Index, sep = ".")
      
      #Norm Plot
      norm.plot <- ggplot2::ggplot(
        data = DT.target,
        mapping = ggplot2::aes(
          x = QC.probe.IDs, y = `Cy3 intensity`, fill = Cyanine)) +
        theme_qc + ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(
            size = 10, hjust = 1, angle = 90, color = "black"),
          # axis.ticks.x = element_blank(),
          axis.title.x = ggplot2::element_text(size = 13),
          strip.text.y = ggplot2::element_text(
            size = 13, angle = 90, color = "black"),
          strip.background = ggplot2::element_rect(
            color = "black", size = 0.5, fill = "white"),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.spacing.x = grid::unit(0, "cm")) +
        ggplot2::labs(
          x = paste(array.type, "probe IDs"), y = "Fluorescence intensity",
          title = paste(cohort, "-", array.type,
                        "Quality control intensities for", target, "probes")) +
        ggplot2::facet_grid(
          Cyanine ~., scales = "free",
          labeller = ggplot2::labeller(.rows = strip.labels)) +
        ggplot2::scale_y_continuous(
          limits = c(0, max(DT.target$`Cy3 intensity`))) +
        ggplot2::scale_x_discrete(
          breaks = unique(DT.target$QC.probe.IDs), labels = probe.labels) +
        ggplot2::scale_fill_manual(values = c("#00ff00", "#ff0000"))
      
      if(nrow(unique(DT.target, by = "Samples")) >= 5){
        #If more than or equal to 5 samples draw boxplots alone
        norm.plot <- norm.plot + ggplot2::geom_boxplot()
      } else {
        #If less than 5 samples draw dotplots
        norm.plot <- norm.plot +
          ggplot2::geom_point(shape = 21, size = 3, mapping = ggplot2::aes(
            fill = DT.target$Cyanine))
      }
      #Plot and return norm plot
      norm.plot
      return(norm.plot)
    }
  }
}

#' Draws a PCA biplots on methylation array quality control targets
#' 
#' @param RnBSet           A \code{RnBSet} basic object for storing methylation
#'                         array data and experimental quality information
#'                         (Bisulfite data not supported).
#'                         \itemize{
#'                          \item{For more information about RnBSet object read
#'                          \link[RnBeads]{RnBSet-class}.}
#'                          \item{To create an RnBSet object run
#'                          \link[RnBeads]{rnb.execute.import}.}
#'                          \item{For additionnal options to import methylation
#'                          array data in the RnBSet see options available in
#'                          \link[RnBeads]{rnb.options}.}
#'                         }
#' @param PCx              An \code{integer} matching the principal component
#'                         values to display on X-axis.
#' @param PCy              An \code{integer} matching the principal component
#'                         values to display on Y-axis.
#' @param point.size       A \code{double} specifying the size of points.
#' @param loadings         A \code{logical} specifying whether the loadings
#'                         should be displayed (TRUE) or not (FALSE).
#' @param loadings.col     A \code{character} specifying a color to be used for
#'                         loadings.
#' @param top.load.by.quad An \code{integer} specifying the top n most important
#'                         loadings to be displayed in the four quadrants of the
#'                         biplot graph (by quadrants). This parameters allows
#'                         to display only the most important loadings, and to
#'                         hide the less important ones, to improve visibility
#'                         when there is too many of them.
#' @return A \code{gg} object of targets PCA biplot.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' library(RnBeads)
#' idat.dir <- "~/data/MethylationEPIC/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' rnb.options(identifiers.column = "barcode")
#' #Draw samples biplot on quality control data
#' target.biplot(RnBSet = rnb.set)
#' @references Pageaud Y. et al., BiocompR - Advanced visualizations for data
#'             comparison.

target.biplot <- function(
  RnBSet, PCx = 1, PCy = 2, point.size = 3, loadings = TRUE,
  loadings.col = "blue", top.load.by.quad = NULL){
  if(methview.qc::get.platform(RnBSet = RnBSet) == "MethylationEPIC"){
    DT.QC.meta <- methview.qc::load.metharray.QC.meta(
      array.meta = "controlsEPIC")
  } else if(methview.qc::get.platform(RnBSet = RnBSet) == "HM450K"){
    DT.QC.meta <- methview.qc::load.metharray.QC.meta(
      array.meta = "controls450")
  }
  #Merge Red and Green intensities matrices with QC probes metadata
  QC.data <- methview.qc::merge.QC.intensities.and.meta(
    RnBSet = RnBSet, DT.QC.meta = DT.QC.meta)
  QC.data <- data.table::rbindlist(
    l = QC.data, use.names = TRUE, idcol = "Channel")
  QC.data[, Target := as.factor(Target)]
  #Compute PCA
  pca_res <- stats::prcomp(x = QC.data[, -c(1:11), ], scale. = FALSE)
  #Plot target biplot
  if(is.null(top.load.by.quad)){
    target <- BiocompR::ggbipca(
      prcomp.res = pca_res, data = QC.data[, c(1:11),], PCx = PCx, PCy = PCy,
      color.data = "Channel", shape.data = "Target", point.size = point.size,
      loadings = loadings, loadings.col = loadings.col)
  } else {
    target <- BiocompR::ggbipca(
      prcomp.res = pca_res, data = QC.data[, c(1:11),], PCx = PCx, PCy = PCy,
      color.data = "Channel", shape.data = "Target", point.size = point.size,
      loadings = loadings, loadings.col = loadings.col,
      top.load.by.quad = top.load.by.quad)
  }
  target <- target + ggplot2::scale_color_manual(values = c("#00ff00", "red"))
  if(methview.qc::get.platform(RnBSet = RnBSet) == "MethylationEPIC"){
    target <- target + ggplot2::scale_shape_manual(
      values = c(66, 66, 69, 72, 25, 78, 88, 88, 88, 88, 82, 83, 83, 8, 84))
  } else if(methview.qc::get.platform(RnBSet = RnBSet) == "HM450K"){
    target <- target + ggplot2::scale_shape_manual(
      values = c(66, 66, 69, 72, 25, 78, 88, 88, 88, 88, 83, 83, 8, 84))
  } else {
    stop("Unknown 'array.type'. Supported array.type are 'HM450K' & 'EPIC'.")
  }
  return(target)
}


#' Draws a customizable PCA biplots on samples methylation array QC data
#' 
#' @param RnBSet           A \code{RnBSet} basic object for storing methylation
#'                         array data and experimental quality information
#'                         (Bisulfite data not supported).
#'                         \itemize{
#'                          \item{For more information about RnBSet object read
#'                          \link[RnBeads]{RnBSet-class}.}
#'                          \item{To create an RnBSet object run
#'                          \link[RnBeads]{rnb.execute.import}.}
#'                          \item{For additionnal options to import methylation
#'                          array data in the RnBSet see options available in
#'                          \link[RnBeads]{rnb.options}.}
#'                         }
#' @param PCx              An \code{integer} matching the principal component
#'                         values to display on X-axis.
#' @param PCy              An \code{integer} matching the principal component
#'                         values to display on Y-axis.
#' @param point.size       A \code{double} specifying the size of points.
#' @param loadings         A \code{logical} specifying whether the loadings
#'                         should be displayed (TRUE) or not (FALSE).
#' @param loadings.col     A \code{character} specifying a color to be used for
#'                         loadings.
#' @param top.load.by.quad An \code{integer} specifying the top n most important
#'                         loadings to be displayed in the four quadrants of the
#'                         biplot graph (by quadrants). This parameters allows
#'                         to display only the most important loadings, and to
#'                         hide the less important ones, to improve visibility
#'                         when there is too many of them
#'                         (Default: top.load.by.quad = 5).
#' @param color.data       A \code{character} specifying the column name in
#'                         'data' to be used to map colors to points. You can
#'                         specify your own custom palette of colors using the
#'                         'scale_color_manual()' function. For more information
#'                         about how to use it see
#'                         \link[ggplot2]{scale_color_manual}.
#' @param shape.data       A \code{character} specifying the column name in
#'                         'data' to be used to map shapes to points. You can
#'                         specify your own custom set of point shapes using the
#'                         'scale_shape_manual()' function. For more information
#'                         about how to use it see
#'                         \link[ggplot2]{scale_shape_manual}.
#' @return A customizable \code{gg} object of samples PCA biplot.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' library(RnBeads)
#' idat.dir <- "~/data/MethylationEPIC/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' rnb.options(identifiers.column = "barcode")
#' #Draw samples biplot on quality control data
#' sampleQC.biplot(RnBSet = rnb.set)
#' @references Pageaud Y. et al., BiocompR - Advanced visualizations for data
#'             comparison.

sampleQC.biplot <- function(
  RnBSet, PCx = 1, PCy = 2, loadings = TRUE, loadings.col = "blue",
  point.size = 2.5, top.load.by.quad = 5, color.data = "ID", shape.data = NULL){
  if(methview.qc::get.platform(RnBSet = RnBSet) == "MethylationEPIC"){
    DT.QC.meta <- methview.qc::load.metharray.QC.meta(
      array.meta = "controlsEPIC")
  } else if(methview.qc::get.platform(RnBSet = RnBSet) == "HM450K"){
    DT.QC.meta <- methview.qc::load.metharray.QC.meta(
      array.meta = "controls450")
  }
  #Merge Red and Green intensities matrices with QC probes metadata
  QC.data <- methview.qc::merge.QC.intensities.and.meta(
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
  if(is.null(top.load.by.quad)){
    if(is.null(shape.data)){
      sample.biplot <- BiocompR::ggbipca(
        prcomp.res = pca_t.res, data = t.QC.dt[, 1:ncol(RnBSet@pheno)],
        PCx = PCx, PCy = PCy, loadings = loadings, loadings.col = loadings.col,
        point.size = point.size, color.data = color.data)
    } else {
      sample.biplot <- BiocompR::ggbipca(
        prcomp.res = pca_t.res, data = t.QC.dt[, 1:ncol(RnBSet@pheno)],
        PCx = PCx, PCy = PCy, loadings = loadings, loadings.col = loadings.col,
        point.size = point.size, color.data = color.data,
        shape.data = shape.data)
    }
  } else {
    if(is.null(shape.data)){
      sample.biplot <- BiocompR::ggbipca(
        prcomp.res = pca_t.res, data = t.QC.dt[, 1:ncol(RnBSet@pheno)],
        PCx = PCx, PCy = PCy, loadings = loadings, loadings.col = loadings.col,
        top.load.by.quad = top.load.by.quad, point.size = point.size,
        color.data = color.data)
    } else {
      sample.biplot <- BiocompR::ggbipca(
        prcomp.res = pca_t.res, data = t.QC.dt[, 1:ncol(RnBSet@pheno)],
        PCx = PCx, PCy = PCy, loadings = loadings, loadings.col = loadings.col,
        top.load.by.quad = top.load.by.quad, point.size = point.size,
        color.data = color.data, shape.data = shape.data)
    }
  }
  return(sample.biplot)
}


#' Plots FFPE negative control probe fluorescence intensities barplots (WARNING: Experimental).
#' 
#' @param RnBSet A \code{RnBSet} basic object for storing methylation array data
#'               and experimental quality information (HM450K and Bisulfite data
#'               not supported).
#'               \itemize{
#'                \item{For more information about RnBSet object read
#'                \link[RnBeads]{RnBSet-class}.}
#'                \item{To create an RnBSet object run
#'                \link[RnBeads]{rnb.execute.import}.}
#'                \item{For additionnal options to import methylation array data
#'                in the RnBSet see options available in
#'                \link[RnBeads]{rnb.options}.}
#'               }
#' @param cohort A \code{character} string to specify the name of the cohort to
#'               be displayed as part of the plot title
#'               (Default: cohort = "RnBSet").
#' @return A \code{gtable} barplot of the QC probe fluorescence intensities.
#' @details
#' In MethylationEPIC array we identified the Negative quality control probe
#' "Negative 1720" as sensitive to Formaldehyd Fixed Paraffin Embed (FFPE)
#' sample treatment.\cr This probe (Negative 1720; ID = 36729435) displays
#' unusual high fluorescence intensities in FFPE samples in both, green and red
#' channels, more pronounced in the green channel.\cr We suspect that a similar
#' Negative quality control probe also exists in HM450K array, but have not yet
#' identified one. 
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' library(RnBeads)
#' idat.dir <- "~/data/MethylationEPIC/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' rnb.options(identifiers.column = "barcode")
#' # Plot FFPE negative control probe
#' neg.FFPE <- methview.qc:::plot.negative.FFPE(RnBSet = rnb.set)
#' #Save plot
#' ggsave(filename = "FFPE_negative_control_probe_MethylationEPIC.pdf",
#'        plot = neg.FFPE, device = "pdf", width = 11, height = 7.3,
#'        path = "~/")
#' @keywords internal

plot.negative.FFPE <- function(RnBSet, cohort = "RnBSet"){
  #Set array type
  if(get.platform(RnBSet = RnBSet) == "MethylationEPIC"){
    #Create the data.table with quality control metadata
    DT.QC.meta <- load.metharray.QC.meta(array.meta = "controlsEPIC")
    # Merge red and green channels intensities with QC metadata
    QC.data <- methview.qc::merge.QC.intensities.and.meta(
      RnBSet = RnBSet, DT.QC.meta = DT.QC.meta)
    # Plot FFPE negative control prob 
    methview.qc::plot.array.QC.probe(
      array.type = "EPIC", probe.ID = "36729435", QC.data = QC.data,
      DT.QC.meta = DT.QC.meta, cohort = cohort)
  } else if(get.platform(RnBSet = RnBSet) == "HM450K"){
    stop("No FFPE negative control probe identified yet in HM450K platform.")
  } else {
    stop(paste(
      "RnBSet platform not supported.",
      "Supported platforms are HM450K and MethylationEPIC.",
      "Please contact developper to request support for your methylation data.")
    )
  }
}

#' Draws and saves all quality control plots available in methview.qc
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
#' @param cohort     A \code{character} string to specify the name of the cohort
#'                   to be displayed as part of the plot title
#'                   (Default: cohort = "RnBSet").
#' @param save.dir   A \code{character} string to specify the path where the
#'                   quality control plots should be saved.
#' @param ncores     An \code{integer} to specify the number of cores/threads to
#'                   be used to parallel-compute probes intensities.
#' @param include.gp A \code{logical} to specify whether the genotyping probes
#'                   heatmap should be plotted too (include.gp = TRUE) or not
#'                   (include.gp = TRUE).\cr If you wish to customize your
#'                   genotyping probes heatmap, include.gp must be set to FALSE.
#'                   For more information see the details section.
#' @param include.ffpe A \code{logical} to specify whether the FFPE negative
#'                     control probe should be plotted a second time separately
#'                     from the rest of the Negative QC probes (Warning: this is
#'                     an experimental function; for more information about the
#'                     FFPE negative control probe see
#'                     \link{plot.negative.FFPE}).
#' @details
#' The genotyping probes heatmap produced by \link{plot.all.qc} when
#' \code{include.gp = TRUE} is an automatic, non-custom heatmap created using
#' \link{snp.heatmap}. You can customize the genotyping probes heatmap using the
#' function \link{snp.heatmap} outside \link{plot.all.qc}. Using
#' \link{snp.heatmap} you can provide custom annotations to the top annotation
#' bar, and personalized a lot more the different components of the plot.
#' @author Yoann Pageaud.
#' @export plot.all.qc
#' @export
#' @examples
#' #Create an RnBSet for MethylationEPIC data
#' library(RnBeads)
#' idat.dir <- "~/data/MethylationEPIC/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
#' rnb.options(identifiers.column = "barcode")
#' #Draw all plots from the quality control data of rnb.set
#' plot.all.qc(RnBSet = rnb.set, save.dir = "~/", ncores = 2)

plot.all.qc <- function(
  RnBSet, cohort = "RnBSet", save.dir, ncores = 1, include.gp = TRUE,
  include.ffpe = FALSE){
  #Set array type
  if(get.platform(RnBSet = RnBSet) == "MethylationEPIC"){
    array.type <- "EPIC"
    DT.QC.meta <- methview.qc::load.metharray.QC.meta(
      array.meta = "controlsEPIC")
  } else if(get.platform(RnBSet = RnBSet) == "HM450K"){
    array.type <- "HM450K"
    DT.QC.meta <- methview.qc::load.metharray.QC.meta(
      array.meta = "controls450")
  } else {
    stop(paste(
      "RnBSet platform not supported.",
      "Supported platforms are HM450K and MethylationEPIC.",
      "Please contact developper to request support for your methylation data.")
    )
  }
  
  #Merge Red and Green intensities matrices with QC probes metadata
  QC.data <- methview.qc::merge.QC.intensities.and.meta(
    RnBSet = RnBSet, DT.QC.meta = DT.QC.meta)
  
  #Create graph directory
  if(!dir.exists(save.dir)){
    dir.create(save.dir)
  }
  #Create QC_Barplots & QC_Targets subdirectories
  if(!dir.exists(file.path(save.dir, "QC_Barplots"))){
    dir.create(file.path(save.dir, "QC_Barplots"))
  }
  if(!dir.exists(file.path(save.dir, "QC_Targets"))){
    dir.create(file.path(save.dir, "QC_Targets"))
  }
  cat("Plotting...\n")
  #Loop over probes target types
  invisible(lapply(X = levels(DT.QC.meta$Target), FUN = function(target){
    cat(paste0("\t", target, "\n"))
    
    #Create cohort & target directory
    if(!dir.exists(file.path(save.dir, "QC_Barplots", cohort))){
      dir.create(path = file.path(save.dir, "QC_Barplots", cohort))
    }
    if(!dir.exists(
      file.path(save.dir, "QC_Barplots", cohort, target))){
      dir.create(path = file.path(save.dir, "QC_Barplots", cohort, target))
    }
    if(!dir.exists(file.path(save.dir, "QC_Targets", cohort))){
      dir.create(path = file.path(save.dir, "QC_Targets", cohort))
    }
    #Create QC barplots for every probe
    invisible(mclapply(
      X = DT.QC.meta[Target == target]$ID, mc.cores = ncores, FUN = function(i){
        #Make QC barplots for every probe
        qc.plot <- methview.qc::plot.array.QC.probe(
          array.type = array.type, probe.ID = i, QC.data = QC.data,
          DT.QC.meta = DT.QC.meta, cohort = cohort)
        
        #Check if none of the values are missing
        if(!(all(is.na(unlist(QC.data$`Cy3 - Electric Lime Green`[
          QC.probe.IDs == i, -c(1:10), ]))) &
          all(is.na(unlist(QC.data$`Cy5 - Dark Red`[
            QC.probe.IDs == i, -c(1:10), ])))
        )){
          if(ncol(QC.data$`Cy5 - Dark Red`[, -c(1:10), ]) > 40 &
             ncol(QC.data$`Cy5 - Dark Red`[, -c(1:10), ]) < 70){
            width.plot <- 45
          } else if(ncol(QC.data$`Cy5 - Dark Red`[, -c(1:10), ]) >= 70 &
                    ncol(QC.data$`Cy5 - Dark Red`[, -c(1:10), ]) < 110){
            width.plot <- 65
          } else if(ncol(QC.data$`Cy5 - Dark Red`[, -c(1:10), ]) >= 110){
            width.plot <- 135
          } else { width.plot <- 25 }
          #Save plot
          ggsave(filename = paste(
            "Barplot", cohort, array.type, "QC", paste0(
              gsub(pattern = " ", replacement = "_", x = paste(
                DT.QC.meta[ID == i]$Target, DT.QC.meta[ID == i]$Description,
                DT.QC.meta[ID == i]$Index)), ".pdf"), sep = "_"),
            plot = qc.plot, device = "pdf",
            path = file.path(save.dir, "QC_Barplots", cohort, target),
            width = width.plot, height = 16, units = "cm", limitsize = FALSE)
        }
      }))
    
    #Create QC plots for the target type
    target.plot <- methview.qc::plot.array.QC.target(
      array.type = array.type, target = target, QC.data = QC.data,
      DT.QC.meta = DT.QC.meta, ncores = ncores, cohort = cohort)
    
    #Plot intensities for probes by target type
    if(target %in% c(
      "Bisulfite Conversion I", "Bisulfite Conversion II", "Extension",
      "Hybridization", "Non-polymorphic", "Specificity I", "Specificity II",
      "Staining", "Target Removal")){
      
      #Save plot
      ggsave(filename = paste(
        "Target", cohort, array.type, "QC", paste0(target, ".pdf"), sep = "_"),
        plot = target.plot, device = "pdf",
        path = file.path(save.dir, "QC_Targets", cohort),
        width = 25, height = 16, units = "cm", limitsize = FALSE)
      
    } else { #Plot Norm & Negative plots
      
      if(target == "Negative") {
        #Negative Plot
        invisible(lapply(X = seq_along(target.plot), FUN = function(i){
          #Save plot
          ggsave(filename = paste(
            "Target", cohort, array.type, "QC", target,
            paste0(names(target.plot)[i], ".pdf"), sep = "_"),
            plot = target.plot[[i]], device = "pdf",
            path = file.path(save.dir, "QC_Targets", cohort),
            width = 25, height = 16, units = "cm", limitsize = FALSE)
        }))
      } else { #Target is Norm
        #Save plot
        ggsave(filename = paste(
          "Target", cohort, array.type, "QC", paste0(target, ".pdf"),
          sep = "_"), plot = target.plot, device = "pdf", path = file.path(
            save.dir, "QC_Targets", cohort), width = 25, height = 16,
          units = "cm", limitsize = FALSE)
      }
    }
  }))
  if(include.ffpe){
    cat("\tFFPE Negative control probe\n")
    # Plot FFPE negative control probe
    neg.FFPE <- methview.qc:::plot.negative.FFPE(
      RnBSet = RnBSet, cohort = cohort)
    ggsave(
      filename = "FFPE_negative_control_probe_MethylationEPIC.pdf",
      plot = neg.FFPE, device = "pdf", width = 11, height = 7.3,
      path = file.path(save.dir, "QC_Barplots", cohort))
  }
  if(include.gp){
    cat("\tGenotyping probes heatmap\n")
    #Plot genotyping probes heatmap
    snp.htmp <- snp.heatmap(
      RnBSet = RnBSet, annot.grps = list("IDs" = RnBSet@pheno[, 1]),
      annot.pal = rainbow(n = length(RnBSet@pheno[, 1])), plot.title =
        paste(cohort, "- Heatmap of", get.platform(RnBSet = RnBSet),
              "genotyping probes"))
    invisible(lapply(X = c("pdf", "png"), FUN = function(frmt){
      ggsave(
        filename = paste0("Heatmap_genotyping_probes_",
                          get.platform(RnBSet = RnBSet), ".", frmt),
        plot = snp.htmp$result.grob, device = frmt, width = 11, height = 11,
        path = save.dir)
    }))
  }
  cat("Done!\n")
}


#' Plots QC deviation heatmaps based on samples fluorescence deviation score.
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
#' @param target     A \code{character} string specifying the QC target type.
#'                   Each 'target' matches a specific step in Illumina array
#'                   methods. Supported values:
#'                   target = c("Bisulfite Conversion I",
#'                   "Bisulfite Conversion II", "Extension", "Hybridization",
#'                   "Negative", "Non-polymorphic", "Norm A", "Norm C",
#'                   "Norm G", "Norm T", "Specificity I", "Specificity II",
#'                   "Staining", "Target Removal", NULL). If target is NULL,
#'                   quality control probes from all targets will be considered
#'                   (Default: target = NULL).
#' @param samples    A \code{character} vector specifying the samples to include
#'                   for the deviation score calculation. You can catch the
#'                   sample IDs you wish to evaluate running
#'                   \code{RnBSet@pheno[,1]}. If samples is NULL, all samples in
#'                   RnBSet will be considered (Default: samples = NULL).
#' @param dist.col   A \code{character} to specify the distance methods to be
#'                   used for clustering data on columns and plot dendrogram on 
#'                   columns (Default: dist.col = "manhattan").
#' @param annot.grps A \code{list} of vectors of groups to which variables
#'                   belongs for the annotation sidebars. Vectors' lengths have
#'                   to match the number of variables.
#' @param annot.pal  A \code{vector} or a list of vectors containing colors as
#'                   characters for the annotation sidebars. The length of
#'                   vectors has to match the number of levels of vectors listed
#'                   in 'annot.grps'.
#'                   \itemize{
#'                    \item{If annot.pal is a list: its length must match the
#'                    length of the list provided to 'annot.grps'.}
#'                    \item{If annot.pal is a vector: make sure that the levels
#'                    content of annotations listed in 'annot.grps' is the same,
#'                    and that no annotation contains less or more levels than
#'                    another one in 'annot.grps'.}
#'                   }
#' @param annot.size A \code{numeric} defining the width of the annotation bars
#'                   (Default: annot.size = 1).
#' @param show.annot A \code{logical} to specify whether annotations should be
#'                   displayed at the top of the heatmap (show.annot = TRUE) or
#'                   not (show.annot = FALSE).
#' @param dend.size  A \code{numeric} defining columns dendrogram size
#'                   (Default: dend.size = 1).
#' @param lgd.width  A \code{numeric} specifying the width of the legend space
#'                   (Default: lgd.space.width = 1).
#' @param ncores     An \code{integer} specifying the number of cores or threads
#'                   to be used for parallel processing.
#' @param draw       A \code{logical} to specify whether the final heatmap
#'                   should be drawn automatically when gg2heatmap() execution
#'                   ends (draw = TRUE), or if it shouldn't (draw = FALSE).
#' @param verbose    A \code{logical} to display information about the
#'                   step-by-step processing of the data if TRUE
#'                   (Default: verbose = FALSE).
#' @return A \code{grob} of the computed heatmap. If draw = FALSE, you can draw
#'         the plot using the grid::grid.draw() function on the grob.
#' @details For more information about how the fluorescence deviation score is
#'          calculated, please refer to the details section of
#'          \link{devscore.fluo}.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' #Create an RnBSet for HM450K data
#' library(RnBeads)
#' idat.dir <- "~/data/my_idat_dir/"
#' sample.annotation <- "~/data/Annotations/sample_sheet.csv"
#' data.source <- c(idat.dir, sample.annotation)
#' rnb.set <- rnb.execute.import(
#'   data.source = data.source, data.type = "idat.dir")
#' #Draw deviation score heatmaps for green and red channels
#' dev.heatmap <- devscore.heatmap(RnBSet = rnb.set)
#' #If 'draw' is set to FALSE you can plot heatmap as following
#' grid::grid.newpage()
#' grid::grid.draw(dev.heatmap)

devscore.heatmap <- function(
  RnBSet, target = NULL, samples = NULL, dist.col = "manhattan",
  annot.grps = NULL, annot.pal = NULL, annot.size = 1,
  show.annot = FALSE, dend.size = 1, lgd.width = 1, ncores = 1,
  draw = TRUE, verbose = FALSE){
  #If no specific samples provided take them all
  if(is.null(samples)){ samples <- as.character(RnBSet@pheno[, 1]) }
  
  if(is.null(target)){
    ls.dt.target <- lapply(
      X = levels(load.metharray.QC.meta("controls450")$Target),
      FUN = function(t){ methview.qc::devscore.fluo(
        RnBSet = RnBSet, samples = samples, target = t, ncores = ncores)})
    DT.target <- data.table::rbindlist(ls.dt.target)
  } else {
    #Compute deviation score of fluorescence
    DT.target <- methview.qc::devscore.fluo(
      RnBSet = RnBSet, samples = samples, target = target, ncores = ncores)
  }
  #Create probes labels
  DT.target[, probe.labels := paste(paste(Description, Index, sep = "."),
                                    " (ID = ", QC.probe.IDs, ")", sep = "")]
  #Convert Cyanine as factor
  DT.target[, Cyanine := as.factor(Cyanine)]
  #Rename levels
  data.table::setattr(x = DT.target$Cyanine, name = "levels",
                      value = c("Green Channel", "Red Channel"))
  #Create molten data.frame for gg2heatmap input
  molten_dt <- DT.target[, c(14, 12, 13, 1:11), ]
  #Set annotation default values
  if(is.null(annot.grps)){
    annot.grps <- list("Groups" = unique(molten_dt[[2]]))
  }
  if(is.null(annot.pal)){
    annot.pal <- grDevices::rainbow(n = length(unique(molten_dt[[2]])))
  }
  #If target is NULL means all targets
  if(is.null(target)){ target <- "All" }
  #Set plot theme
  if(target %in% c(
    "Bisulfite Conversion I", "Bisulfite Conversion II", "Extension",
    "Hybridization", "Non-polymorphic", "Specificity I", "Specificity II",
    "Staining", "Target Removal")){
    plot_theme = theme(
      axis.title.y.left = ggplot2::element_text(size = 14),
      axis.text.y.left = ggplot2::element_text(size = 12, colour = "black"),
      axis.ticks.y.left = ggplot2::element_line(color = "black"),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = element_text(
        size = 8, angle = -90, hjust = 0, vjust = 0.5, colour = "black"),
      strip.text.y = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1 ))
  } else if(target %in% c("Norm A", "Norm G")){
    plot_theme = theme(
      axis.title.y.left = ggplot2::element_text(size = 14),
      axis.text.y.left = ggplot2::element_text(size = 8, colour = "black"),
      axis.ticks.y.left = ggplot2::element_line(color = "black"),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = element_text(
        size = 8, angle = -90, hjust = 0, vjust = 0.5, colour = "black"),
      strip.text.y = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1 ))
  } else if(target %in% c("Norm C", "Norm T", "Negative", "All")){
    plot_theme = theme(
      axis.title.y.left = ggplot2::element_text(size = 14),
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = element_text(
        size = 8, angle = -90, hjust = 0, vjust = 0.5, colour = "black"),
      strip.text.y = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1 ))
  }
  #Row string
  if(target == "All"){ row.string <- "QC probes measures" } else {
    row.string <- paste(target, "QC probes measures")
  }
  #Draw heatmap
  res.htmp <- BiocompR::gg2heatmap(
    m = molten_dt, dendrograms = c(FALSE, TRUE), dend.size = c(0, dend.size),
    facet = NULL, row.type = row.string,
    plot.title = "HM450K fluorescence deviation score",
    dist.method = c("none", dist.col), lgd.space.width = lgd.width,
    y.lab = paste(target, "quality control probes"), annot.grps = annot.grps,
    annot.pal = annot.pal, annot.size = annot.size, show.annot = show.annot,
    scale_fill_grad = scale_fill_gradientn(
      colors = c("mediumblue", "mediumblue", "green3", "green3", "red", "red"),
      values = rescale(c(-100, -30, -20, 20, 30, 100)), limits = c(-100, 100),
      breaks = seq(-100, 100, by = 20),
      labels = c("-100%\nor lower", paste0(seq(-80, 80, by = 20), "%"),
                 "100%\nor higher"), oob = squish, na.value = "black"),
    guide_custom_bar = guide_colorbar(
      title = "Fluorescence\ndeviation", title.vjust = 0.86,
      ticks.linewidth = 2, barwidth = 25),
    theme_heatmap = plot_theme, split.by.rows = "Cyanine", na.handle = "keep",
    ncores = ncores, verbose = verbose, draw = draw)
  #Return Grob of the final heatmap
  return(res.htmp$result.grob)
}