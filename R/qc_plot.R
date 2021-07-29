
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
#' @param axis.text.x        An \code{element_text} object to setup X axis text
#'                           (Default: axis.text.x = element_text(size = 10,
#'                           angle = -45, hjust = 0, vjust = 0.5,
#'                           color = "black")).
#' @param axis.text.y.right  An \code{element_text} object to setup right Y axis
#'                           text (Default: axis.text.y.right =
#'                           element_text(size = 7, color = "black")).
#' @param axis.title.y.right An \code{element_text} object to setup right Y axis
#'                           title
#'                           (Default: axis.title.y.right =
#'                           element_text(size = 11)).
#' @param lgd.text           An \code{element_text} object to setup legend
#'                           labels (Default: lgd.text =
#'                           element_text(size = 10)).
#' @param lgd.space.width    A \code{numeric} specifying the width of the legend
#'                           space (Default: lgd.space.width = 1).
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
#'   annot.grps = list("Donors" = rnb.set@pheno$ID),
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
  RnBSet, dist.method = "manhattan", annot.grps, annot.pal, heatmap.pal = c(
    "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0","#FDDBC7", "#F4A582", "#D6604D",
    "#B2182B"), x.lab = "Samples", plot.title = NULL,
  axis.text.x = element_text(
    size = 10, angle = -45, hjust = 0, vjust = 0.5, color = "black"),
  axis.text.y.right = element_text(size = 7, color = "black"),
  axis.title.y.right = element_text(size = 11),
  lgd.text = element_text(size = 10), lgd.space.width = 1, annot.size = 1,
  dend.size = 1){

  #Get the 65 or 59 genotyping (rs) probes
  rs.probes <- rownames(RnBSet@sites)[
    grepl(pattern = "rs", x = rownames(RnBSet@sites))]
  #Extract methylation matrix from RnBSet
  meth.mat <- RnBeads::meth(RnBSet, row.names = TRUE)
  rs.meth.mat <- meth.mat[rs.probes, ]

  if(length(rs.probes) == 59){
    array.type <- "MethylationEPIC"
    if(is.null(plot.title)){
      plot.title <- "Heatmap of MethylationEPIC genotyping probes"
    }
  } else if(length(rs.probes) == 65){
    array.type <- "HM450K"
    if(is.null(plot.title)){
      plot.title <- "Heatmap of HM450K genotyping probes"
    }
  } else {
    stop("RnBSet platform not supported. Supported platforms are HM450K and MethylationEPIC.")
  }

  #Plot SNP CpG heatmap using genotyping probes from methylation array data
  snp.heatmap <- BiocompR::gg2heatmap(
    m = rs.meth.mat, dist.method = dist.method, row.type = "genotyping probes",
    y.lab = paste(array.type, "genotyping probes"), x.lab = x.lab,
    axis.text.y.right = axis.text.y.right,
    axis.ticks.y.right = ggplot2::element_line(color = "black"),
    axis.text.x = axis.text.x, axis.title.y.right = axis.title.y.right,
    annot.grps = annot.grps, annot.pal = annot.pal,
    lgd.scale.name = "Methylation", lgd.text = lgd.text,
    annot.size = annot.size, dend.size = dend.size,
    heatmap.pal = heatmap.pal, dendrograms = TRUE, y.axis.right = TRUE,
    plot.title = plot.title, lgd.space.width = lgd.space.width)

  return(snp.heatmap)
}


#' Loads QC ggplot2 default theme.
#'
#' @return A ggplot2 \code{theme} object used as a default theme for QC
#'         plots.
#' @author Yoann Pageaud.
#' @export
#' @examples theme_qc <- HM450.QCView::load.metharray.QC.theme()
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
#' probe.plot <- HM450.QCView:::plot.array.QC.probe(
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
  theme_qc <- HM450.QCView::load.metharray.QC.theme()
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
    DT.probe.ratio <- HM450.QCView::compute.intensity.ratio(
      DT.probe.ratio = DT.probe.ratio)

    #Create DT.expected.intensity
    DT.expected.intensity <- HM450.QCView::get.expected.intensity(
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
#' @export
#' @examples
#' #' #Create an RnBSet for MethylationEPIC data
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
#' target.plot <- HM450.QCView:::plot.array.QC.target(
#'   array.type = "EPIC", target = "Staining", QC.data = dt.mrg,
#'   DT.QC.meta = dt.meta, ncores = 2)
#' #Save plot in a PDF file
#' ggsave(filename = "QC_staining.pdf", plot = target.plot, device = "pdf",
#'        path = "~/")

plot.array.QC.target <- function(
  array.type = "HM450K", target, QC.data, DT.QC.meta, cohort = "RnBSet",
  ncores = 1){
  #Create QC boxplots for all probes target types
  DT.target.Cy3 <- data.table::melt.data.table(
    data = QC.data$`Cy3 - Electric Lime Green`[Target == target],
    measure.vars = colnames(QC.data$`Cy3 - Electric Lime Green`)[-c(1:10)],
    variable.name = "Samples", value.name = "Cy3 intensity")
  DT.target.Cy5 <- data.table::melt.data.table(
    data = QC.data$`Cy5 - Dark Red`[Target == target],
    measure.vars = colnames(QC.data$`Cy5 - Dark Red`)[-c(1:10)],
    variable.name = "Samples", value.name = "Cy5 intensity")
  #Load metharray Quality Control theme
  theme_qc <- HM450.QCView::load.metharray.QC.theme()
  #Rbind data.tables
  ls.dt.target <- list(DT.target.Cy3, DT.target.Cy5)
  names(ls.dt.target) <- names(QC.data)
  DT.target <- data.table::rbindlist(
    l = ls.dt.target, idcol = "Cyanine", use.names = FALSE)

  #Check expected intensities for each probes
  ls.exp.intens <- mclapply(
    X = unique(DT.target$QC.probe.IDs), mc.cores = ncores, FUN = function(i){
      HM450.QCView::get.expected.intensity(
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
