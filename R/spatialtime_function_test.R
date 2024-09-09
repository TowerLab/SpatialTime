# libraries ----
library(Seurat)
library(tidyverse)
library(Hmisc)
library(jsonlite)
library(monocle)
library(devtools)
library(testthat)
library(usethis)
library(pkgdown)

pkgdown::build_site()
usethis::use_pkgdown_github_pages()

# Import dataset ----
fracture <- readRDS("R/fracture_original_annotation.rds")
is(fracture, "Seurat")
Idents(fracture)


#' @file - Seurat Object
#' @cluster - cluster names in Idents
#' @export.all - Export coordinates files
#' @slice.n - slice number
#' @dir.out - directory output

subsetLabels <- function(file = "", cluster = NULL, export.all = T,
                         slice.n = "slice1", dir.out = "") {

  if (file.exists(file)) {
    data <- readRDS(file)

  } else {
    stop("File do not exist.")
  }

  if (export.all) {

    # Check if its Seurat data
    if (is.null(cluster)) {
      stop("Cluster IDs were not specified.")
    }

    if (length(cluster) > 3) {
      stop("Number of clusters should be less than 3 per times.")
    }

    if (all(cluster %in% Idents(data))) {
      clstrs <- list()

      for (i in seq_along(cluster)) {
        clstrs[[i]] <- subset(data, idents = cluster[i])
      }
    }

    for (z in seq_along(clstrs)) {
      write.csv(clstrs[[z]]@images[[slice.n]]@coordinates,
                paste0(dir.out, cluster[z], "_coordinates.csv"))

      return(clstrs)

    }

  } else {

    write.csv(data@images[[slice.n]]@coordinates, paste0(dir.out, "all_coordinates.csv"))

  }
}

# if "export.all = T", then "cluster" parameter is not necessary
subsetLabels("R/fracture_original_annotation.rds", cluster = c("MarrowFibrin", "Periosteum", "Cortical"),
             slice.n = "slice1", export.all = T, dir.out = "R/excel/") # Add error handling for terms that do not exist or poorly written



################################## WORKING #####################################
subsetLabels <- function(file = "", cluster = NULL, export.all = TRUE,
                         slice.n = "slice1", dir.out = "") {

  if (!file.exists(file)) {
    stop("File does not exist.")
  }

  data <- readRDS(file)

  if (export.all) {
    if (is.null(cluster)) {
      stop("Cluster IDs were not specified.")
    }

    if (length(cluster) > 3) {
      stop("Number of clusters should be less than or equal to 3.")
    }

    if (!all(cluster %in% Idents(data))) {
      stop("One or more specified clusters do not exist in the data.")
    }

    clstrs <- list()

    for (i in seq_along(cluster)) {
      clstrs[[i]] <- subset(data, idents = cluster[i])
      write.csv(clstrs[[i]]@images[[slice.n]]@coordinates,
                file.path(dir.out, paste0(cluster[i], "_coordinates.csv")))
    }

    return(clstrs)

  } else {
    write.csv(data@images[[slice.n]]@coordinates,
              file.path(dir.out, "all_coordinates.csv"))
  }
}

subsetLabels("R/fracture_original_annotation.rds", cluster = c("MarrowFibrin", "Periosteum", "Cortical"),
             slice.n = "slice1", export.all = T, dir.out = "R/excel/")

#' @file - Fiji file
#' @factor - scale factor
#' @colors - Line color
#' @tissue - Tissue ID

SpatialCalc <- function(file = "", factor = 1, colors = NULL, tissue = NULL) {

  x <- read_csv(file)

  for (i in 1:length(colors)) {
    x <- x %>%
      arrange(desc(get(colors[i]))) %>%
      mutate(!!paste0(tissue[i], "_row") := case_when(get(colors[i]) == 255 ~ Y/factor, TRUE ~ 0),
             !!paste0(tissue[i], "_col") := case_when(get(colors[i]) == 255 ~ X/factor, TRUE ~ 0))
  }

  dist <- x %>%
    filter_at(vars(ends_with("_row"), ends_with("_col")), any_vars(. > 0)) %>%
    as.data.frame()

  return(dist)
}

colors <- c("Red", "Green", "Blue") # matchLines("fiji_output.csv")
tissue <- c("MarrowFibrin", "Periosteum", "Cortical")

data <- SpatialCalc("R/fiji_output.csv", factor = 0.087658, colors = colors, tissue = tissue) # round up factor
View(data)

#' @file - cluster coordinates file
#' @id - cluster name

SpatialTime <- function(file = "", id = NULL) {

  tissue_to <- read.csv(file)
  tissue_to <- tissue_to %>%
    mutate(barcode = X) %>%
    select(imagerow, imagecol, barcode)

  tissue_from <- data %>%
    select(contains(id)) %>%
    filter_if(is.numeric, all_vars((.) != 0))

  st_calc <- tissue_to %>%
    rowwise() %>%
    mutate(st_abs = min(sqrt((imagerow - tissue_from[[paste0(id, "_row")]]) ^ 2 + (imagecol - tissue_from[[paste0(id, "_col")]]) ^ 2))) %>%
    ungroup() %>%
    mutate(st_rel = (st_abs - min(st_abs)) / (max(st_abs) - min(st_abs)))

  return(st_calc)
}

# Calculate SpatialTime ----
sp <- SpatialTime("R/excel/all_coordinates.csv", id = "MarrowFibrin")
View(sp)


#' @file - Seurat object
#' @st.calc - Spatialtime values
#' @spatial.by - Absolute or relative values
#' @slice - Select tissue slice
#' @return_obj - Return object

SpatialVis <- function(file = NULL, st.calc = NULL, spatial.by = c("abs", "rel"), slice = "slice1", return_obj = F) {

  if (is.null(file) || is.null(st.calc)) {
    stop("Both 'file' and 'st.calc' must be provided.")
  }

  spatial.by <- match.arg(spatial.by)

  myBarcode <- rownames(file@meta.data)
  TissueID <- st.calc[match(myBarcode, st.calc$barcode), ]

  if (spatial.by == "abs") {
    file$st <- TissueID$st_abs
  } else if (spatial.by == "rel") {
    file$st <- TissueID$st_rel
  }

  file$st[is.na(file$st)] <- 0

  if (return_obj == T) {
    return(file)
  }

  SpatialFeaturePlot(file, features = "st", images = slice)

}

vis <- SpatialVis(fracture, sp, spatial.by = c("rel"), return_obj = T)
View(vis@meta.data)


# Add module scores ----
osteo <- c("Alpl","Bglap","Bglap2","Col1a1","Col1a2","Dmp1","Ibsp","Mef2c","Postn",
           "Runx2","Sp7","Sparc","Phex","Satb2","Pth1r","Ostn","Car3")

angio <- c("Cdh5","Ddit3","S1pr1","Ankrd17","Lyl1","Mmp2","Reck","Acvrl1")

gene_modules <- list(Module1 = osteo, Module2 = angio)
module <- AddModuleScore(vis, features = gene_modules, assay = "SCT", name = "cluster")
View(module@meta.data)


#saveRDS(vis, "fracture_original_annotation_to_monocle.rds")

#' @file - Seurat object with module annotation
#' @assay - Seurat Assay
#' @min_expr - minimum gene expresion
#' @min_cells - Minimum cells expression
#' @mean_expr - Mean gene expression

PseudoTime <- function(file = NULL, assay = "RNA", min_expr = 0.1, min_cells = 2, mean_expr = 0.1) {
  if (is.null(file)) {
    stop("Input file is NULL. Please provide a valid Seurat object.")
  }

  if (!assay %in% names(file@assays)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }

  data <- as(as.matrix(file[[assay]]@data), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = file@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)

  HSMM <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)

  HSMM <- detectGenes(HSMM, min_expr = min_expr)
  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= min_cells))

  disp_table <- dispersionTable(HSMM)
  ordering_genes <- subset(disp_table, mean_expression >= mean_expr & 2*dispersion_empirical >= dispersion_fit)$gene_id
  HSMM <- setOrderingFilter(HSMM, ordering_genes)

  HSMM <- reduceDimension(HSMM, max_components=2)
  HSMM <- orderCells(HSMM, reverse=FALSE)

  file$ps <- HSMM@phenoData@data[["Pseudotime"]]

  return(file)
}

ps <- PseudoTime(file = module, assay = "SCT", min_expr = 0.1, min_cells = 5, mean_expr = 0.1) # vis
View(ps@meta.data)


# Get specific gene expression values
gene1 <- FetchData(ps, vars = "PDGFRA")
ps$PDGFRA <- gene1

gene2 <- FetchData(ps, vars = "HBA1")
ps$HBA1 <- gene2
View(ps@meta.data)

gene3 <- FetchData(ps, vars = "SOX9")
ps$SOX9 <- gene3
View(ps@meta.data)


# Visualization ----
ggplot(ps@meta.data) +
  geom_smooth(aes(st, HBA1, color = "HBA1"), method = "loess", span = 1.5, se = FALSE) +
  geom_smooth(aes(st, PDGFRA, color = "PDGFRA"), method = "loess", span = 1.5, se = FALSE) +
  geom_smooth(aes(st, SOX9, color = "SOX9"), method = "loess", span = 1.5, se = FALSE) +
  scale_color_manual(values = c("HBA1" = "#999999", "PDGFRA" = "#E69F00", "sox9" = "#56B4E9")) +
  theme_classic()

ggplot(ps@meta.data) +
  geom_smooth(aes(st, cluster1, color = "cluster1"), method = "loess", span = 1.5, se = FALSE) +
  theme_classic()



################################################################################

GeneVis <- function(data = NULL, genes = NULL, colors = NULL, st = NULL, span = NULL, se = FALSE) {

  for (gene in genes) {
    gene_data <- FetchData(data, vars = gene)
    data[[gene]] <- gene_data

    #data <- data@meta.data %>%
    #  select(st, all_of(genes))

    #m <- melt(data, id.vars = "st")

    #ggplot(m, aes(st,value, col=variable)) +
    #geom_smooth(span = 1.5, se = FALSE) +
    #theme_classic()
  }

  return(data)
}

ds <- GeneVis(ps, genes = c("HBA1","PDGFRA","SOX9"))
View(ds@meta.data)

genes <- c("HBA1","PDGFRA","SOX9")
class(genes)

ds <- ds@meta.data %>%
  select(st, all_of(genes))

m <- melt(metadata, id.vars = "st")
View(m)

ggplot(m, aes(st,value, col=variable)) +
  geom_smooth(span = 1.5, se = FALSE) +
  theme_classic()

s <-
  "A       B        C       G       Xax
0.451   0.333   0.034   0.173   0.22
0.491   0.270   0.033   0.207   0.34
0.389   0.249   0.084   0.271   0.54
0.425   0.819   0.077   0.281   0.34
0.457   0.429   0.053   0.386   0.53
0.436   0.524   0.049   0.249   0.12
0.423   0.270   0.093   0.279   0.61
0.463   0.315   0.019   0.204   0.23
"
d <- read.delim(textConnection(s), sep="")

library(ggplot2)
library(reshape2)
d <- melt(d, id.vars="Xax")

# Everything on the same plot
ggplot(d, aes(Xax,value, col=variable)) +
  geom_point() +
  stat_smooth()

View(metadata)
SpatialFeaturePlot(ps, features = "cluster2")
