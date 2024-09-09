# In this script I run spatialtime by subsetting clusters. Don't use line references ----

# libraries ----
library(Seurat)
library(tidyverse)
library(Hmisc)
library(jsonlite)
library(monocle)
library(fs)


usethis::use_test("subsetLabels")
devtools::test()

# Import dataset ----
#vis <- readRDS("D10_young_DRG_mapped_stvalues.rds")
#View(fracture@meta.data)
#SpatialFeaturePlot(fracture, features = "PEP_rel", images = "slice1.3")

fracture <- readRDS("R/fracture_original_annotation.rds") # D10_young_DRG_mapped_st.rds
View(fracture@meta.data)

is(fracture, "Seurat")
#Idents(fracture) <- fracture$LTMR_status
SpatialDimPlot(fracture, images = "slice1")

#' @file - Seurat Object
#' @cluster - cluster names in Idents
#' @export.all - Export coordinates files
#' @slice.n - slice number
#' @dir.out - directory output

subsetLabels <- function(data = NULL, cluster = NULL, export.all = TRUE,
                         dir.out = "") {

  if (is.null(data)) {
    stop("Data is not provided.")
  }

  if (export.all) {

    # Check if it's Seurat data
    if (is.null(cluster)) {
      stop("Cluster IDs were not specified.")
    }

    if (length(cluster) > 3) {
      stop("Number of clusters should be less than or equal to 3 at a time.")
    }

    if (all(cluster %in% Idents(data))) {
      clstrs <- list()

      for (i in seq_along(cluster)) {
        clstrs[[i]] <- subset(data, idents = cluster[i])
      }

      for (z in seq_along(clstrs)) {
        for (y in seq_along(names(clstrs[[z]]@images))) {

          file_name <- paste0(dir.out, cluster[z], "_", names(clstrs[[z]]@images)[y], "_coordinates.csv")
          write.csv(clstrs[[z]]@images[[y]]@coordinates, file_name)
        }
      }

      return(clstrs)

    } else {

      stop("One or more specified clusters are not present in the data.")
    }

  } else {

    for (y in seq_along(names(data@images))) {

      file_name <- paste0(dir.out, names(data@images)[y], "_coordinates.csv")
      write.csv(data@images[[y]]@coordinates, file_name)
    }

    return(data)
  }
}

# if "export.all = T", then "cluster" parameter is not necessary
subsetLabels(fracture, cluster = c("Periosteum"), export.all = F, dir.out = "R/All_peri/") # Add error handling for terms that do not exist or poorly written
# LTMRpos R/LTMR/
# Function to merge all exported coordinates together ----

#' @file - Direcotry with slices .csv file
#' @pattern - Regex pattern to select specific files within directory

CoordMerge <- function(files = "", pattern = NULL) {

  csv_files <- str_sort(fs::dir_ls(files), pattern = NULL)
  rdr <- readr::read_csv(csv_files)

  return(rdr)
}


peri <- CoordMerge(file = "R/Periosteum/")
View(peri)

drg <- CoordMerge(file = "R/DRG/")
View(drg)

ltmr <- CoordMerge(file = "R/LTMR/")
View(ltmr)

motor <- CoordMerge(file = "R/MOTOR/")
View(motor)

pep <- CoordMerge(file = "R/PEP/")
View(pep)

np <- CoordMerge(file = "R/NP/")

all <- CoordMerge(file = "R/All/")
View(all)

all_peri <- CoordMerge(file = "R/All_peri/")
View(all)

# Reorganize the function to get coordinates without lines drawn ----
#' @file - cluster coordinates file
#' @id - cluster name

SpatialTime <- function(file = NULL, file2 = NULL, id = NULL) {

  tissue_to <- file %>%
    mutate(barcode = ...1) %>%
    select(imagerow, imagecol, barcode) #%>%
    #filter(grepl("1_1$", barcode))

  tissue_from <- file2 %>%
    mutate(barcode = ...1) %>%
    select(imagerow, imagecol, barcode) #%>%
    #filter(grepl("1_1$", barcode))

  st_calc <- tissue_to %>%
    rowwise() %>%
    mutate(st_abs = min(sqrt((imagerow - tissue_from[["imagerow"]]) ^ 2 + (imagecol - tissue_from[["imagecol"]]) ^ 2))) %>%
    ungroup() %>%
    mutate(st_rel = (st_abs - min(st_abs)) / (max(st_abs) - min(st_abs)))

  return(st_calc)
}

# Usage
sp <- SpatialTime(all_peri, peri)
View(sp)


#' @file - Seurat object
#' @st.calc - Spatialtime values
#' @spatial.by - Absolute or relative values
#' @slice - Select tissue slice
#' @return_obj - Return object

SpatialVis <- function(file = NULL, st.calc = NULL, spatial.by = c("abs", "rel"), id = "", slice = "slice1", return_obj = FALSE) {

  if (is.null(file) || is.null(st.calc)) {
    stop("Both 'file' and 'st.calc' must be provided.")
  }

  spatial.by <- match.arg(spatial.by, several.ok = TRUE)

  myBarcode <- rownames(file@meta.data)
  TissueID <- st.calc[match(myBarcode, st.calc$barcode), ]

  for (type in spatial.by) {
    column_name <- paste0("st_", type)
    file[[column_name]] <- TissueID[[column_name]]
    file[[column_name]][is.na(file[[column_name]])] <- 0

    file@meta.data <- file@meta.data %>%
      rename_with(~ paste0(id, "_", type), all_of(column_name))
  }

  if (return_obj) {
    return(file)
  }

  SpatialFeaturePlot(file, features = paste0(id, "_", spatial.by), images = slice)
}

# Example usage:
vis <- SpatialVis(fracture, sp, spatial.by = c("abs", "rel"), id = "Peri", return_obj = TRUE)
View(vis@meta.data)


saveRDS(vis, "R/D10_young_DRG_mapped_st.rds")


vis@meta.data <- vis@meta.data %>% rename_at("st", ~"LTMR_abs_slice1")
SpatialFeaturePlot(vis, features = "LTMR_abs_slice1", images = "slice1")
saveRDS(vis, "D10_young_DRG_mapped_st.rds")

SpatialFeaturePlot(vis, features = "Peri_abs", images = "slice1")



# Add module scores ----
osteo <- c("Alpl","Bglap","Bglap2","Col1a1","Col1a2","Dmp1","Ibsp","Mef2c","Postn",
           "Runx2","Sp7","Sparc","Phex","Satb2","Pth1r","Ostn","Car3")

angio <- c("Cdh5","Ddit3","S1pr1","Ankrd17","Lyl1","Mmp2","Reck","Acvrl1")

gene_modules <- list(Module1 = osteo, Module2 = angio)
module <- AddModuleScore(vis, features = gene_modules, assay = "SCT", name = "cluster")
View(module@meta.data)

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

saveRDS(ps, "R/D10_young_DRG_mapped_st.rds")
saveRDS(ps, "D10_young_DRG_mapped_stvalues_monocle.rds")
saveRDS(ps, "fracture_original_annotation_for_monocle.rds")


gene1 <- FetchData(ps, vars = "PDGFRA")
ps$SOX9 <- gene1

# Visualization ----
metadata <- ps@meta.data

ggplot(metadata, aes(Peri_rel, SOX9)) +
  geom_smooth(method = "loess", span = 2, se = F) +
  theme_classic()

SpatialFeaturePlot(ps, features = "cluster2")
