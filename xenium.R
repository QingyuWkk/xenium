library(Seurat)
library(ggplot2)

path <- "C:/Users/AlainaKK/Downloads/Xenium/output-XETG00074__0010499__IC071922__20231001__221931"
# Load the Xenium data
xenium_obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium_obj <- subset(xenium_obj, subset = nCount_Xenium > 0)

# Extract Cell Centroids and Molecule Coordinates
# The counts matrix (“matrix”): This contains expression data for cells and features.
# Cell centroids in pixel coordinate space (“centroids”): Provides cell centroid coordinates (x, y).
# Molecule pixel coordinates (“microns”): Gives the pixel coordinates of individual molecules.
xenium_read <- ReadXenium(path, outs = c("matrix", "microns"), type = "centroids", mols.qv.threshold = 20)

# Run UMAP to reduce dimensions for visualization
xenium_obj <- SCTransform(xenium_obj, assay = "Xenium")
xenium_obj <- RunPCA(xenium_obj, npcs = 30, features = rownames(xenium_obj))
xenium_obj <- RunUMAP(xenium_obj, dims = 1:30)

# Identify clusters of cells
xenium_obj <- FindNeighbors(xenium_obj, reduction = "pca", dims = 1:30)
xenium_obj <- FindClusters(xenium_obj, resolution = 0.3)

# Function to compute each cell's area and correlate it with clusters
cell_area_and_cluster <- function(xenium_obj, cluster = "seurat_clusters",plot_options= "Violin", add_p_val = F,exc_cell_type = F) {
  # This function gets processed Xenium objects and give several plot outputs
  
  # get area for each cell
  xenium_obj[['cell_area']] <- sapply(xenium_obj@images[["fov"]]@boundaries[["segmentation"]]@polygons,function(x) x@area)
  
  # one can change a cell_type or another metadata column data with seurat_clusters
  ## Does it mean to add annotation for the features or to add new clustering method??
  xenium_obj@meta.data[["seurat_clusters"]] <- cluster
  
  # We have violin, boxplot, violin + boxplot options
  if (plot_options == "Violin") {
    ggplot2::ggplot(xenium_obj@meta.data, mapping=aes(x=seurat_clusters,y=cell_area,fill=seurat_clusters)) + geom_violin(trim=FALSE)
  } else if (plot_options == "Boxplot") {
    ggplot2::ggplot(xenium_obj@meta.data, mapping=aes(x=seurat_clusters,y=cell_area,fill=seurat_clusters)) + geom_boxplot(width=0.2,position = position_dodge(0.9))
  } else if (plot_options == "ViolinBoxplot") {
    ggplot2::ggplot(xenium_obj@meta.data, mapping=aes(x=seurat_clusters,y=cell_area,fill=seurat_clusters)) + geom_violin(trim=FALSE,color="white") + geom_boxplot(width=0.1,position = position_dodge(0.9))
  }
  
  # If add_p_val is True add pairwise comparison
  if (add_p_val) {
    pairwise.t.test(x = xenium_obj@meta.data$cell_area, g = xenium_obj@meta.data$seurat_clusters, p.adjust.method = "BH")
  }
  
  # One can exclude several clusters / Cell Types by using exc_cell_type
  if (exc_cell_type != FALSE) {
    xenium_obj@meta.data <- subset(xenium_obj@meta.data, subset = !(seurat_clusters %in% exc_cell_type))
  }
  
  # Return the modified Seurat object, which now includes cell area metadata
  return(xenium_obj)
}


