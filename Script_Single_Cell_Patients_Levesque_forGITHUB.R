# Thelma Poluha reanalyses Mitch Levesque collaboration
library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(GSEABase)
library(GSA)
library(dplyr)
library(MAST)
library(StanHeaders)
library(rstan)
library(Rcpp)
library(RcppEigen)
library(future)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(tidyr)
library(tibble)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(future)
library(readr)
library(AUCell)





# Directory containing the processed files
directory <- "~/Single_Cell_Patients_Levesque"


# List of file names
files <- c("HG2019_2182_raw.rds", "HG2019_2532_raw.rds", "HG2019_2753_raw.rds", 
           "HG2019_4830_raw.rds", "HG2019_867_raw.rds", "HG2022_1844_raw.rds", 
           "HG2022_3845_raw.rds", "N1_raw.rds", "N3_raw.rds") 

# Directory where the files are located
directory <- "~/Single_Cell_Patients_Levesque"  


# For loop to process each file
for (file in files) {
  
  # 1. Load the file
  file_path <- file.path(directory, file)  # Full path to the file
  obj <- readRDS(file_path)  # Load the .rds file
  
  # Check if 'scDblFinder.class' is present in the file
  if ("scDblFinder.class" %in% colnames(colData(obj))) {
    print(paste("Processing scDblFinder.class for file:", file))
    
    # 1a. Assign scDblFinder.class if it exists in the metadata
    obj$scDblFinder.class <- colData(obj)$scDblFinder.class
    
    # 1b. Display a table of scDblFinder classifications
    print(table(obj$scDblFinder.class))
  } else {
    warning(paste("Warning: scDblFinder.class not found in file:", file))
  }
  
  # 2. Convert Ensembl IDs to gene symbols
  gene_symbols <- rowData(obj)$Symbol  # Extract gene symbols
  if (is.null(gene_symbols)) {
    stop(paste("Error: The file", file, "does not have gene symbols in rowData."))
  }
  rownames(obj) <- gene_symbols  # Set gene symbols as rownames
  
  # 3. Handle duplicate gene names
  if (sum(duplicated(rownames(obj))) > 0) {  # If there are duplicates
    rownames(obj) <- make.unique(rownames(obj))  # Make names unique if duplicates are found
  }
  
  # 4. Convert to Seurat object
  obj <- as.Seurat(obj, counts = "counts")  # Convert to Seurat object
  
  # Make sure scDblFinder.class is in the meta data after conversion to Seurat
  if (!"scDblFinder.class" %in% colnames(obj@meta.data)) {
    obj@meta.data$scDblFinder.class <- colData(obj)$scDblFinder.class
  }
  
  # 5. Change project name based on the file name (without .rds extension)
  obj$orig.ident <- gsub("\\.rds$", "", file)  # Set project name to file name without .rds
  
  # 6. Filter the metadata
  cols_to_select <- c("orig.ident", "nCount_originalexp", "nFeature_originalexp", "Sample", "Barcode", "scDblFinder.class")
  
  # Select the available columns (handle missing columns with warnings)
  missing_cols <- setdiff(cols_to_select, colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    warning(paste("Warning: The following columns are missing in file", file, ":", 
                  paste(missing_cols, collapse = ", ")))
  } else {
    obj@meta.data <- obj@meta.data %>%
      dplyr::select(all_of(cols_to_select))
  }
  
  # 7. Add mitochondrial gene percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # 8. Visualize quality control (QC) plots
  VlnPlot(obj, features = c("nCount_originalexp", "nFeature_originalexp", "percent.mt"), 
          ncol = 3, pt.size = 0)
  
  # 9. Subset the object to keep only singlets
  obj <- subset(obj, subset = scDblFinder.class == "singlet")
  
  # 10. Save or assign the processed object
  assign(gsub("\\.rds$", "", file), obj)  # Assign the object to a variable based on the file name
  
  # Optional: Save the processed object for future use
  saveRDS(obj, file = file.path(directory, paste0("processed_", file)))  # Save the processed object
}





########### Load the processed Seurat objects ################
# Since you've already processed the objects and saved them, let's load them directly

# List of processed files
processed_files <- c("processed_HG2019_2182_raw.rds", "processed_HG2019_2532_raw.rds", 
                     "processed_HG2019_2753_raw.rds", "processed_HG2019_4830_raw.rds", 
                     "processed_HG2019_867_raw.rds", "processed_HG2022_1844_raw.rds", 
                     "processed_HG2022_3845_raw.rds", "processed_N1_raw.rds","processed_N3_raw.rds")


# Load the Seurat objects
seurat_objects <- lapply(processed_files, function(file) readRDS(file.path(directory, file)))

# Assign sample names to each object based on the file names (without the "processed_" prefix)
names(seurat_objects) <- gsub("processed_", "", processed_files)

names(seurat_objects) <- gsub("_raw.rds", "", names(seurat_objects))

########### Check and update mitochondrial gene percentage if needed ################
# This step might have already been done during processing. 
# If needed, you can calculate or verify the mitochondrial gene percentage:

for (i in 1:length(seurat_objects)) {
  seurat_objects[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
}

########### Merging the Seurat Objects ################
# Here, we will merge all the loaded Seurat objects into a single one for downstream analysis
merged_obj <- Reduce(function(x, y) merge(x, y, project = "Merged_Seurat_Project"), seurat_objects)

# Check the summary of the merged object
summary(merged_obj)
print(merged_obj)
head(merged_obj@meta.data)


table(merged_obj$orig.ident)
#saveRDS(merged_obj, file = file.path(directory, paste0("merged_obj_unfiltered.rds")))
merged_obj = readRDS(file = file.path(directory, paste0("merged_obj_unfiltered.rds")))
#view(merged_obj@meta.data) cannot be executed because the object is very big for R studio to handle


########### Quality Control (QC) Plots ################
# This might have already been done during the processing, but let's visualize it for the merged object
VlnPlot(merged_obj, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3, pt.size = 0.0)
VlnPlot(merged_obj, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0.0)

# see paper for mt.rna QC
# we decide to set it at 25% to be less stringent (nevi samples seem to have a higher fraction)
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8599307/
VlnPlot(merged_obj, features = c("percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0.0)
VlnPlot(merged_obj, features = c("percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0.0, y.max = 25)

# for nFeatures we take out after 8000


##################PCA BEFORE NORMALIZATION##################

# Using Seurat's built-in human cell cycle genes
s.genes <- cc.genes.updated.2019$s.genes  # S phase genes
g2m.genes <- cc.genes.updated.2019$g2m.genes  # G2M phase genes

# Applying cell cycle scoring on the merged Seurat object
merged_obj <- CellCycleScoring(merged_obj, 
                               s.features = s.genes, 
                               g2m.features = g2m.genes, 
                               set.ident = TRUE)

# Normalize the data with SCTransform, using the originalexp assay and regressing out mitochondrial percentage, S phase, and G2M phase scores
merged_obj <- SCTransform(merged_obj, 
                          assay = "originalexp",  # Specify the correct assay
                          verbose = FALSE, 
                          vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

# Perform PCA
merged_obj <- RunPCA(merged_obj, verbose = TRUE)



# Scatterplot of the first two principal components
DimPlot(merged_obj, reduction = "pca", dims = c(1, 2), group.by = "orig.ident")
















########### Subsetting Based on Quality Control ################
# This step subsets cells based on RNA features and mitochondrial gene percentage.
# We filter out cells with too few or too many detected features and too much mitochondrial content.

merged_obj_filtered <- subset(merged_obj, subset = nFeature_originalexp > 500 & nFeature_originalexp < 8000 & percent.mt < 25)
merged_obj_filtered <- SCTransform(merged_obj_filtered, assay = "originalexp", verbose = TRUE, vars.to.regress = c("percent.mt"))

########### SCTransform Normalization ################
# If SCTransform normalization hasn't already been done on your `processed_` files, let's do it now.

#merged_obj <- SCTransform(merged_obj, verbose = TRUE, vars.to.regress = c("percent.mt"))




####################################################
################### SCT with CELL CYCLE ############
####################################################

# Using Seurat's built-in human cell cycle genes
s.genes <- cc.genes.updated.2019$s.genes  # S phase genes
g2m.genes <- cc.genes.updated.2019$g2m.genes  # G2M phase genes

# Applying cell cycle scoring on the merged Seurat object
merged_obj_filtered <- CellCycleScoring(merged_obj_filtered, 
                                        s.features = s.genes, 
                                        g2m.features = g2m.genes, 
                                        set.ident = TRUE)

# Normalize the data with SCTransform, using the originalexp assay and regressing out mitochondrial percentage, S phase, and G2M phase scores
merged_obj_filtered <- SCTransform(merged_obj_filtered, 
                                   assay = "originalexp",  # Specify the correct assay
                                   verbose = FALSE, 
                                   vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

# Perform PCA
merged_obj_filtered <- RunPCA(merged_obj_filtered, verbose = TRUE)



# Scatterplot of the first two principal components
DimPlot(merged_obj_filtered, reduction = "pca", dims = c(1, 2), group.by = "orig.ident")


# Feature plot of PC_1 and PC_2 scores
FeaturePlot(merged_obj_filtered, features = c("PC_1", "PC_2"), reduction = "pca")





p1 <- DimPlot(merged_obj, reduction = "pca", dims = c(1, 2), group.by = "orig.ident")
p2 <- DimPlot(merged_obj_filtered, reduction = "pca", dims = c(1, 2), group.by = "orig.ident")
CombinePlots(list(p1, p2))






# Run clustering on PCs
merged_obj_filtered <- FindNeighbors(merged_obj_filtered, dims = 1:15)
merged_obj_filtered <- FindClusters(merged_obj_filtered, resolution = 0.7)

# Visualize clusters in PCA space
DimPlot(merged_obj_filtered, reduction = "pca", group.by = "seurat_clusters")


# Save the filtered Seurat object (after cell cycle scoring and clustering)
#saveRDS(merged_obj_filtered, file = file.path(directory, paste0("merged_obj_filtered_after_CC_and_clustering.rds")))



p1 <- DimPlot(merged_obj, reduction = "pca", dims = c(1, 2), group.by = "orig.ident")
p2 <- DimPlot(merged_obj_filtered, reduction = "pca", dims = c(1, 2), group.by = "orig.ident")
CombinePlots(list(p1, p2))



# Harmony integration to remove batch effects based on 'orig.ident'
merged_obj_filtered <- RunHarmony(merged_obj_filtered, group.by.vars = "orig.ident", assay.use = "SCT")

# Checkpoint: Save the processed object
#saveRDS(merged_obj_filtered, file = file.path(directory, "merged_obj_filtered_after_CC_reg.rds"))

# PCA colored by cell cycle phase
DimPlot(merged_obj_filtered_15PCs, group.by = "Phase", reduction = "pca", pt.size = 0.02) + ggtitle("PCA - Cell Cycle Phases")


# List objects in the environment
ls()

# Try with 15 PCs
merged_obj_filtered_15PCs <- FindNeighbors(merged_obj_filtered, dims = 1:15, reduction = "harmony")
merged_obj_filtered_15PCs <- FindClusters(merged_obj_filtered_15PCs, resolution = 0.7, reduction = "harmony")
merged_obj_filtered_15PCs <- RunUMAP(merged_obj_filtered_15PCs, dims = 1:15, reduction = "harmony")


# Save the processed object
#saveRDS(merged_obj_filtered_15PCs, file = file.path(directory, "merged_obj_15PCs.rds"))

# Load the saved object if needed later
merged_obj_filtered_15PCs <- readRDS(file.path(directory, "merged_obj_15PCs.rds"))





############### Distribution of Pc1 by batch ############### 

RidgePlot(merged_obj_filtered_15PCs, features = "PC_1", group.by = "orig.ident") +
  ggtitle("Distribution of PC1 by Batch")







# Visualize PCA with cells colored by cell cycle phase
pca_plot <- DimPlot(merged_obj_filtered_15PCs, 
                    reduction = "pca", 
                    group.by = "Phase", 
                    pt.size = 0.5) + 
  ggtitle("PCA - Cell Cycle Phases 15 PCA") +
  theme_minimal() +  # for a clean, minimal background
  labs(color = "Cell Cycle Phase") +  # add a label to the color legend
  theme(legend.position = "right")  # position the legend on the right

# Display the plot
print(pca_plot)






####################################################
################### HEATMAP + ELBOW PLOTS ##########
####################################################

# Generate a Harmony heatmap
#pdf(file.path(directory, "Harmony_Heatmap_all_after_CC_reg.pdf"), width = 7, height = 7)
#harmony_embeddings <- Embeddings(merged_obj_filtered, 'harmony')
#col_fun <- colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
#Heatmap(harmony_embeddings, 
#cluster_rows = TRUE, 
#cluster_columns = FALSE,  
#clustering_distance_columns = "euclidean",
#clustering_method_columns = "complete", 
#show_column_names = TRUE,
#show_row_names = FALSE,
#name = "Harmony_embedding")
#dev.off()

# Elbow plot to decide the number of PCs to use
ElbowPlot(merged_obj_filtered, ndims = 50)

# Calculate percentage of variance explained by each principal component
eigValues <- (merged_obj_filtered[["pca"]]@stdev)^2  # Eigenvalues from the PCA
varExplained <- eigValues / sum(eigValues)  # Percentage of variance explained by each PC
plot(varExplained, type = "b", xlab = "Principal Component", ylab = "Variance Explained")





# Set ComplexHeatmap options to suppress messages
ht_opt$message = FALSE

# Define color function for the heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))

# Select the top 50 most variable genes from harmony_embeddings
top_genes <- head(order(apply(harmony_embeddings, 1, var), decreasing = TRUE), 100)
harmony_embeddings_top50 <- harmony_embeddings[top_genes, ]

# Generate a clearer heatmap
Heatmap(harmony_embeddings_top50, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,  # Cluster columns to group similar cells together
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        show_column_names = FALSE,  # Hide column names to reduce clutter
        show_row_names = TRUE,      # Show row names for the selected genes
        row_names_gp = gpar(fontsize = 8),  # Adjust row label font size for clarity
        name = "Harmony_embedding", 
        col = col_fun)







####################################################
################### UMAP AND CLUSTERING ###########
####################################################

# Try different numbers of PCs for clustering and UMAP projection
# First, with 15 PCs
merged_obj_filtered_15PCs <- FindNeighbors(merged_obj_filtered, dims = 1:15, reduction = "harmony")
merged_obj_filtered_15PCs <- FindClusters(merged_obj_filtered_15PCs, resolution = 0.7, reduction = "harmony") #higher resolution numbrer = more clusters
merged_obj_filtered_15PCs <- RunUMAP(merged_obj_filtered_15PCs, dims = 1:15, reduction = "harmony")


# UMAP plots for the merged data
DimPlot(merged_obj_filtered_15PCs, group.by = "seurat_clusters", pt.size = 0.02, label = TRUE)
DimPlot(merged_obj_filtered_15PCs, group.by = "orig.ident", pt.size = 0.02, label = FALSE, split.by = "orig.ident")
DimPlot(merged_obj_filtered_15PCs, group.by = "orig.ident", pt.size = 0.02, label = FALSE)
FeaturePlot(merged_obj_filtered_15PCs, features = c("LPAR1"), pt.size = 0.02, label = T, order = T)
DotPlot(merged_obj_filtered_15PCs, features = c("PMEL", "MLANA","MITF"))

# Save the processed object for future use
saveRDS(merged_obj_filtered_15PCs, file = file.path(directory, "(merged_obj_filtered_15PCs.rds"))









# Get unique sample identifiers
samples <- unique(merged_obj_filtered_15PCs$orig.ident)

# Loop through each sample and create UMAP plots
for (sample in samples) {
  # Subset the data for the specific sample
  sample_data <- subset(merged_obj_filtered_15PCs, orig.ident == sample)
  
  # Generate the UMAP plot
  p <- DimPlot(sample_data, pt.size = 0.02) + ggtitle(sample)
  
  # Save the plot as a PNG file
  ggsave(paste0("UMAP_", sample, ".png"), plot = p)
}









######### MARKERS ############

#plan()
future::plan("multisession", workers = 20)
options(future.globals.maxSize= 891289600)
markers_melanoma_vs_nevi <- FindAllMarkers(merged_obj_filtered_15PCs, only.pos = T, verbose = T)
plan("sequential")

merged_obj_filtered_15PCs@misc = markers_melanoma_vs_nevi
write_csv(markers_melanoma_vs_nevi, file = file.path(directory, "all_markers_clusters.csv"))

top20_genes = markers_melanoma_vs_nevi %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
top20_genes = top20_genes[order(top20_genes$avg_log2FC, decreasing = T),]
top20_genes = top20_genes[order(top20_genes$cluster),]


write_csv2(top20_genes, file = file.path(directory, "top20_markers_clusters.csv"))

# Open the markers table in a new tab in RStudio
View(markers_melanoma_vs_nevi)

# Alternatively, for the top 20 markers per cluster
View(top20_genes)



# Define the directory where you want to save the files
output_directory <- "~/Single_Cell_Patients_Levesque"

# Save the markers table as a CSV file
write.csv(markers_melanoma_vs_nevi, file = file.path(output_directory, "markers_melanoma_vs_nevi.csv"), row.names = FALSE)

# Save the top 20 markers per cluster as a CSV file
write.csv(top20_genes, file = file.path(output_directory, "top20_genes.csv"), row.names = FALSE)



# subset smaller set of cells
downsampled.obj <- merged_obj_filtered_15PCs[, sample(colnames(merged_obj_filtered_15PCs), size = 20000, replace=F)]
DoHeatmap(downsampled.obj, features = top20_genes$gene)




#20 genes is too much for the computer power we have and the heatmap does not show.
#Heatmap
DoHeatmap(merged_obj_filtered_15PCs, features = top20_genes$gene)






# Extract top 50 genes per cluster
top50_genes_per_cluster <- markers_melanoma_vs_nevi %>%
  group_by(cluster) %>%
  top_n(50, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Directory to save the heatmaps
output_dir <- "~/Single_Cell_Heatmaps"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Generate and save one heatmap per cluster
clusters <- unique(top50_genes_per_cluster$cluster)

for (cluster_id in clusters) {
  # Filter genes for the current cluster
  cluster_genes <- top50_genes_per_cluster %>%
    filter(cluster == cluster_id) %>%
    pull(gene)
  
  # Generate the heatmap
  heatmap_plot <- DoHeatmap(
    merged_obj_filtered_15PCs,
    features = cluster_genes,
    group.by = "cluster"
  ) +
    ggtitle(paste("Heatmap for Cluster", cluster_id))
  
  # Save the heatmap to a file
  ggsave(
    filename = file.path(output_dir, paste0("Heatmap_Cluster_", cluster_id, ".png")),
    plot = heatmap_plot,
    width = 10,
    height = 8
  )
  
  # Optionally display in RStudio (comment out if not needed)
  print(heatmap_plot)
}




#We try with the top 10 genes but it does not work either. 
# Select the top 10 genes based on avg_log2FC from the top20_genes dataframe
top10_genes <- top20_genes %>%
  arrange(desc(avg_log2FC)) %>%  # Sort by average log2FC in descending order
  head(10)  # Take the top 10 genes

# Generate the heatmap using the top 10 genes
DoHeatmap(merged_obj_filtered_15PCs, features = top10_genes$gene)







###################################################################
################## AUCell for cell types ##########################
###################################################################
signatures_list <- list(
  #https://www.nature.com/articles/s42003-020-0922-4
  TCELLS = list(c("CD3E", "CD3G", "CD3D", "LCK")),
  TREGS= list(c("FOXP3", "CTLA4")),
  BCELLS= list(c("CD19", "CD22","CD72")),
  MACROPHAGES_DCs = list(c("LYZ", "AIF1", "HLA-DRA", "CD68", "ITGAX")),
  VASCULAR_EC = list(c("SELE", "CLDN5", "VWF", "CDH5")),
  LYMPHATIC_EC = list(c("CLDN5", "LYVE1", "PROX1")),
  PERICYTES = list(c("ACTA2", "RGS5", "PDGFRB")),
  MELANOCYTES = list(c("PMEL", "MLANA", "TYRP1", "DCT")),
  UNDIFFERENTIATED_KERATINOCYTES = list(c("KRT5", "KRT14", "TP63", "ITGB1", "ITGA6")),
  DIFFERENTIATED_KERATINOCYTES = list(c("KRT1", "KRT10", "SBSN", "KRTDAP")),
  FIBROBLASTS = list(c("LUM", "DCN", "VIM", "PDGFRA", "COL1A2"))
)


saveRDS(signatures_list, file = file.path(directory, "(signatures_list.rds"))


#Ranking for clusters
# AUCell plots and scores
#cells_rankings <- AUCell_buildRankings(merged_obj_filtered_15PCs@assays[["originalexp"]]@counts)

#Seguimos estos pasos para obtener el alineamiento, y nos va sacando imágenes de cada tipo celular. Obtenemos dos
#Una es los histogramas, donde nos muestra los AUC. Debería ser bimodal, siendo la parte derecha y más oscura nuestra
#población celular
# plots and scores
setwd("~/Single_Cell_Patients_Levesque")
dir.create("Feature_Plots")
dir.create("VLN_plots")
dir.create("Dot_plots")


#### FEATURE PLOT
# Cambia el directorio de trabajo
setwd("~/Single_Cell_Patients_Levesque/Feature_Plots/")


# Make sure to load the paletteer package at the start of your script
library(paletteer)

#### FEATURE PLOT
# Cambia el directorio de trabajo
setwd("~/Single_Cell_Patients_Levesque/Feature_Plots/")

for (i in 1:length(signatures_list)) {
  vector_signature <- as.vector(signatures_list[[i]][[1]])
  available_genes <- intersect(vector_signature, rownames(merged_obj_filtered_15PCs))
  if (length(available_genes) > 0) {
    geneSets <- GeneSet(available_genes, setName = names(signatures_list)[i])
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
    
    # Guardar histograma con umbral
    jpeg(filename = paste0("AUC_histogram_", names(signatures_list)[i], ".jpeg"))
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, nCores = 1, assign = TRUE)
    dev.off()
    
    # Obtener la puntuación en los metadatos
    signature_AUC <- getAUC(cells_AUC)
    signature_AUC <- t(signature_AUC)
    colnames(signature_AUC) <- names(signatures_list)[i]
    merged_obj_filtered_15PCs@meta.data <- cbind(merged_obj_filtered_15PCs@meta.data, signature_AUC)
    
    # Use paletteer_c to get the color scale
    p <- FeaturePlot(merged_obj_filtered_15PCs, features = names(signatures_list)[i], label = TRUE) +
      scale_color_gradientn(colors = paletteer_c("grDevices::Zissou 1", 30)) +  # Use the new color palette
      theme_minimal()
    
    # Save the plot
    ggsave(filename = paste0("AUC_plot_", names(signatures_list)[i], ".jpeg"), plot = p, width = 300, height = 200, units = "mm", dpi = 320)
  } else {
    message(paste("No genes available for signature:", names(signatures_list)[i]))
  }
}





#paletteer_c("grDevices::Zissou 1", 30)




#### VLN 
# Cambiar el directorio de trabajo
setwd("~/Single_Cell_Patients_Levesque/VLN_plots/")

for (i in 1:length(signatures_list)) {
  vector_signature <- as.vector(signatures_list[[i]][[1]])
  available_genes <- intersect(vector_signature, rownames(merged_obj_filtered_15PCs))
  if (length(available_genes) > 0) {
    # Crear y guardar el VlnPlot
    p <- VlnPlot(merged_obj_filtered_15PCs, features = available_genes, pt.size = 0)
    ggsave(filename = paste0("VLN_plot_", names(signatures_list)[i], ".jpeg"), plot = p, width = 900, height = 200, units = "mm", dpi = 320)
  } else {
    message(paste("No genes available for signature:", names(signatures_list)[i], "- Skipping this plot."))
  }
}



#### DOTPLOT---> Usa para firmas particulares, para el primer approach evitar
setwd("~/Single_Cell_Patients_Levesque/Dot_plots/")

for (i in 1:length(signatures_list)) {
  vector_signature <- as.vector(signatures_list[[i]][[1]])
  available_genes <- intersect(vector_signature, rownames(merged_obj_filtered_15PCs))
  if (length(available_genes) > 0) {
    # Crear y guardar el DotPlot
    p <- DotPlot(merged_obj_filtered_15PCs, features = available_genes)
    ggsave(filename = paste0("DotPlot_plot_", names(signatures_list)[i], ".jpeg"), plot = p, width = 200, height = 200, units = "mm", dpi = 320)
  } else {
    message(paste("No genes available for signature:", names(signatures_list)[i], "- Skipping this plot."))
  }
}





###### PRELIMINARY ANNOTATION ######

# Define the annotations for each cluster (0-20)
cluster_annotations <- c(
  "MELANOCYTIC",     #0
  "VASCULAR EC",     #1
  "T CELLS",         #2 
  "MELANOCYTIC",     #3
  "PERICYTES",       #4
  "MELANOCYTIC",     #5
  "MELANOCYTIC",     #6
  "MACROPHAGES/DCs", #7
  "T CELLS",         #8
  "FIBROBLASTS",     #9
  "MELANOCYTIC",     #10
  "KERATINOCYTES",   #11
  "MELANOMA",        #12
  "MELANOCYTIC",     #13
  "MELANOCYTIC",     #14
  "MELANOMA",        #15
  "MELANOCYTIC",     #16
  "CELL ADHESION",   #17
  "LYMPHATIC EC",    #18
  "NEVI SPECIFIC",   #19
  "MELANOMA"         #20
)

# Assign names to the vector based on cluster numbers (0 to 20)
names(cluster_annotations) <- 0:20

# Assign cluster annotations to the Seurat object metadata
# Ensure Idents are converted to character to match the annotation vector
merged_obj_filtered_15PCs$cluster_annotation <- cluster_annotations[as.character(Idents(merged_obj_filtered_15PCs))]

# Verify the annotation assignment
print(head(merged_obj_filtered_15PCs@meta.data))
print(table(merged_obj_filtered_15PCs$cluster_annotation))

# Visualize UMAP with annotated clusters (default colors)
DimPlot(merged_obj_filtered_15PCs, reduction = "umap", 
        group.by = "cluster_annotation", 
        order = T,
        label = TRUE,          # Adds labels to clusters
        label.size = 4,        # Adjust the label size for readability
        repel = TRUE,          # Repels labels to avoid overlapping
        pt.size = 0.4) +       # Adjust point size
  ggtitle("UMAP - Cluster Annotations") +
  theme_minimal() +
  theme(legend.position = "right")

# Save the annotated Seurat object for future use
saveRDS(merged_obj_filtered_15PCs, file = "annotated_merged_obj_filtered_15PCs.rds")


#Heatmap
DoHeatmap(merged_obj_filtered_15PCs, features = top20_genes$gene)







# Step 1: Map conditions to your samples and add the 'condition' column to metadata
condition_mapping <- c(
  "HG2019_2182_raw" = "melanoma",
  "HG2019_2532_raw" = "melanoma",
  "HG2019_2753_raw" = "melanoma",
  "HG2019_4830_raw" = "melanoma",
  "HG2019_867_raw" = "melanoma",
  "HG2022_1844_raw" = "nevi",
  "HG2022_3845_raw" = "nevi",
  "N1_raw" = "nevi",
  "N3_raw" = "melanoma"
)

# Add the condition column based on orig.ident
merged_obj_filtered_15PCs$condition <- condition_mapping[merged_obj_filtered_15PCs$orig.ident]

# Step 2: Verify that the condition column was added correctly
print(table(merged_obj_filtered_15PCs$condition))

# Step 3: Create a UMAP plot split by condition
DimPlot(merged_obj_filtered_15PCs, 
        reduction = "umap", 
        group.by = "seurat_clusters",  
        split.by = "condition",  
        label = FALSE,  
        pt.size = 0.4) +  
  ggtitle("UMAP - Clusters by Condition (Nevi vs Melanomas)") +
  theme_minimal() +
  theme(legend.position = "right")


# Define the file path for saving the object
output_file <- "~/Single_Cell_Patients_Levesque/merged_obj_filtered_15PCs_with_condition.rds"

# Save the updated Seurat object
saveRDS(merged_obj_filtered_15PCs, file = output_file)

# Confirm save
cat("Seurat object saved to:", output_file, "\n")





# Extract UMAP coordinates and cluster information
umap_coords <- merged_obj_filtered_15PCs[["umap"]]@cell.embeddings
cluster_info <- merged_obj_filtered_15PCs@meta.data$seurat_clusters

# Combine the coordinates with cluster information
umap_data <- data.frame(UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2], Cluster = cluster_info)

# Calculate the centroid of each cluster by averaging the UMAP coordinates
centroids <- umap_data %>%
  group_by(Cluster) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# Create the UMAP plot without labels initially
umap_plot <- DimPlot(merged_obj_filtered_15PCs, 
                     reduction = "umap", 
                     group.by = "seurat_clusters",  # Color by cluster number
                     split.by = "condition",  # Split the UMAP by condition
                     label = FALSE,  # Do not display automatic cluster labels
                     pt.size = 0.4) +  # Adjust point size for visibility
  ggtitle("UMAP - Clusters by Condition (Nevi vs Melanomas)") +
  theme_minimal() +
  theme(legend.position = "right")

# Add cluster numbers at the centroids
umap_plot <- umap_plot + 
  geom_text(data = centroids, 
            aes(x = UMAP_1, y = UMAP_2, label = as.factor(Cluster)),
            size = 5,  # Adjust the size of the cluster numbers
            color = "black", 
            fontface = "bold", 
            alpha = 0.8,  # Transparency of text
            check_overlap = TRUE)  # Avoid overlapping text labels

# Display the UMAP plot
print(umap_plot)



# Create separate plots for melanomas and nevi
# Melanoma
melanoma_obj <- subset(merged_obj_filtered_15PCs, subset = condition == "Melanoma")
umap_melanoma <- DimPlot(melanoma_obj, 
                         reduction = "umap", 
                         group.by = "cluster_annotation", 
                         order = T,
                         label = TRUE, 
                         repel = TRUE, 
                         pt.size = 0.4) +
  ggtitle("UMAP - Melanomas") +
  theme_minimal()

# Nevi
nevi_obj <- subset(merged_obj_filtered_15PCs, subset = condition == "Nevi")
umap_nevi <- DimPlot(nevi_obj, 
                     reduction = "umap", 
                     group.by = "cluster_annotation", 
                     order = T,
                     label = TRUE, 
                     repel = TRUE, 
                     pt.size = 0.4) +
  ggtitle("UMAP - Nevi") +
  theme_minimal()

# Display the separate plots
print(umap_melanoma)
print(umap_nevi)





FeaturePlot(merged_obj_filtered_15PCs, split.by = "condition", order = T, features = "MDK", label = T)
FeaturePlot(merged_obj_filtered_15PCs, split.by = "condition", order = T, features = "PTN", label = T)
FeaturePlot(merged_obj_filtered_15PCs, split.by = "condition", order = T, features = "DDX46", label = T)
FeaturePlot(merged_obj_filtered_15PCs, split.by = "condition", order = T, features = "SMC4", label = T)
FeaturePlot(merged_obj_filtered_15PCs, split.by = "condition", order = T, features = "STAU1", label = T)






################ UMAP BY CELL CYCLE PHASE ################


# Generate UMAP colored by cell cycle phases
DimPlot(merged_obj_filtered_15PCs, 
        reduction = "umap", 
        group.by = "Phase",  # "Phase" contains the cell cycle phase information
        label = TRUE,        # Add labels to the clusters
        pt.size = 0.4) +     # Adjust point size for visualization
  ggtitle("UMAP - Cell Cycle Phases") +
  theme_minimal() +
  theme(legend.position = "right")  # Position the legend on the right


#Even though the code was generated for each cell cycle phase separately, it was not very helpful in terms of interpretation.

# Generate UMAP for G1 phase
DimPlot(merged_obj_filtered_15PCs, 
        reduction = "umap", 
        cells.highlight = WhichCells(merged_obj_filtered_15PCs, expression = Phase == "G1"),
        label = FALSE,        # Turn off cluster labels
        pt.size = 0.4,        # Adjust point size
        cols.highlight = c("#F8766D"),  # Default pink (minimal theme)
        cols = "lightgray") + # Light gray for non-highlighted cells
  ggtitle("UMAP - G1 Phase Cells") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend

# Generate UMAP for S phase
DimPlot(merged_obj_filtered_15PCs, 
        reduction = "umap", 
        cells.highlight = WhichCells(merged_obj_filtered_15PCs, expression = Phase == "S"),
        label = FALSE,        # Turn off cluster labels
        pt.size = 0.4,        # Adjust point size
        cols.highlight = c("#00BFC4"),  # Default blue (minimal theme)
        cols = "lightgray") + # Light gray for non-highlighted cells
  ggtitle("UMAP - S Phase Cells") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend

# Generate UMAP for G2M phase
DimPlot(merged_obj_filtered_15PCs, 
        reduction = "umap", 
        cells.highlight = WhichCells(merged_obj_filtered_15PCs, expression = Phase == "G2M"),
        label = FALSE,        # Turn off cluster labels
        pt.size = 0.4,        # Adjust point size
        cols.highlight = c("#7CAE00"),  # Default green (minimal theme)
        cols = "lightgray") + # Light gray for non-highlighted cells
  ggtitle("UMAP - G2M Phase Cells") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend







##########################################################################################################
#################### DIFFERENT VISUALIZATION METHODS TO BETTER UNDERSTAND THE RESULTS ####################
##########################################################################################################



#################### STACKED BARPLOT #######################################

# Normalize by condition
DimPlot(merged_obj_filtered_15PCs, split.by = "orig.ident")

# Tabulate cells by preliminary annotation and orig.ident (acting as condition)
x <- as.data.frame(table(merged_obj_filtered_15PCs$cluster_annotation, 
                         merged_obj_filtered_15PCs$condition))

# Filter out any unwanted condition, if needed (e.g., "NonImplanted" in the original example)
# Uncomment and modify the line below if you need to filter any specific condition
# x <- x %>% filter(Var2 != "Condition_to_exclude")

# Rename columns for clarity
colnames(x) <- c("Preliminary_annotation", "condition", "Frequency")

# Calculate the total cells by condition and obtain the relative frequency
x <- x %>% 
  group_by(condition) %>% 
  mutate(TotalCells = sum(Frequency),
         RelativeFrequency = Frequency / TotalCells)

# Create the stacked bar plot with normalized proportions
ggplot(x, aes(x = condition, y = RelativeFrequency, fill = Preliminary_annotation)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Condition", y = "Relative Frequency", fill = "Preliminary Annotation") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# Create a contingency table for the analysis
contingency_table <- as.data.frame.matrix(
  table(merged_obj_filtered_15PCs$cluster_annotation, merged_obj_filtered_15PCs$condition)
)

# Add row names as a column for better readability
contingency_table <- contingency_table %>%
  rownames_to_column(var = "Preliminary_annotation")

# Perform chi-squared test for each annotation
statistics_results <- contingency_table %>%
  rowwise() %>%
  mutate(
    Chi_sq_p_value = chisq.test(c(melanoma, nevi))$p.value
  ) %>%
  ungroup()

# Adjust p-values for multiple testing using Bonferroni correction
statistics_results <- statistics_results %>%
  mutate(Adjusted_p_value = p.adjust(Chi_sq_p_value, method = "bonferroni"))

# Save the statistics table
write.csv2(statistics_results, "Cell_Type_Statistics.csv", row.names = FALSE)

# Display the results
statistics_results







######################################################
#################### BY CONDITION ####################
######################################################


#################### UMAP BY CONDITION ####################


# Update the meta data in merged_obj_filtered_15PCs
merged_obj_filtered_15PCs@meta.data <- merged_obj_filtered_15PCs@meta.data %>%
  mutate(condition = case_when(
    Sample %in% c("HG2019_2182", "HG2019_2532", "HG2019_2753", "HG2019_4830", "HG2019_867", "N3") ~ "Melanoma",
    Sample %in% c("HG2022_1844", "HG2022_3845", "N1") ~ "Nevi",
    TRUE ~ NA_character_
  ))

# Verify the result
print(merged_obj_filtered_15PCs@meta.data)



# UMAP plot with cells colored by condition
DimPlot(merged_obj_filtered_15PCs, group.by = "condition", reduction = "umap") +
  ggtitle("UMAP Plot Colored by Condition") +
  theme_minimal()




#################### VIOLIN PLOT BY CONDITION #######################################


# Violin plot for quality metrics by condition
VlnPlot(merged_obj_filtered_15PCs, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), 
        group.by = "condition", ncol = 3, pt.size = 0) +
  ggtitle("Violin Plots of QC Metrics by Condition")



#################### DOT PLOT BY CONDITION #######################################


# Dot plot of marker genes by condition
DotPlot(merged_obj_filtered_15PCs, features = c("PMEL", "MLANA", "MITF"), group.by = "condition") +
  labs(title = "Dot Plot of Marker Genes by Condition") +
  theme_minimal()


#################### HEATMAP BY CONDITION TOP 20 GENES #######################################


# Assuming `top20_genes` contains the top marker genes
DoHeatmap(merged_obj_filtered_15PCs, features = top20_genes$gene, group.by = "condition") +
  ggtitle("Heatmap of Top Marker Genes by Condition")








#########################################################
############## DE ANALYSIS ##############################
#########################################################

#This is to get all the differentially expressed genes together (not separated by clusters).

# Normalize RNA assay, if not already done
DefaultAssay(merged_obj_filtered_15PCs) <- "SCT"

# Set "condition" as the active identity class
Idents(merged_obj_filtered_15PCs) <- "condition"

# List to store DE results
de_results_list <- list()

# Start timer
start_time <- Sys.time()

# Differential Expression Analysis between "Nevi" and "Melanoma"
tryCatch({
  print("Processing differential expression analysis between Nevi and Melanoma")
  
  # Perform DE analysis
  markers <- FindMarkers(merged_obj_filtered_15PCs, ident.1 = "Nevi", ident.2 = "Melanoma", 
                         group.by = "condition", assay = "SCT", slot = "data")
  
  # Only proceed if FindMarkers was successful
  if (!is.null(markers)) {
    # Add FDR and gene column
    markers$FDR <- p.adjust(markers$p_val, method = 'fdr')
    markers$gene <- rownames(markers)
    
    # Store the result in the list
    de_results_list[["Nevi_vs_Melanoma"]] <- markers
  } else {
    # If FindMarkers failed, store a message or an empty data frame
    de_results_list[["Nevi_vs_Melanoma"]] <- data.frame(
      message = "FindMarkers failed for Nevi vs Melanoma comparison"
    )
  }
}, error = function(e) {
  message("Error in FindMarkers: ", e$message)
})

# End timer
end_time <- Sys.time()
print(end_time - start_time)

# Create results directory if it doesn't exist
output_directory <- "~/Single_Cell_Patients_Levesque/DE_Results"  # Adjust this to your directory path
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Save the DE results as a CSV file
write.csv2(de_results_list[["Nevi_vs_Melanoma"]], 
           file = file.path(output_directory, "DE_Results_Nevi_vs_Melanoma.csv"),
           row.names = FALSE)
gc()  # Free up memory






###############################################################################################
################################# DE ANALYSIS FOR EACH CELL TYPE ##############################
###############################################################################################

###### FOR EACH cell type!!
#The ident.1 = condition of interest
#The ident.1 is the condition that you compare it with

###################In this case the fold change of 2x (for example) means that x gene is upregulated 2 times in melanomas compared to the nevi###################

# Normalizar los datos de la muestra
merged_obj_filtered_15PCs <- NormalizeData(merged_obj_filtered_15PCs)


# Crear lista para almacenar resultados de cada cluster
merged_obj_filtered_15PCs_list <- list()

# Bucle para comparar cada cluster entre las condiciones "Melanoma" y "Nevi"
start_time <- Sys.time()

for (cell_type in unique(merged_obj_filtered_15PCs$cluster_annotation)) {
  # Filtrar el objeto para el cluster actual
  comparison_subset <- subset(merged_obj_filtered_15PCs, cluster_annotation == cell_type)
  print(paste("Processing cell_type:", cell_type))  # Progreso
  # Realizar el análisis de expresión diferencial entre "Melanoma" y "Nevi"
  markers <- tryCatch({
    FindMarkers(comparison_subset, ident.1 = "Melanoma", ident.2 = "Nevi", 
                group.by = "condition", assay = "originalexp", slot = "data")
  }, error = function(e) {
    message("Error en FindMarkers: ", e$message)
    return(NULL)
  })
  # Solo proceder si FindMarkers fue exitoso
  if (!is.null(markers)) {
    # Añadir columnas de FDR, cluster y gene
    markers$FDR <- p.adjust(markers$p_val, method = 'fdr')
    markers$cluster <- cell_type
    markers$gene <- rownames(markers)
    # Almacenar los resultados en la lista
    merged_obj_filtered_15PCs_list[[as.character(cell_type)]] <- markers
  } else {
    # Si FindMarkers falla, almacenar un mensaje o un data frame vacío
    merged_obj_filtered_15PCs_list[[as.character(cell_type)]] <- data.frame(
      message = "FindMarkers failed for this cluster"
    )
  }
}

end_time <- Sys.time()
print(end_time - start_time)

#Guarda en sheets todos los que haya en la lista
library(writexl)
#Change one name because it gives error
names(merged_obj_filtered_15PCs_list) <- gsub("/", "_", names(merged_obj_filtered_15PCs_list))

write_xlsx(merged_obj_filtered_15PCs_list, "~/Single_Cell_Patients_Levesque/DE_CELL_TYPES_BY_CONDITION.xlsx")
