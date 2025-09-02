```R

# -----------------------------
# 1. Subset: Select B cells only
# -----------------------------
Sobj_BCELL <- subset(Sobj_SINGLET_FIL, subset = Monaco_Cluster == "B cells")

# -----------------------------
# 2. QC and Normalization
# -----------------------------
# Calculate ribosomal gene percentage (RPS / RPL genes)
Sobj_BCELL$RbPercent <- PercentageFeatureSet(Sobj_BCELL, pattern = "^RP[SL]")

# Normalize data with SCTransform (glmGamPoi method) 
# Regress out ribosomal content (RbPercent)
Sobj_BCELL <- SCTransform(
  Sobj_BCELL, 
  method = "glmGamPoi", 
  verbose = TRUE, 
  vars.to.regress = c("RbPercent")
)                          

# -----------------------------
# 3. Dimensionality Reduction
# -----------------------------
# Run PCA using variable features
Sobj_BCELL <- RunPCA(
  Sobj_BCELL, 
  features = VariableFeatures(Sobj_BCELL)
)  

# Determine the number of informative PCs
ElbowPlot(Sobj_BCELL, ndims = 30)
y <- 7  # Number of PCs to use

# -----------------------------
# 4. Clustering and Embedding
# -----------------------------
# Build KNN graph and perform clustering (Leiden algorithm)
Sobj_BCELL <- FindNeighbors(Sobj_BCELL, dims = 1:y)
Sobj_BCELL <- FindClusters(
  Sobj_BCELL, 
  random.seed = 1, 
  resolution = 0.1, 
  algorithm = 4
)

# Run UMAP for visualization
Sobj_BCELL <- RunUMAP(Sobj_BCELL, dims = 1:y)

# UMAP plot colored by Seurat clusters
DimPlot(Sobj_BCELL, group.by = "seurat_clusters", pt.size = 1) +
  ggtitle("UMAP - Color: seurat_clusters") +
  theme_minimal()

# -----------------------------
# 5. Cluster Annotation (SingleR)
# -----------------------------
# Load Monaco Immune reference
ref <- celldex::MonacoImmuneData()

# Extract normalized data from SCT assay
sce_data <- GetAssayData(Sobj_BCELL, assay = "SCT", layer = "data")

# Cluster identities
clusters <- Sobj_BCELL$seurat_clusters              

# Run SingleR annotation
Monaco_cluster_BCELL <- SingleR(
  test = sce_data,
  ref = ref,
  labels = ref$label.fine,
  clusters = clusters
)

# Check assigned labels
Monaco_cluster_BCELL$labels

# -----------------------------
# 6. Manual Cluster Annotation (example mapping)
# -----------------------------
# Assign biological labels to Seurat clusters
Sobj_BCELL$Monaco_Cluster_BCells <- case_when(
  Sobj_BCELL$seurat_clusters == 1 ~ "Switched memory B cells",
  Sobj_BCELL$seurat_clusters == 2 ~ "Naive B cells",
  Sobj_BCELL$seurat_clusters == 3 ~ "Non-switched memory B cells",
  Sobj_BCELL$seurat_clusters == 4 ~ "Plasmablasts",
  TRUE ~ "Other group"
)
