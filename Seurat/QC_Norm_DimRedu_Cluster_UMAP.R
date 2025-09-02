```R

# -----------------------------
# Quality Control and Filtering
# -----------------------------
# Calculate percentage of mitochondrial genes
Sobj_SINGLET$MtPercent <- PercentageFeatureSet(Sobj_SINGLET, pattern = "^MT-")

# Filter cells based on RNA features, counts, and mitochondrial content
Sobj_SINGLET_FIL <- subset(
  Sobj_SINGLET,
  subset = nFeature_RNA > 400 &
           nFeature_RNA < 5500 &
           nCount_RNA > 1000 &
           nCount_RNA < 40000 &
           MtPercent < 20
)

# -----------------------------
# Set options for large datasets
# -----------------------------
options(future.globals.maxSize = 50 * 1024^3)  # Increase memory limit
options(glmGamPoi.parallel = TRUE)             # Enable parallel processing for SCTransform

# -----------------------------
# Normalization using SCTransform
# -----------------------------
Sobj_SINGLET_FIL <- SCTransform(
  Sobj_SINGLET_FIL,
  method = "glmGamPoi",
  verbose = TRUE
)

# -----------------------------
# Dimensionality Reduction
# -----------------------------
Sobj_SINGLET_FIL <- RunPCA(
  Sobj_SINGLET_FIL,
  features = VariableFeatures(Sobj_SINGLET_FIL)
)

# Plot elbow plot to determine number of PCs
ElbowPlot(Sobj_SINGLET_FIL, ndims = 30)
NPCs <- 10  # Number of PCs to use downstream

# -----------------------------
# Clustering
# -----------------------------
Sobj_SINGLET_FIL <- FindNeighbors(Sobj_SINGLET_FIL, dims = 1:NPCs)
Sobj_SINGLET_FIL <- FindClusters(
  Sobj_SINGLET_FIL,
  random.seed = 1,
  resolution = 1.5,
  algorithm = 4
)

# -----------------------------
# UMAP Visualization
# -----------------------------
Sobj_SINGLET_FIL <- RunUMAP(Sobj_SINGLET_FIL, dims = 1:NPCs)

DimPlot(
  Sobj_SINGLET_FIL,
  group.by = "seurat_clusters",
  pt.size = 1,
  label = TRUE
) +
  ggtitle("UMAP - Color: Seurat Cluster") +
  theme_minimal()
