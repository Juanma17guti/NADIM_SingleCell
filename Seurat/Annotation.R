```R

# -----------------------------
# SingleR Annotation
# -----------------------------
ref <- celldex::MonacoImmuneData()
sce_data <- GetAssayData(Sobj_SINGLET_FIL, assay = "SCT", layer = "data")
clusters <- Sobj_SINGLET_FIL$seurat_clusters              

Monaco_cluster <- SingleR(
  test = sce_data,
  ref = ref,
  labels = ref$label.main,
  clusters = clusters
)

# Add automatic SingleR annotations
Sobj_SINGLET_FIL$SingleR_labels <- Monaco_cluster$labels[match(
  Sobj_SINGLET_FIL$seurat_clusters,
  rownames(Monaco_cluster)
)]

# -----------------------------
# Manual curation (if needed)
# -----------------------------
Sobj_SINGLET_FIL$Monaco_Cluster <- case_when(
  Sobj_SINGLET_FIL$seurat_clusters == 29 ~ "Dendritic cells",
  Sobj_SINGLET_FIL$seurat_clusters == 34 ~ "Progenitors",
  TRUE ~ Sobj_SINGLET_FIL$SingleR_labels  # fallback: keep SingleR labels
)
