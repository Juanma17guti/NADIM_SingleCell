```R

# -----------------------------
# Load libraries
# -----------------------------
library(Seurat)
library(stringr)

# -----------------------------
# Initialize list to store Seurat objects for all samples
# -----------------------------
seurat_list <- list()

# -----------------------------
# Example metadata vectors (replace with your own sample metadata)
# -----------------------------
SAM <- c("Sample1", "Sample2", "Sample3")              # Sample names
NaM <- c("Patient1", "Patient2", "Patient3")          # Patient IDs
RES <- c("CPR", "NCPR", "CPR")                        # Response status
SEX <- c("M", "F", "M")                               # Sex
PRO <- c("0", "1", "0")                               # Progression status

# -----------------------------
# Define results directory (replace with your path)
# -----------------------------
results_dir <- "/path/to/scrublet/filtered_samples"
sample_folders <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)

# -----------------------------
# Loop through each sample folder
# -----------------------------
for (folder in sample_folders) {
  sample_name <- basename(folder)
  
  # Define directories for singlets and doublets
  singlet_dir <- file.path(folder, paste0(sample_name, "_filtered_matrix_no_doublets"))
  doublet_dir <- file.path(folder, paste0(sample_name, "_filtered_matrix_only_doublets"))
  
  # -----------------------------
  # Load singlet cells and create Seurat object
  # -----------------------------
  if (dir.exists(singlet_dir)) {
    singlet_counts <- Read10X(data.dir = singlet_dir)
    singlet_seurat <- CreateSeuratObject(singlet_counts, project = sample_name)
    singlet_seurat$patient <- sample_name
    singlet_seurat$doublet_status <- "singlet"
  }
  
  # -----------------------------
  # Load doublet cells and create Seurat object
  # -----------------------------
  if (dir.exists(doublet_dir)) {
    doublet_counts <- Read10X(data.dir = doublet_dir)
    doublet_seurat <- CreateSeuratObject(doublet_counts, project = sample_name)
    doublet_seurat$patient <- sample_name
    doublet_seurat$doublet_status <- "doublet"
  }
  
  # -----------------------------
  # Merge singlet and doublet Seurat objects for this sample
  # -----------------------------
  combined_sample <- merge(singlet_seurat, y = doublet_seurat)
  
  # -----------------------------
  # Add example metadata (visit, patient ID, response, progression, sex)
  # Replace with your own metadata vectors
  # -----------------------------
  EXTRA <- str_sub(sample_name, -3)      # Example: extract last 3 chars of sample name
  combined_sample$Visit <- EXTRA
  
  POSITION <- which(sapply(SAM, function(x) x == sample_name))
  combined_sample$ExID <- NaM[POSITION]
  combined_sample$Response <- RES[POSITION]
  combined_sample$Progression <- PRO[POSITION]
  combined_sample$Sex <- SEX[POSITION]
  
  # Store the combined Seurat object in the list
  seurat_list[[sample_name]] <- combined_sample
}

# -----------------------------
# Merge all samples into a single Seurat object
# -----------------------------
combined_all <- merge(seurat_list[[1]], y = seurat_list[-1], 
                      add.cell.ids = names(seurat_list), project = "ExampleProject")

# -----------------------------
# Subset only singlet cells
# -----------------------------
Sobj_SINGLET <- subset(combined_all, subset = doublet_status == "singlet")
