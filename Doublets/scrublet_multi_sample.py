```python

# -----------------------------
# Dependencies
# -----------------------------
!pip install scanpy matplotlib scipy numpy pandas scrublet umap-learn

# -----------------------------
# Import libraries
# -----------------------------
import os
import gzip
import numpy as np
import pandas as pd
import scipy.io
import scrublet as scr
import matplotlib.pyplot as plt
from scipy.io import mmwrite
from umap import UMAP

# -----------------------------
# Define input/output directories (replace with your paths)
# -----------------------------
input_root = r"/path/to/your/samples"         # Folder containing all sample folders
output_root = r"/path/to/save/results"        # Folder to save all results
os.makedirs(output_root, exist_ok=True)       # Create output folder if it doesn't exist

# -----------------------------
# List to store summary data
# -----------------------------
summary_data = []  # Will store sample name, number of doublets, number of singlets

# -----------------------------
# Function to load a 10X Genomics matrix from a folder
# -----------------------------
def load_10x_matrix(matrix_dir):
    """Load a 10X sparse count matrix and metadata from a folder."""
    
    # Load and transpose sparse matrix
    matrix = scipy.io.mmread(os.path.join(matrix_dir, 'matrix.mtx.gz')).T.tocsc()
    
    # Load cell barcodes
    barcodes_path = os.path.join(matrix_dir, 'barcodes.tsv.gz')
    with gzip.open(barcodes_path, 'rt') as f:
        barcodes = np.array([line.strip() for line in f])
    
    # Load gene/features information
    features_path = os.path.join(matrix_dir, 'features.tsv.gz') \
        if os.path.exists(os.path.join(matrix_dir, 'features.tsv.gz')) \
        else os.path.join(matrix_dir, 'genes.tsv.gz')
    with gzip.open(features_path, 'rt') as f:
        features = [line.strip().split('\t') for line in f]
    feature_ids = [f[0] for f in features]
    feature_names = [f[1] if len(f) > 1 else f[0] for f in features]
    
    return matrix, barcodes, feature_names, feature_ids

# -----------------------------
# Function to save a matrix in 10X-like format
# -----------------------------
def save_matrix(out_dir, matrix, barcodes, gene_ids, gene_names):
    """Save a filtered 10X-like count matrix with barcodes and gene names."""
    os.makedirs(out_dir, exist_ok=True)
    mmwrite(os.path.join(out_dir, "matrix.mtx"), matrix.T)
    with open(os.path.join(out_dir, "barcodes.tsv"), 'w') as f:
        f.writelines([bc + '\n' for bc in barcodes])
    with open(os.path.join(out_dir, "genes.tsv"), 'w') as f:
        for i in range(len(gene_names)):
            f.write(f"{gene_ids[i]}\t{gene_names[i]}\n")

# -----------------------------
# Function to retrieve doublet rate for a sample
# -----------------------------
def get_doublet_rate(sample_name, samples, rates):
    """Return the expected doublet rate for a given sample."""
    if not samples or not rates:
        raise ValueError("Sample list or rate list is empty.")
    if sample_name not in samples:
        raise ValueError(f"Sample {sample_name} not found in the list.")
    idx = samples.index(sample_name)
    return rates[idx]

# -----------------------------
# Example sample names and corresponding expected doublet rates
# Replace with your own sample list and rates
# -----------------------------
samples = ["Sample1", "Sample2", "Sample3"]
expected_rates = [0.08, 0.12, 0.05]  # Example doublet rates

# -----------------------------
# Loop through each sample folder
# -----------------------------
for folder in os.listdir(input_root):
    sample_path = os.path.join(input_root, folder, "filtered_feature_bc_matrix")
    matrix_file = os.path.join(sample_path, "matrix.mtx.gz")
    if not os.path.exists(matrix_file):
        continue  # Skip if matrix file does not exist
    
    print(f"Processing sample: {folder}")
    
    # Load matrix and metadata
    matrix, barcodes, gene_names, gene_ids = load_10x_matrix(sample_path)
    
    # Get expected doublet rate
    doublet_rate = get_doublet_rate(folder, samples, expected_rates)
    
    # Run Scrublet for doublet detection
    scrub = scr.Scrublet(matrix, expected_doublet_rate=doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_cells=3, n_prin_comps=30)
    
    # -----------------------------
    # Create output folder for sample
    # -----------------------------
    sample_output_dir = os.path.join(output_root, folder)
    os.makedirs(sample_output_dir, exist_ok=True)
    
    # -----------------------------
    # Plot and save histogram of doublet scores
    # -----------------------------
    hist_path = os.path.join(sample_output_dir, f"{folder}_histogram.png")
    plt.figure()
    scrub.plot_histogram()
    plt.title(f"{folder} - Scrublet Histogram")
    plt.savefig(hist_path)
    plt.close()
    
    # -----------------------------
    # Save doublet call table
    # -----------------------------
    result_df = pd.DataFrame({
        'barcode': barcodes,
        'doublet_score': doublet_scores,
        'predicted_doublet': predicted_doublets
    })
    output_csv_path = os.path.join(sample_output_dir, f"{folder}_doublet_calls.csv")
    result_df.to_csv(output_csv_path, index=False)
    
    # -----------------------------
    # Add sample prefix to barcodes
    # -----------------------------
    barcodes_with_prefix = [f"{folder}_{bc}" for bc in barcodes]
    
    # -----------------------------
    # Split matrix into singlets and doublets
    # -----------------------------
    matrix_no_doublets = matrix[~predicted_doublets, :]
    barcodes_no_doublets = np.array(barcodes_with_prefix)[~predicted_doublets]
    matrix_only_doublets = matrix[predicted_doublets, :]
    barcodes_only_doublets = np.array(barcodes_with_prefix)[predicted_doublets]
    
    # Save filtered matrices
    save_matrix(os.path.join(sample_output_dir, "filtered_no_doublets"), 
                matrix_no_doublets, barcodes_no_doublets, gene_ids, gene_names)
    save_matrix(os.path.join(sample_output_dir, "filtered_only_doublets"), 
                matrix_only_doublets, barcodes_only_doublets, gene_ids, gene_names)
    
    # -----------------------------
    # Collect summary data
    # -----------------------------
    num_doublets = np.sum(predicted_doublets)
    num_singlets = len(predicted_doublets) - num_doublets
    summary_data.append([folder, num_doublets, num_singlets])
    
    # -----------------------------
    # UMAP embedding based on Scrublet manifold
    # -----------------------------
    embedding = UMAP(n_neighbors=10, min_dist=0.3).fit_transform(scrub.manifold_obs_)
    
    # Plot UMAP colored by doublet score
    umap_path = os.path.join(sample_output_dir, f"{folder}_UMAP.png")
    plt.figure(figsize=(8, 8))
    plt.scatter(embedding[:, 0], embedding[:, 1], c=scrub.doublet_scores_obs_, cmap='viridis', s=5)
    plt.colorbar(label='Doublet Score')
    plt.title(f'{folder} - UMAP Colored by Doublet Score')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.tight_layout()
    plt.savefig(umap_path, dpi=300)
    plt.close()
    
    print(f"UMAP plot saved at: {umap_path}")
    print(f"✓ Sample {folder} processed.\n")

# -----------------------------
# Save summary CSV for all samples
# -----------------------------
summary_df = pd.DataFrame(summary_data, columns=['Sample', 'Num_Doublets', 'Num_Singlets'])
summary_csv_path = os.path.join(output_root, 'scrublet_summary.csv')
summary_df.to_csv(summary_csv_path, index=False)

print("All Scrublet analyses completed.")
