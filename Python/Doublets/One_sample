```python

# -----------------------------
# Dependencies
# -----------------------------
!pip install scanpy matplotlib scipy numpy pandas scrublet

# -----------------------------
# Import libraries
# -----------------------------
import numpy as np
import scipy.io
import scrublet as scr
import pandas as pd
from scipy.io import mmwrite

# -----------------------------
# Load data
# -----------------------------
counts_matrix = scipy.io.mmread("path/to/matrix.mtx").T.tocsc()  # Load and transpose sparse count matrix
barcodes = pd.read_csv("path/to/barcodes.tsv", header=None)[0]   # Load cell barcodes

# -----------------------------
# Doublet detection using Scrublet
# -----------------------------
scrub = scr.Scrublet(counts_matrix)                              # Initialize Scrublet object
doublet_scores, predicted_doublets = scrub.scrub_doublets()       # Predict doublets

# Plot histogram of doublet scores
scrub.plot_histogram()

# -----------------------------
# Count doublets and singlets
# -----------------------------
n_doublets = np.sum(predicted_doublets)
n_singlets = np.sum(~predicted_doublets)
print(f"Detected {n_doublets} doublets and {n_singlets} singlets")

# -----------------------------
# Filter out doublets
# -----------------------------
filtered_matrix = counts_matrix[~predicted_doublets, :]           # Keep only singlet cells
filtered_barcodes = barcodes[~predicted_doublets]

# -----------------------------
# Save filtered data
# -----------------------------
mmwrite("filtered_matrix.mtx", filtered_matrix.T)                # Save filtered matrix
filtered_barcodes.to_csv("filtered_barcodes.tsv", index=False, header=False)  # Save filtered barcodes
