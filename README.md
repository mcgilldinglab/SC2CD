# SC2CD
SC2CD: Co-optimizes Spatial Clustering and Cell Decomposition through Bidirectional Information Flow in Single-Cell Spatial Transcriptomics

![figure1](https://github.com/user-attachments/assets/2d9ae649-d3d4-42de-921c-db68ab79c64a)

## Introduction
Single-cell spatial transcriptomics enables the study of cellular organization. However, computational techniques for spatial clustering and cell-type decomposition are often required, as spot-based platforms such as Visium lack single-cell resolution. Although these two tasks are inherently interconnected, current approaches typically address them separately, limiting information exchange and reducing overall accuracy. We introduce SC2CD, an iterative framework that co-optimizes spatial clustering and cell decomposition through bidirectional information sharing. SC2CD integrates deep learning-based spatial clustering with matrix factorization-based decomposition, dynamically refining spatial graphs using decomposition-informed similarities while leveraging clustering results to improve local deconvolution. Unlike conventional methods that treat these tasks in isolation, SC2CD supports iterative refinement, enabling mutual enhancement. SC2CD consistently outperforms existing approaches, producing more precise spatial clusters and higher-resolution cell-type identification. By unifying spatial clustering and cell-type deconvolution within a co-optimization framework, SC2CD improves spatial transcriptomic resolution and facilitates a more comprehensive understanding of complex tissue architecture.

## Implementation

Step 1: Create a conda environment with necessary packages installed:

Download setup_sc2cd_env.sh and run the following in your terminal:

```
# Default installation
bash setup_sc2cd_env.sh

# Optional: Create an environment with a custom name
bash setup_sc2cd_env.sh <your_env_name>
```

Step 2: Download and save the Python package 'SC2CD_py', and the R package 'R' in a folder on your computer

Step 3: Add the data ([https://doi.org/10.5281/zenodo.17601120](https://zenodo.org/records/17601120)) and 'SC2CD_Tutorial.ipynb' to the folder in Step 2.

The tutorial is working now!

