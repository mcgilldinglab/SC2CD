# SC2CD
SC2CD: Co-optimizes Spatial Clustering and Cell Decomposition through Bidirectional Information Flow in Single-Cell Spatial Transcriptomics

![figure1](https://github.com/user-attachments/assets/2d9ae649-d3d4-42de-921c-db68ab79c64a)

## Introduction
Single-cell spatial transcriptomics enables the study of cellular organization. However, computational techniques for spatial clustering and cell-type decomposition are often required, as spot-based platforms such as Visium lack single-cell resolution. Although these two tasks are inherently interconnected, current approaches typically address them separately, limiting information exchange and reducing overall accuracy. We introduce SC2CD, an iterative framework that co-optimizes spatial clustering and cell decomposition through bidirectional information sharing. SC2CD integrates deep learning-based spatial clustering with matrix factorization-based decomposition, dynamically refining spatial graphs using decomposition-informed similarities while leveraging clustering results to improve local deconvolution. Unlike conventional methods that treat these tasks in isolation, SC2CD supports iterative refinement, enabling mutual enhancement. SC2CD consistently outperforms existing approaches, producing more precise spatial clusters and higher-resolution cell-type identification. By unifying spatial clustering and cell-type deconvolution within a co-optimization framework, SC2CD improves spatial transcriptomic resolution and facilitates a more comprehensive understanding of complex tissue architecture.

## Implementation

Step 1: Create a conda environment with rpy2 installed
```
conda create -n scanpy_env
conda activate scanpy_env
conda install numpy=1.19
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg
pip install scanpy==1.8.1

conda create -n scanpy_env python=3.8
conda activate scanpy_env
conda install -c conda-forge r-base
conda install -c conda-forge rpy2
which R
export R_HOME=/path/to/R
export LD_LIBRARY_PATH=$R_HOME/lib:$LD_LIBRARY_PATH

Open your Python, check whether rpy2 works

import rpy2.robjects as robjects
robjects.r('''
   print("Hello from R")
''')
```

Step 2: Download and save the Python package 'GraphST_new', and the R package 'R' in the same folder on your computer

Step 3: Add the data (available in **Need a link**) and 'SC2CD_Tutorial.ipynb' to the folder in Step 2.

The tutorial is working now!

