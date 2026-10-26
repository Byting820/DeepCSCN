# Reproducible DeepCSCN environment

The R package contains the GeneCluster Python source under `inst/python/` and
calls it in-process with `reticulate`. No absolute data path, repository path,
or CUDA device number is embedded in the package.

## Conda: portable CPU environment

```bash
conda env create -f environment.yml
conda activate deepcscn
R CMD INSTALL .
Rscript reproducibility/verify_environment.R
Rscript reproducibility/smoke_test.R
```

The CPU environment is the reference environment for reproducibility because
it does not depend on a host CUDA driver. The package defaults to
`device = "auto"`; set `device = "cpu"` when an exact CPU run is required.

## Conda: NVIDIA GPU environment

The GPU specification pins PyTorch 2.2.2 and CUDA runtime 12.1. A compatible
NVIDIA driver must be installed on the host.

```bash
conda env create -f environment-gpu.yml
conda activate deepcscn-gpu
R CMD INSTALL .
Rscript reproducibility/verify_environment.R
```

Inside R, either activate the environment before starting R or select it before
any Python import:

```r
library(DeepCSCN)
configure_deepcscn_python(condaenv = "deepcscn-gpu")
check_deepcscn_python()
```

## Docker

Build and run the CPU reference image:

```bash
docker build -t deepcscn:0.2.0 .
docker run --rm -it -v "$PWD/work:/work" deepcscn:0.2.0
```

Build and run the GPU image:

```bash
docker build -f Dockerfile.gpu -t deepcscn:0.2.0-gpu .
docker run --rm -it --gpus all -v "$PWD/work:/work" deepcscn:0.2.0-gpu
```

## Files recorded for each GeneCluster run

`run_genecluster()` creates a self-contained output directory containing:

- `input_expression.tsv`: the exact gene-by-cell matrix sent by R;
- `features.tsv` and `features.npy`: the final 512-dimensional embeddings;
- `clusters.tsv`: final DeepCluster pseudo-clusters;
- `training_history.tsv`: epoch loss, K-means inertia, silhouette, and learning rate;
- `model.pt`: PyTorch state, parameter values, gene/cell order, and epoch count;
- `run_metadata.json`: all hyperparameters, package versions, resolved device,
  timestamps, and the SHA-256 checksum of the input matrix.

Random seeds are set for R, Python, NumPy, scikit-learn, and PyTorch. PyTorch
deterministic algorithms are requested by default. In deterministic mode, PyTorch
and BLAS/OpenMP numerical reductions are limited to one thread so that tiny
parallel-reduction differences are not amplified by DeepCluster pseudo-labels. Exact floating-point values
can still vary across CPU instruction sets, GPU models, CUDA drivers, and
PyTorch builds; the metadata file makes those sources of variation auditable.
