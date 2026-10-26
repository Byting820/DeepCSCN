"""Reproducible Python backend for the DeepCSCN GeneCluster model.

This module is imported by the DeepCSCN R package through reticulate.  It keeps
GeneCluster's CNN, spatial-pyramid-pooling (SPP), DeepCluster pseudo-label, and
focal-loss workflow while removing machine-specific paths and fixed CUDA device
assignments from the original research script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os

# Set deterministic numerical-library defaults before importing NumPy/PyTorch.
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import platform
import random
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import sklearn
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler
from threadpoolctl import threadpool_limits


@dataclass(frozen=True)
class TrainConfig:
    n_clusters: int = 15
    epochs: int = 100
    batch_size: int = 32
    learning_rate: float = 0.01
    momentum: float = 0.9
    weight_decay: float = 1e-5
    focal_gamma: float = 2.0
    seed: int = 0
    device: str = "auto"
    num_workers: int = 0
    patience: int = 10
    deterministic: bool = True


class GeneExpressionDataset(Dataset):
    """Genes are samples; ordered cell measurements are the 1-D input."""

    def __init__(self, expression: np.ndarray, labels: Optional[np.ndarray] = None):
        self.expression = torch.as_tensor(expression[:, np.newaxis, :], dtype=torch.float32)
        self.labels = None if labels is None else torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.expression.shape[0])

    def __getitem__(self, index: int):
        if self.labels is None:
            return self.expression[index], index
        return self.expression[index], self.labels[index]


class SpatialPyramidPool1D(nn.Module):
    """One-dimensional SPP used by the original GeneCluster architecture."""

    def __init__(self, levels: Sequence[int] = (1, 3, 6, 8, 10, 16)):
        super().__init__()
        self.levels = tuple(int(level) for level in levels)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        length = int(x.size(-1))
        pooled: List[torch.Tensor] = []
        for level in self.levels:
            window = int(np.ceil(length / level))
            target_length = window * level
            pad = target_length - length
            left = pad // 2
            right = pad - left
            padded = F.pad(x, (left, right), mode="constant", value=0.0)
            value = F.max_pool1d(padded, kernel_size=window, stride=window)
            pooled.append(value.flatten(start_dim=1))
        return torch.cat(pooled, dim=1)


class GeneClusterCNN(nn.Module):
    """CNN-SPP encoder with the 512-dimensional GeneCluster representation."""

    def __init__(self, input_channels: int = 1, embedding_dim: int = 512):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv1d(input_channels, 64, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(64),
            nn.ReLU(inplace=True),
            nn.Conv1d(64, 128, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(128),
            nn.ReLU(inplace=True),
            nn.Conv1d(128, 256, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm1d(256),
            nn.ReLU(inplace=True),
        )
        self.spp = SpatialPyramidPool1D()
        self.embedding = nn.Sequential(
            nn.Linear(256 * sum(self.spp.levels), embedding_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.2),
            nn.Linear(embedding_dim, embedding_dim),
            nn.ReLU(inplace=True),
        )
        self.top_layer: Optional[nn.Linear] = None
        self._initialize_weights()

    def _initialize_weights(self) -> None:
        for module in self.modules():
            if isinstance(module, nn.Conv1d):
                nn.init.kaiming_normal_(module.weight, mode="fan_out", nonlinearity="relu")
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.BatchNorm1d):
                nn.init.ones_(module.weight)
                nn.init.zeros_(module.bias)
            elif isinstance(module, nn.Linear):
                nn.init.normal_(module.weight, mean=0.0, std=0.01)
                nn.init.zeros_(module.bias)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.embedding(self.spp(self.features(x)))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        embedding = self.encode(x)
        if self.top_layer is None:
            return embedding
        return self.top_layer(embedding)


class FocalLoss(nn.Module):
    def __init__(self, gamma: float = 2.0):
        super().__init__()
        self.gamma = float(gamma)

    def forward(self, logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        log_probability = F.log_softmax(logits, dim=1)
        probability = torch.exp(log_probability)
        weighted_log_probability = (1.0 - probability).pow(self.gamma) * log_probability
        return F.nll_loss(weighted_log_probability, target)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _set_seed(seed: int, deterministic: bool) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    if deterministic:
        # Multithreaded reductions can vary in the last bits across processes;
        # DeepCluster can amplify those differences through pseudo-labels.
        torch.set_num_threads(1)
        try:
            torch.set_num_interop_threads(1)
        except RuntimeError:
            # PyTorch allows this setting only before inter-op work starts.
            pass
    if hasattr(torch, "use_deterministic_algorithms"):
        try:
            torch.use_deterministic_algorithms(deterministic, warn_only=True)
        except TypeError:
            # PyTorch 1.8 exposes the flag without the warn_only argument.
            torch.use_deterministic_algorithms(deterministic)
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = deterministic
        torch.backends.cudnn.benchmark = not deterministic


def _resolve_device(requested: str) -> torch.device:
    requested = str(requested).strip().lower()
    if requested == "auto":
        return torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if requested.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError(
            "A CUDA device was requested, but torch.cuda.is_available() is False. "
            "Use device='cpu' or install the GPU environment and NVIDIA runtime."
        )
    try:
        device = torch.device(requested)
    except (RuntimeError, ValueError) as exc:
        raise ValueError("device must be 'auto', 'cpu', 'cuda', or 'cuda:<index>'") from exc
    if device.type == "cuda":
        index = 0 if device.index is None else device.index
        if index >= torch.cuda.device_count():
            raise RuntimeError(
                f"CUDA device {index} was requested, but only {torch.cuda.device_count()} device(s) are visible."
            )
    return device


def _read_expression(input_path: Path) -> pd.DataFrame:
    expression = pd.read_csv(input_path, sep="\t", header=0, index_col=0, compression="infer")
    if expression.empty:
        raise ValueError("The input expression matrix is empty.")
    if expression.index.has_duplicates:
        duplicated = expression.index[expression.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Gene names must be unique; duplicated examples: {duplicated}")
    if expression.columns.has_duplicates:
        duplicated = expression.columns[expression.columns.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Cell names must be unique; duplicated examples: {duplicated}")
    try:
        expression = expression.astype(np.float32)
    except (TypeError, ValueError) as exc:
        raise ValueError("All expression values must be numeric.") from exc
    values = expression.to_numpy(copy=False)
    if not np.isfinite(values).all():
        raise ValueError("The expression matrix contains NA, NaN, or infinite values.")
    if expression.shape[0] < 3:
        raise ValueError("At least three genes are required.")
    if expression.shape[1] < 2:
        raise ValueError("At least two cells are required.")
    return expression


def _extract_features(
    model: GeneClusterCNN,
    expression: np.ndarray,
    batch_size: int,
    device: torch.device,
    num_workers: int,
) -> np.ndarray:
    model.eval()
    dataset = GeneExpressionDataset(expression)
    loader = DataLoader(
        dataset,
        batch_size=min(batch_size, len(dataset)),
        shuffle=False,
        num_workers=num_workers,
        pin_memory=device.type == "cuda",
    )
    output = np.empty((len(dataset), 512), dtype=np.float32)
    with torch.no_grad():
        for samples, indices in loader:
            features = model.encode(samples.to(device, non_blocking=True))
            output[indices.numpy()] = features.detach().cpu().numpy()
    return output


def _preprocess_for_clustering(features: np.ndarray, seed: int) -> Tuple[np.ndarray, int]:
    n_components = min(20, features.shape[0] - 1, features.shape[1])
    if n_components < 1:
        raise ValueError("Not enough observations for PCA whitening.")
    with threadpool_limits(limits=1):
        transformed = PCA(n_components=n_components, whiten=True, random_state=seed).fit_transform(features)
        norms = np.linalg.norm(transformed, axis=1, keepdims=True)
        transformed = transformed / np.maximum(norms, np.finfo(np.float32).eps)
    return transformed.astype(np.float32, copy=False), n_components


def _cluster_features(
    features: np.ndarray, n_clusters: int, seed: int
) -> Tuple[np.ndarray, float, float, int]:
    transformed, n_components = _preprocess_for_clustering(features, seed)
    kmeans = KMeans(
        n_clusters=n_clusters,
        init="k-means++",
        n_init=10,
        max_iter=50,
        random_state=seed,
        algorithm="lloyd",
    )
    with threadpool_limits(limits=1):
        labels = kmeans.fit_predict(transformed).astype(np.int64)
        if 1 < np.unique(labels).size < labels.size:
            silhouette = float(silhouette_score(transformed, labels, metric="euclidean"))
        else:
            silhouette = float("nan")
    return labels, float(kmeans.inertia_), silhouette, n_components


def _balanced_sampler(labels: np.ndarray, seed: int) -> WeightedRandomSampler:
    counts = np.bincount(labels)
    weights = np.asarray([1.0 / counts[label] for label in labels], dtype=np.float64)
    generator = torch.Generator()
    generator.manual_seed(int(seed))
    return WeightedRandomSampler(
        weights=torch.as_tensor(weights, dtype=torch.double),
        num_samples=len(labels),
        replacement=True,
        generator=generator,
    )


def _train_epoch(
    model: GeneClusterCNN,
    expression: np.ndarray,
    labels: np.ndarray,
    config: TrainConfig,
    device: torch.device,
    backbone_optimizer: torch.optim.Optimizer,
) -> float:
    model.top_layer = nn.Linear(512, config.n_clusters).to(device)
    nn.init.normal_(model.top_layer.weight, mean=0.0, std=0.01)
    nn.init.zeros_(model.top_layer.bias)
    top_optimizer = torch.optim.SGD(
        model.top_layer.parameters(),
        lr=config.learning_rate,
        momentum=config.momentum,
        weight_decay=config.weight_decay,
    )
    criterion = FocalLoss(gamma=config.focal_gamma)
    dataset = GeneExpressionDataset(expression, labels)
    loader = DataLoader(
        dataset,
        batch_size=min(config.batch_size, len(dataset)),
        sampler=_balanced_sampler(labels, config.seed),
        num_workers=config.num_workers,
        pin_memory=device.type == "cuda",
    )
    model.train()
    total_loss = 0.0
    total_samples = 0
    for samples, targets in loader:
        samples = samples.to(device, non_blocking=True)
        targets = targets.to(device, non_blocking=True)
        backbone_optimizer.zero_grad(set_to_none=True)
        top_optimizer.zero_grad(set_to_none=True)
        logits = model(samples)
        loss = criterion(logits, targets)
        loss.backward()
        backbone_optimizer.step()
        top_optimizer.step()
        batch_n = int(samples.size(0))
        total_loss += float(loss.detach().cpu()) * batch_n
        total_samples += batch_n
    return total_loss / max(total_samples, 1)


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, Path):
        return str(value)
    return value


def train_genecluster(
    input_path: str,
    output_dir: str,
    n_clusters: int = 15,
    epochs: int = 100,
    batch_size: int = 32,
    learning_rate: float = 0.01,
    momentum: float = 0.9,
    weight_decay: float = 1e-5,
    focal_gamma: float = 2.0,
    seed: int = 0,
    device: str = "auto",
    num_workers: int = 0,
    patience: int = 10,
    deterministic: bool = True,
    overwrite: bool = False,
    verbose: bool = True,
) -> Dict[str, object]:
    """Train GeneCluster and write the representation and audit files.

    Parameters are deliberately plain Python types so that reticulate can call
    the function without custom conversion code.
    """

    config = TrainConfig(
        n_clusters=int(n_clusters),
        epochs=int(epochs),
        batch_size=int(batch_size),
        learning_rate=float(learning_rate),
        momentum=float(momentum),
        weight_decay=float(weight_decay),
        focal_gamma=float(focal_gamma),
        seed=int(seed),
        device=str(device),
        num_workers=int(num_workers),
        patience=int(patience),
        deterministic=bool(deterministic),
    )
    if config.epochs < 1:
        raise ValueError("epochs must be at least 1")
    if config.batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    if config.learning_rate <= 0:
        raise ValueError("learning_rate must be positive")
    if config.focal_gamma < 0:
        raise ValueError("focal_gamma must be non-negative")
    if config.num_workers < 0:
        raise ValueError("num_workers must be non-negative")

    input_file = Path(input_path).expanduser().resolve()
    if not input_file.is_file():
        raise FileNotFoundError(f"Expression file not found: {input_file}")
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    expected = [
        destination / "features.tsv",
        destination / "features.npy",
        destination / "clusters.tsv",
        destination / "training_history.tsv",
        destination / "model.pt",
        destination / "run_metadata.json",
    ]
    if not overwrite:
        existing = [str(path) for path in expected if path.exists()]
        if existing:
            raise FileExistsError(
                "Output files already exist. Set overwrite=True to replace them: " + ", ".join(existing)
            )

    started_at = _utc_now()
    started_clock = time.time()
    _set_seed(config.seed, config.deterministic)
    torch_device = _resolve_device(config.device)
    expression_df = _read_expression(input_file)
    expression = expression_df.to_numpy(dtype=np.float32, copy=True)
    if not 2 <= config.n_clusters < expression.shape[0]:
        raise ValueError(
            f"n_clusters must be between 2 and n_genes - 1 ({expression.shape[0] - 1})."
        )

    if verbose:
        print(
            f"DeepCSCN GeneCluster: {expression.shape[0]} genes x {expression.shape[1]} cells; "
            f"device={torch_device}; seed={config.seed}",
            flush=True,
        )

    model = GeneClusterCNN().to(torch_device)
    backbone_parameters: Iterable[torch.nn.Parameter] = list(model.features.parameters()) + list(
        model.embedding.parameters()
    )
    backbone_optimizer = torch.optim.SGD(
        backbone_parameters,
        lr=config.learning_rate,
        momentum=config.momentum,
        weight_decay=config.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.ExponentialLR(backbone_optimizer, gamma=0.9)

    history: List[Dict[str, object]] = []
    best_loss = float("inf")
    bad_epochs = 0
    stopped_early = False

    for epoch in range(config.epochs):
        features = _extract_features(
            model, expression, config.batch_size, torch_device, config.num_workers
        )
        labels, inertia, silhouette, pca_components = _cluster_features(
            features, config.n_clusters, config.seed + epoch
        )
        loss = _train_epoch(
            model, expression, labels, config, torch_device, backbone_optimizer
        )
        current_lr = float(backbone_optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": epoch + 1,
                "loss": loss,
                "kmeans_inertia": inertia,
                "silhouette": silhouette,
                "learning_rate": current_lr,
                "pca_components": pca_components,
            }
        )
        if verbose:
            print(
                f"epoch={epoch + 1:03d} loss={loss:.6f} "
                f"silhouette={silhouette:.4f} lr={current_lr:.6g}",
                flush=True,
            )
        if loss < best_loss - 1e-8:
            best_loss = loss
            bad_epochs = 0
        else:
            bad_epochs += 1
        scheduler.step()
        if config.patience >= 0 and bad_epochs > config.patience:
            stopped_early = True
            if verbose:
                print(f"Early stopping after epoch {epoch + 1}.", flush=True)
            break

    final_features = _extract_features(
        model, expression, config.batch_size, torch_device, config.num_workers
    )
    epochs_run = len(history)
    final_labels, final_inertia, final_silhouette, final_pca_components = _cluster_features(
        final_features, config.n_clusters, config.seed + epochs_run
    )

    feature_columns = [f"feature_{index + 1}" for index in range(final_features.shape[1])]
    feature_table = pd.DataFrame(final_features, index=expression_df.index, columns=feature_columns)
    feature_path = destination / "features.tsv"
    feature_npy_path = destination / "features.npy"
    cluster_path = destination / "clusters.tsv"
    history_path = destination / "training_history.tsv"
    model_path = destination / "model.pt"
    metadata_path = destination / "run_metadata.json"

    feature_table.to_csv(feature_path, sep="\t", index=True, index_label="gene")
    np.save(feature_npy_path, final_features)
    pd.DataFrame(
        {"gene": expression_df.index.astype(str), "cluster": final_labels.astype(int)}
    ).to_csv(cluster_path, sep="\t", index=False)
    pd.DataFrame(history).to_csv(history_path, sep="\t", index=False)
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "config": asdict(config),
            "genes": expression_df.index.astype(str).tolist(),
            "cells": expression_df.columns.astype(str).tolist(),
            "epochs_run": epochs_run,
        },
        model_path,
    )

    metadata: Dict[str, object] = {
        "backend": "DeepCSCN GeneCluster Python backend",
        "backend_version": "1.0.0",
        "started_at_utc": started_at,
        "finished_at_utc": _utc_now(),
        "elapsed_seconds": time.time() - started_clock,
        "input_path": str(input_file),
        "input_sha256": _sha256(input_file),
        "n_genes": int(expression_df.shape[0]),
        "n_cells": int(expression_df.shape[1]),
        "config": asdict(config),
        "resolved_device": str(torch_device),
        "epochs_run": epochs_run,
        "stopped_early": stopped_early,
        "best_training_loss": best_loss,
        "final_kmeans_inertia": final_inertia,
        "final_silhouette": final_silhouette,
        "final_pca_components": final_pca_components,
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "torch": torch.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "cuda_available": bool(torch.cuda.is_available()),
            "torch_cuda": torch.version.cuda,
            "cudnn": torch.backends.cudnn.version() if torch.backends.cudnn.is_available() else None,
        },
        "outputs": {
            "features_tsv": str(feature_path),
            "features_npy": str(feature_npy_path),
            "clusters_tsv": str(cluster_path),
            "training_history_tsv": str(history_path),
            "model_checkpoint": str(model_path),
        },
    }
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(_json_safe(metadata), handle, indent=2, ensure_ascii=False)

    result = {
        "features_path": str(feature_path),
        "features_npy_path": str(feature_npy_path),
        "clusters_path": str(cluster_path),
        "history_path": str(history_path),
        "model_path": str(model_path),
        "metadata_path": str(metadata_path),
        "n_genes": int(expression_df.shape[0]),
        "n_cells": int(expression_df.shape[1]),
        "epochs_run": epochs_run,
        "device": str(torch_device),
        "final_silhouette": final_silhouette,
    }
    return _json_safe(result)


def dependency_versions() -> Dict[str, object]:
    """Return the Python dependency versions used by the active interpreter."""

    return _json_safe(
        {
            "python_executable": sys.executable,
            "python": platform.python_version(),
            "torch": torch.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "cuda_available": bool(torch.cuda.is_available()),
            "torch_cuda": torch.version.cuda,
            "cuda_device_count": int(torch.cuda.device_count()),
            "deterministic_thread_limit": 1,
        }
    )


def self_test(output_dir: str) -> Dict[str, object]:
    """Run a one-epoch CPU smoke test on a deterministic small matrix."""

    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(123)
    n_genes, n_cells, n_modules = 18, 24, 3
    module_signal = rng.normal(size=(n_modules, n_cells)).astype(np.float32)
    expression = np.empty((n_genes, n_cells), dtype=np.float32)
    for gene in range(n_genes):
        module = gene // (n_genes // n_modules)
        expression[gene] = module_signal[module] + 0.2 * rng.normal(size=n_cells)
    expression -= expression.min()
    input_path = destination / "self_test_expression.tsv"
    pd.DataFrame(
        expression,
        index=[f"gene_{index + 1}" for index in range(n_genes)],
        columns=[f"cell_{index + 1}" for index in range(n_cells)],
    ).to_csv(input_path, sep="\t", index=True, index_label="gene")
    return train_genecluster(
        input_path=str(input_path),
        output_dir=str(destination / "output"),
        n_clusters=n_modules,
        epochs=1,
        batch_size=6,
        seed=123,
        device="cpu",
        num_workers=0,
        patience=-1,
        deterministic=True,
        overwrite=True,
        verbose=False,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the DeepCSCN GeneCluster backend")
    parser.add_argument("--input", help="Gene-by-cell tab-delimited expression matrix")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--n-clusters", type=int, default=15)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--momentum", type=float, default=0.9)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--focal-gamma", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser


def main() -> None:
    args = _parser().parse_args()
    if args.self_test:
        result = self_test(args.output)
    else:
        if not args.input:
            raise SystemExit("--input is required unless --self-test is used")
        result = train_genecluster(
            input_path=args.input,
            output_dir=args.output,
            n_clusters=args.n_clusters,
            epochs=args.epochs,
            batch_size=args.batch_size,
            learning_rate=args.learning_rate,
            momentum=args.momentum,
            weight_decay=args.weight_decay,
            focal_gamma=args.focal_gamma,
            seed=args.seed,
            device=args.device,
            num_workers=args.num_workers,
            patience=args.patience,
            overwrite=args.overwrite,
        )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
