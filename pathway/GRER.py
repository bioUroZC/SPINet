# GRER.py
# Graph-regularized Expression Refinement (dropout-only PPI injection)

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from typing import List, Tuple, Dict, Optional


def build_ppi_row_normalized_adjacency(
    gene_names: List[str],
    edges: List[Tuple[str, str]],
    undirected: bool = True,
    dtype=np.float32,
) -> sp.csr_matrix:
    """
    Build row-normalized adjacency A_hat = D^{-1} A (genes x genes).
    """
    gene_to_idx: Dict[str, int] = {g: i for i, g in enumerate(gene_names)}
    n = len(gene_names)

    rows, cols = [], []
    for u, v in edges:
        if u not in gene_to_idx or v not in gene_to_idx:
            continue
        i, j = gene_to_idx[u], gene_to_idx[v]
        rows.append(i); cols.append(j)
        if undirected:
            rows.append(j); cols.append(i)

    data = np.ones(len(rows), dtype=dtype)
    A = sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=dtype)

    row_sums = np.asarray(A.sum(axis=1)).ravel()
    inv = np.zeros_like(row_sums, dtype=dtype)
    nz = row_sums > 0
    inv[nz] = 1.0 / row_sums[nz]

    A_hat = sp.diags(inv, format="csr") @ A
    return A_hat


def dropout_only_ppi_residual_injection_from_sparse(
    X_csr: sp.csr_matrix,      # rows x genes (rows = cells here)
    A_hat: sp.csr_matrix,      # genes x genes
    beta: float = 0.1,
    block_size: int = 256,
    dtype=np.float32,
) -> np.ndarray:
    """
    Dropout-only injection:
      if X == 0: X_new = beta * (X @ A_hat.T)
      else:      X_new = X
    Returns dense np.ndarray.
    """
    X_csr = X_csr.tocsr().astype(dtype)
    n_rows, n_genes = X_csr.shape
    if A_hat.shape != (n_genes, n_genes):
        raise ValueError(f"A_hat shape {A_hat.shape} != ({n_genes},{n_genes})")

    A_t = A_hat.T.tocsr().astype(dtype)
    X_corr = np.zeros((n_rows, n_genes), dtype=dtype)

    for start in range(0, n_rows, block_size):
        end = min(start + block_size, n_rows)
        Xb = X_csr[start:end]

        neigh = (Xb @ A_t).toarray()
        Xb_dense = Xb.toarray()

        mask = (Xb_dense == 0)
        Xb_dense[mask] = beta * neigh[mask]
        X_corr[start:end] = Xb_dense

    return X_corr


def grer_expression_df(
    expression_data: pd.DataFrame,   # genes x cells
    links_df: pd.DataFrame,          # PPI dataframe
    gene_col1: str = "protein1",
    gene_col2: str = "protein2",
    beta: float = 0.1,
    block_size: int = 256,
    dtype=np.float32,
    undirected: bool = True,
    report: bool = False,
) -> pd.DataFrame:
    """
    Apply GRER to expression_data (genes x cells) using PPI from links_df.

    - expression_data: index=genes, columns=cells
    - links_df: must contain two gene columns (default protein1/protein2)
    - Only fills zeros (dropout-only). Non-zeros unchanged.
    - Returns refined expression_data with same shape (genes x cells).
    """
    # Basic checks
    if gene_col1 not in links_df.columns or gene_col2 not in links_df.columns:
        raise ValueError(f"links_df missing required columns: {gene_col1}, {gene_col2}")

    if expression_data.index.has_duplicates:
        raise ValueError("expression_data index has duplicated gene names.")
    if expression_data.columns.has_duplicates:
        raise ValueError("expression_data columns has duplicated cell names.")

    genes = expression_data.index.astype(str).tolist()
    cells = expression_data.columns.astype(str).tolist()

    # Build edge list (only keep valid string pairs; NaN dropped)
    edges = list(
        links_df[[gene_col1, gene_col2]]
        .dropna()
        .astype(str)
        .itertuples(index=False, name=None)
    )

    # Build A_hat on the genes present in expression_data
    A_hat = build_ppi_row_normalized_adjacency(
        gene_names=genes,
        edges=edges,
        undirected=undirected,
        dtype=dtype,
    )

    # Convert to cells x genes for computation
    X_cg = expression_data.T.astype(dtype)
    X_csr = sp.csr_matrix(X_cg.values)

    if report:
        zr_before = 1.0 - X_csr.nnz / (X_csr.shape[0] * X_csr.shape[1])
        print(f"[GRER] Zero rate BEFORE (global): {zr_before:.4f}")
        print(f"[GRER] A_hat nnz: {A_hat.nnz}")

    X_corr = dropout_only_ppi_residual_injection_from_sparse(
        X_csr=X_csr,
        A_hat=A_hat,
        beta=beta,
        block_size=block_size,
        dtype=dtype,
    )

    if report:
        zr_after = float(np.mean(X_corr == 0))
        print(f"[GRER] Zero rate AFTER (global): {zr_after:.4f}")
        print(f"[GRER] Absolute reduction: {zr_before - zr_after:.4f}")

    refined = pd.DataFrame(X_corr, index=cells, columns=genes).T
    return refined
