import time
import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import csr_matrix


def prepare_data_RWR(links, expression_data):
    # Now get common genes with updated expression_data
    genes_in_links = pd.unique(links[['protein1', 'protein2']].values.ravel())
    common_genes = expression_data.index.intersection(genes_in_links)
    # Filter links based on updated gene list
    link_filtered = links[
        links['protein1'].isin(common_genes) & links['protein2'].isin(common_genes)
        ].drop_duplicates(subset=['protein1', 'protein2'])
    links = link_filtered
    print(links.shape)
    # Filter expression matrix again (just to be safe)
    used_genes = pd.unique(link_filtered[['protein1', 'protein2']].values.ravel())
    expression_data = expression_data.loc[expression_data.index.isin(used_genes)]
    # === Construct Protein Interaction Network ===

    G = nx.Graph()
    for _, row in links.iterrows():
        G.add_edge(row['protein1'], row['protein2'], weight=row['score'])

    # Store real original weights for each edge
    real_original_weights = {
        (row['protein1'], row['protein2']): row['score']
        for _, row in links.iterrows()
    }

    return G, real_original_weights, expression_data




def RWR_single_sample(sample_id, G, real_original_weights, expression_data, nodes, node_to_index):
    print(f"\nProcessing Sample: {sample_id}")

    # Step 1: Reset edge weights to real original values
    for (u, v), weight in real_original_weights.items():
        if G.has_edge(u, v):  # Check to avoid missing edges
            G[u][v]['weight'] = weight

    sample_expression = expression_data[sample_id]

    # === Step 2: Build Transition Matrix ===
    n = len(nodes)

    row, col, edge_weights = [], [], []
    for u, v in G.edges():
        weight = G[u][v]['weight']
        row.append(node_to_index[u])
        col.append(node_to_index[v])
        edge_weights.append(weight)

        row.append(node_to_index[v])
        col.append(node_to_index[u])
        edge_weights.append(weight)

    # Create and normalize the sparse transition matrix
    T_sparse = csr_matrix((edge_weights, (row, col)), shape=(n, n))
    row_sums = np.array(T_sparse.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    T_sparse = T_sparse.multiply(1 / row_sums[:, None])  # Normalize rows

    # === Step 3: Random Walk Calculation ===
    start_time = time.time()
    expr_vec = sample_expression.reindex(nodes).fillna(0.0).values.astype(float)

    s = expr_vec.sum()
    if s == 0:
        expr_vec[:] = 1.0
        s = expr_vec.sum()
    P0 = expr_vec / s

    rwr_alpha = 0.5

    if not np.isfinite(P0).all():
        P0[:] = 1.0
        P0 = P0 / P0.sum()

    P = P0.copy()
    max_iter = 50  # Limit the number of iterations
    tol = 1e-4  # Relaxed convergence tolerance

    for step in range(max_iter):
        # Sparse matrix-vector multiplication
        P_new = rwr_alpha * P0 + (1 - rwr_alpha) * T_sparse.T @ P
        # Dense norm for convergence check
        if np.linalg.norm(P - P_new) < tol:
            print(f"Converged at Step {step + 1}")
            P = P_new
            break
        P = P_new

    # Normalize P
    P_normalized = np.log1p(P)
    P_normalized = P_normalized / (P_normalized.sum() + 1e-12)

    step_time = time.time() - start_time
    print(f"Step 3 (Random Walk Calculation) Time: {step_time:.4f} seconds")

    # Step 4:  Save Current Sample
    node_scores = pd.DataFrame({
        "Node": nodes,
        "Score": P_normalized
    })
    node_scores["Sample"] = sample_id

    return node_scores

