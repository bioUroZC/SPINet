import os
import time
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys

sys.path.append(r"/proj/c.zihao/work2/function/")
import GRER
import RWRnode

# ============================================
# Config
# ============================================

cancer_type = "GSE100501"
out_dir = f"/proj/c.zihao/work2/1pathway/2pertu/{cancer_type}/pathway/3GRER/"
read_dir = f"/proj/c.zihao/work2/1pathway/2pertu/{cancer_type}/out/"
links_path = "/proj/c.zihao/work1/2survival/links.csv"

# Ensure output directory exists
os.makedirs(out_dir, exist_ok=True)

# ============================================
# Discover datasets
# ============================================

subdirs = [
    d for d in os.listdir(read_dir)
    if os.path.isdir(os.path.join(read_dir, d))
]
print("Datasets:", subdirs)
dataset_ids = subdirs


links = pd.read_csv(links_path, index_col=0)
print(links.head())
print("links.shape:", links.shape)
print("Links Data Loaded.")



# ============================================
# Process each dataset
# ============================================

for dataset_name in dataset_ids:
    print("#==========================================================")
    print("Dataset:", dataset_name)

    base_dir = os.path.join(read_dir, dataset_name)
    txt_file = os.path.join(base_dir, "logNorm.txt")

    # ---------- Load expression ----------
    expression_data = pd.read_csv(txt_file, sep="\t", index_col=0)
    print(expression_data.iloc[:3, :3])

    #expression_data = expression_data.iloc[:, :20]
    print("expression_data.shape:", expression_data.shape)
    print("Successfully read expression_data")
    
    refined_data = GRER.grer_expression_df(
        expression_data=expression_data,
        links_df=links,
        gene_col1="protein1",
        gene_col2="protein2",
        report=True
    )

    expression_data = refined_data

    # ---------- Prepare network ----------
    t0 = time.time()
    G, real_original_weights, expression_data = RWRnode.prepare_data_RWR(
        links=links,
        expression_data=expression_data
    )
    print(f"prepare_data_RWR time: {time.time() - t0:.4f} seconds")
    print(f"Network created. nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")

    # ---------- Fixed node order ----------
    nodes = sorted(G.nodes())
    node_to_index = {node: i for i, node in enumerate(nodes)}
    samples = list(expression_data.columns)

    # ---------- Pre-allocate node matrix ----------
    node_mat = np.empty((len(nodes), len(samples)), dtype=np.float32)

    # ---------- Per-sample run ----------
    codestart_time = time.time()

    for j, sample_id in enumerate(tqdm(samples, desc=f"Processing {dataset_name}", unit="sample")):

        node_scores = RWRnode.RWR_single_sample(
            sample_id=sample_id,
            G=G,
            real_original_weights=real_original_weights,
            expression_data=expression_data,
            nodes=nodes,
            node_to_index=node_to_index
        )

        # ===== Node: write column j =====
        node_mat[:, j] = node_scores["Score"].to_numpy(dtype=np.float32)

    elapsed_time = time.time() - codestart_time
    print(f"[{dataset_name}] Execution time: {elapsed_time:.4f} seconds")

    # ============================================
    # Build wide tables
    # ============================================

    # ---------- Save node wide ----------
    node_scores_wide = pd.DataFrame(node_mat, index=nodes, columns=samples)
    node_scores_wide.reset_index().rename(columns={"index": "Node"}).to_csv(
        os.path.join(out_dir, f"{dataset_name}_node_score_wide.txt"),
        sep="\t",
        index=False
    )

    # ---------- Runtime log ----------
    log_out = os.path.join(out_dir, f"{dataset_name}_runtime_log.txt")
    with open(log_out, "w", encoding="utf-8") as f:
        f.write(f"Execution time: {elapsed_time:.4f} seconds\n")
    print("Saved runtime log to:", log_out)
