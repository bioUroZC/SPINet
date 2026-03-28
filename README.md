---

# 🧬 SPINet

**Network-informed diffusion framework for cell-specific pathway activity and protein interaction inference from single-cell RNA-seq**

---



<img width="2047" height="1091" alt="eb2d9a91-12b8-461e-ae24-60b0a014a323" src="https://github.com/user-attachments/assets/56ad8759-b8c5-4f99-b93a-6c6cfa0e8e7e" />


## 📌 Overview

**SPINet** is a network-based computational framework designed to infer:

* **Cell-specific pathway activities**
* **Cell-specific protein-protein interaction (PPI) strengths**

from **single-cell RNA-seq (scRNA-seq)** data.

It addresses a key challenge in single-cell analysis:

> ⚠️ *Dropout-induced sparsity obscures weak but coordinated biological signals.*

SPINet integrates:

* **Network diffusion (Random Walk with Restart, RWR)**
* **Graph Neural Networks (GNNs)**
* **Self-supervised representation learning**

to leverage prior biological knowledge (PPI networks) and recover functional signals at single-cell resolution.

---

## 🚀 Key Features

* 🔬 **Dropout-robust pathway inference**
* 🧠 **Diffusion-informed node activity scoring**
* 🔗 **Cell-specific PPI network inference**
* 🤖 **Graph Attention Network (GATv2) encoder**
* 🔁 **Self-supervised graph autoencoder training**
* 📊 **Improved AUC-ROC / AUPRC performance**
* 🧪 **Validated across perturbation, drug response, and multimodal datasets**
* ⚡ **Scalable to large single-cell datasets**

---

## 🧩 Method Overview

SPINet consists of two major components:

### 1. Diffusion-based Node Activity Inference

* Construct a **weighted PPI network** (e.g., STRING)
* Perform **Random Walk with Restart (RWR)** per cell
* Generate **diffusion-informed gene activity scores**

### 2. Cell-specific Interaction Inference via GNN

* Construct node features:

  * Gene expression
  * Diffusion scores
  * Network degree
* Train a **Graph Attention Network (GATv2) encoder**
* Use **self-supervised reconstruction (autoencoder)**
* Infer **edge-level interaction strength**:

  ```
  Score(i, j) = sigmoid( Z_i^T Z_j )
  ```

---

## 🏗️ Framework Architecture

```
scRNA-seq → PPI Network → RWR Diffusion → Node Features
                                      ↓
                              GNN Encoder (GATv2)
                                      ↓
                               Latent Embedding
                                      ↓
                          Decoder (Feature Reconstruction)
                                      ↓
                      Cell-specific PPI Interaction Scores
                                      ↓
                         Pathway Activity Inference
```

---

## 📂 Repository Structure

```
SPINet/
│
├── data/                  # Input datasets (scRNA-seq, gene sets, PPI)
├── preprocessing/        # QC, normalization, integration
├── network/              # PPI construction and processing
├── diffusion/            # RWR implementation
├── model/
│   ├── encoder.py        # GATv2 encoder
│   ├── decoder.py        # MLP decoder
│   ├── training.py       # self-supervised training
│
├── inference/            # Edge scoring and pathway activity
├── benchmarking/         # Evaluation scripts
├── utils/                # Helper functions
│
├── notebooks/            # Example workflows
├── results/              # Output results
└── README.md
```

---

## ⚙️ Installation

```bash
git clone https://github.com/your-username/SPINet.git
cd SPINet

conda create -n spinet python=3.10
conda activate spinet

pip install -r requirements.txt
```

---

## 📥 Input Requirements

### Required inputs:

* **scRNA-seq expression matrix**

  * format: cells × genes
* **PPI network**

  * e.g., STRING high-confidence interactions
* **Gene sets**

  * KEGG / Reactome / Hallmark

---

## ▶️ Usage

### Step 1: Preprocessing

```bash
python preprocessing/run_qc.py
```

### Step 2: Diffusion (RWR)

```bash
python diffusion/run_rwr.py
```

### Step 3: Train GNN

```bash
python model/training.py
```

### Step 4: Interaction Inference

```bash
python inference/run_edge_scoring.py
```

### Step 5: Pathway Activity

```bash
python inference/run_pathway_scoring.py
```

---

## 📊 Output

* 📈 **Pathway activity matrix** (pathway × cell)
* 🔗 **Cell-specific PPI interaction matrix** (edge × cell)
* 🧬 **Latent gene embeddings**
* 📉 **Benchmarking metrics (AUC, AUPRC)**

---

## 🧪 Benchmarking

SPINet has been evaluated against:

* AUCell
* UCell
* JASMINE
* ssGSEA
* GSVA
* PLAGE
* Z-score baseline

### Results:

* ✅ Higher sensitivity to weak pathways
* ✅ Robust under increasing dropout
* ✅ Better cell-type discrimination
* ✅ Improved network inference accuracy

---

## 🔬 Applications

* Cell-type-specific pathway analysis
* Drug response and resistance modeling
* Immune cell functional profiling
* Network rewiring under perturbations
* Multimodal integration (e.g., CITE-seq)

---

## 📦 Dependencies

* Python ≥ 3.9
* PyTorch
* PyTorch Geometric
* Scanpy / Seurat-compatible preprocessing
* NetworkX
* NumPy / SciPy / Pandas

---

## 📖 Citation

If you use SPINet, please cite:

```
[Your paper citation here]
```

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

---

## 📬 Contact

For questions or collaborations:

* Email: [your_email@domain.com](mailto:your_email@domain.com)
* GitHub Issues

---

## ⭐ Acknowledgements

* STRING database
* MSigDB
* Single-cell community tools

---

## 📄 License

This project is licensed under the MIT License.

---
