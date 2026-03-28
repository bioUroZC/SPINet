#===============================================================
# 0) Clean & libraries
#===============================================================
rm(list = ls())

library(Matrix)
library(GSEABase)
library(sparseMatrixStats)
library(UCell)
library(AUCell)
library(GSVA)
library(data.table)
library(GSEABase)

#===============================================================
# 1) Paths: MTX + GMT
#===============================================================

cancer_type = "GSE100501"
out_dir  <-  paste0("/proj/c.zihao/work2/1pathway/2pertu/", cancer_type, "/pathway/1base/")
read_dir <- paste0("/proj/c.zihao/work2/1pathway/2pertu/", cancer_type, "/out/")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


gmt_file1 <- "/proj/c.zihao/work2/2pathway/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt"
gmt1 <- GSEABase::getGmt(gmt_file1)
gmt_file2 <- "/proj/c.zihao/work2/2pathway/h.all.v2025.1.Hs.symbols.gmt"
gmt2 <- GSEABase::getGmt(gmt_file2)
gmt_file3 <- "/proj/c.zihao/work2/2pathway/c2.cp.reactome.v2025.1.Hs.symbols.gmt"
gmt3 <- GSEABase::getGmt(gmt_file3)

gs1 <- as(gmt1, "list")
gs2 <- as(gmt2, "list")
gs3 <- as(gmt3, "list")

names(gs1) <- paste0("KEGG_", names(gs1))
names(gs2) <- paste0("HALLMARK_",  names(gs2))
names(gs3) <- paste0("REACTOME_",  names(gs3))

gmt <- GeneSetCollection(c(gs1, gs2, gs3))


class(gmt)
print(length(gmt))
head(names(gmt))
length(GSEABase::geneIds(gmt[[1]]))



sample_id = basename(list.dirs(read_dir, recursive = FALSE, full.names = TRUE))

print(sample_id)

for(i in 1:length(sample_id)) {
  
  data_dir <- file.path(read_dir, sample_id[i])
  sid <- basename(data_dir)
  cat("\n==================== ", sid, " ====================\n")
  
  setwd(data_dir)
  data   <- fread('logNorm.txt')
  data <- as.data.frame(data)
  data[1:5,1:5]
  names(data)[1] <- 'gene'
  rownames(data) <- data$gene
  data$gene <- NULL
  data[1:5,1:5]

  data  <- as.matrix(data)
  data   <- as(data, "CsparseMatrix")
  
  cat("data   dim:", dim(data),   " nnz:", length(data@x), "\n")
  
  #===============================================================
  # 3) Read GMT and filter gene sets by intersection with expression genes
  #===============================================================
  
  geneset_list <- setNames(lapply(gmt, GSEABase::geneIds), names(gmt))
  
  # unify symbols: '_' -> '-'  (keep your original rule)
  geneset_list <- lapply(geneset_list, function(x) gsub("_", "-", x))
  
  # intersect with expression genes and filter small sets
  geneset_list <- lapply(geneset_list, function(x) intersect(x, rownames(data)))
  geneset_list <- geneset_list[lengths(geneset_list) >= 10]
  
  cat("Gene sets kept (>=10 genes after intersect):", length(geneset_list), "\n")
  stopifnot(length(geneset_list) > 0)
  
  #===============================================================
  # 4) Method 1: z-score 
  #===============================================================

  expr_dense <- as.matrix(data)
  expr_dense[1:5,1:5]
  

  cell_mean <- colMeans(expr_dense, na.rm = TRUE)                 
  cell_sd   <- apply(expr_dense, 2, sd, na.rm = TRUE)             
  eps <- 1e-8
  cell_sd[cell_sd < eps] <- eps
  
  
  pathway_mean_list <- lapply(geneset_list, function(gs) {
    colMeans(expr_dense[gs, , drop = FALSE], na.rm = TRUE)
  })
  
  pathway_mean <- do.call(rbind, pathway_mean_list)         
  rownames(pathway_mean) <- names(geneset_list)
  pathway_z <- sweep(pathway_mean, 2, cell_mean, FUN = "-")
  combined_z_mat <- sweep(pathway_z,    2, cell_sd,   FUN = "/")
  combined_z_mat[1:5,1:5]
  cat("combined_z_mat:", dim(combined_z_mat), "\n")
  
  #===============================================================
  # 5) Method 2: UCell (works on sparse)
  #===============================================================
  ucell_scores <- ScoreSignatures_UCell(
    matrix   = data,
    features = geneset_list
  )
  # UCell output: cells x signatures (usually)
  cat("ucell_scores:", dim(ucell_scores), "\n")
  
  # If you want geneSet x cell, transpose:
  ucell_mat <- t(as.matrix(ucell_scores))
  cat("ucell_mat (geneSet x cell):", dim(ucell_mat), "\n")
  
  #===============================================================
  # 6) Method 3: AUCell (works on sparse)
  #===============================================================
  cells_rankings <- AUCell_buildRankings(data, plotStats = FALSE)
  auc <- AUCell_calcAUC(geneset_list, cells_rankings)
  aucell_mat <- as.matrix(getAUC(auc))   # geneSet x cell
  cat("aucell_mat:", dim(aucell_mat), "\n")
  
  #===============================================================
  # 7) Method 4: GSVA / ssGSEA (optional if GSVA installed)
  #===============================================================
  
  gsva_mat <- GSVA::gsva(
    expr_dense,
    geneset_list,
    method = "gsva",
    kcdf = "Gaussian",
    verbose = TRUE
  )
  
  ssgsea_mat <- GSVA::gsva(
    expr_dense,
    geneset_list,
    method = "ssgsea",
    kcdf = "Gaussian",
    verbose = TRUE
  )
  
  cat("gsva_mat:", dim(gsva_mat), "\n")
  cat("ssgsea_mat:", dim(ssgsea_mat), "\n")
  
  #===============================================================
  # 8) Method 5: scSE 
  #===============================================================
  run_scSE <- function(expr_counts, geneset_list) {
    stopifnot(all(colnames(expr_counts) != ""))
    
    umi_total <- sparseMatrixStats::colSums2(expr_counts, na.rm = TRUE)
    umi_total[umi_total == 0] <- NA
    
    # Precompute gene index for speed
    gene_index <- setNames(seq_len(nrow(expr_counts)), rownames(expr_counts))
    
    res <- lapply(geneset_list, function(gs) {
      idx <- gene_index[intersect(gs, names(gene_index))]
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) return(rep(NA, ncol(expr_counts)))
      
      umi_gs <- sparseMatrixStats::colSums2(expr_counts, rows = idx, na.rm = TRUE)
      (umi_gs / umi_total) * 100
    })
    
    mat <- do.call(rbind, res)
    rownames(mat) <- names(geneset_list)
    colnames(mat) <- colnames(expr_counts)
    mat
  }
  
  scSE_mat <- run_scSE(data, geneset_list)
  cat("scSE_mat:", dim(scSE_mat), "\n")
  
  #===============================================================
  # 9) Method 6: JASMINE (requires dense)
  #===============================================================
  RankCalculation <- function(x, genes) {
    subdata <- x[x != 0]
    if (length(subdata) == 0) return(0)
    ranks <- rank(subdata)
    ranks_sig <- ranks[names(ranks) %in% genes]
    if (length(ranks_sig) == 0) return(0)
    mean(ranks_sig, na.rm = TRUE) / length(subdata)
  }
  
  ORCalculation <- function(data_mtx, genes) {
    GE  <- data_mtx[rownames(data_mtx) %in% genes, , drop = FALSE]
    NGE <- data_mtx[!rownames(data_mtx) %in% genes, , drop = FALSE]
    
    SigExp  <- colSums(GE != 0)
    NSigExp <- colSums(NGE != 0)
    
    SigNE <- nrow(GE) - SigExp
    SigNE[SigNE == 0] <- 1
    NSigExp[NSigExp == 0] <- 1
    
    NSigNE <- nrow(data_mtx) - (SigExp + NSigExp)
    NSigNE <- NSigNE - SigNE
    
    (SigExp * NSigNE) / (SigNE * NSigExp)
  }
  
  normalize01 <- function(x) {
    if (all(x == 0)) return(x)
    (x - min(x)) / (max(x) - min(x))
  }
  
  JASMINE_score <- function(data_mtx, genes) {
    genes <- intersect(genes, rownames(data_mtx))
    if (length(genes) < 2) return(rep(NA, ncol(data_mtx)))
    
    RM <- apply(data_mtx, 2, RankCalculation, genes = genes)
    RM <- normalize01(RM)
    
    OR <- ORCalculation(data_mtx, genes)
    OR <- normalize01(OR)
    
    (RM + OR) / 2
  }
  
  run_jasmine <- function(expr_dense, geneset_list) {
    res <- lapply(geneset_list, function(gs) JASMINE_score(expr_dense, gs))
    mat <- do.call(rbind, res)
    rownames(mat) <- names(geneset_list)
    colnames(mat) <- colnames(expr_dense)
    mat
  }
  
  jasmine_mat <- run_jasmine(expr_dense, geneset_list)
  cat("jasmine_mat:", dim(jasmine_mat), "\n")
  
  #===============================================================
  # 10) Method 7: PLAGE (requires dense)
  #===============================================================
  run_plage <- function(expr_dense, geneset_list) {
    res <- lapply(geneset_list, function(gs) {
      gs <- intersect(gs, rownames(expr_dense))
      if (length(gs) < 2) return(rep(NA, ncol(expr_dense)))
      
      submat <- expr_dense[gs, , drop = FALSE]
      
      gene_sd <- apply(submat, 1, sd)
      submat <- submat[gene_sd > 0, , drop = FALSE]
      if (nrow(submat) < 2) return(rep(NA, ncol(expr_dense)))
      
      pca <- prcomp(t(submat), center = TRUE, scale. = TRUE)
      pca$x[, 1]
    })
    
    mat <- do.call(rbind, res)
    rownames(mat) <- names(geneset_list)
    colnames(mat) <- colnames(expr_dense)
    mat
  }
  
  plage_mat <- run_plage(expr_dense, geneset_list)
  cat("plage_mat:", dim(plage_mat), "\n")
  
  #===============================================================
  # 11) (Optional) Save results to disk
  #===============================================================
  
  file_to_save <- paste0(sid, "_pathway_scores.rds")
  
  saveRDS(
    list(
      zscore = combined_z_mat,
      ucell      = ucell_mat,
      aucell     = aucell_mat,
      gsva       = gsva_mat,
      ssgsea     = ssgsea_mat,
      scSE       = scSE_mat,
      jasmine    = jasmine_mat,
      plage      = plage_mat,
      genesets   = geneset_list
    ),
    file = file.path(out_dir, file_to_save)
  )
  
  cat("Done. Saved: pathway_scores.rds in ", out_dir, "\n")
  print('===============================================================')
  
}

print("All samples done.")
