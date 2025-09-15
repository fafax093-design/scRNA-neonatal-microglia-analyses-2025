
##############################################

### Description: Single-cell RNA-seq analysis pipeline
### Author: 
### Date: 
##############################################

### Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggsci)
library(gplots)

###-----------------------------
### 1. Load raw data and create Seurat objects
###-----------------------------
# 读取所有样本数据，并命名清晰
Lps_5_1 <- readRDS("5L_1_seurat.rds")
Lps_5_2 <- readRDS("5L_2_seurat.rds")
NS_1 <- readRDS("NS_1_seurat.rds")
NS_2 <- readRDS("NS_2_seurat.rds")
LPS_5_P7_1 <- Read10X("J:\\zjj\\7d\\data\\matrix\\2308154_LPS_1_P7\\filtered_feature_bc_matrix")
LPS_5_P7_2 <- Read10X("J:\\zjj\\7d\\data\\matrix\\2308154_LPS_2_P7\\filtered_feature_bc_matrix")
NS_P7_1 <- Read10X("J:\\zjj\\7d\\data\\matrix\\2308154_NS_1_P7\\filtered_feature_bc_matrix")
NS_P7_2<- Read10X("J:\\zjj\\7d\\data\\matrix\\2308154_NS_2_P7\\filtered_feature_bc_matrix")
LPS_5_LPS_P12_1 <- Read10X("J:/zjj/12d/data/matrix/2308153_5LPS_1_P12/filtered_feature_bc_matrix")
LPS_5_LPS_P12_2 <- Read10X("J:/zjj/12d/data/matrix/2308153_5LPS_2_P12/filtered_feature_bc_matrix")
LPS_10_LPS_24h <- Read10X("J:/zjj/12d/data/matrix/2308153_10LPS_24h/filtered_feature_bc_matrix")
NS_P12_1 <- Read10X("J:/zjj/12d/data/matrix/2308153_NS_1_P12/filtered_feature_bc_matrix")
NS_P12_2 <- Read10X("J:/zjj/12d/data/matrix/2308153_NS_2_P12/filtered_feature_bc_matrix")

# 创建Seurat对象
LPS_5_P7_1 <- CreateSeuratObject(counts = LPS_5_P7_1, project = "LPS_5_P7_1", min.cells = 3, min.features = 200)
LPS_5_P7_2 <- CreateSeuratObject(counts = LPS_5_P7_2, project = "LPS_5_P7_2", min.cells = 3, min.features = 200)
NS_P7_1 <- CreateSeuratObject(counts = NS_P7_1, project = "NS_P7_1", min.cells = 3, min.features = 200)
NS_P7_2 <- CreateSeuratObject(counts = NS_P7_2, project = "NS_P7_2", min.cells = 3, min.features = 200)
LPS_5_LPS_P12_1 <- CreateSeuratObject(counts = LPS_5_LPS_P12_1, project = "LPS_5_LPS_P12_1", min.cells = 3, min.features = 200)
LPS_5_LPS_P12_2<- CreateSeuratObject(counts = LPS_5_LPS_P12_2, project = "LPS_5_LPS_P12_2", min.cells = 3, min.features = 200)
LPS_10_LPS_24h <- CreateSeuratObject(counts = LPS_10_LPS_24h, project = "LPS_10_LPS_24h", min.cells = 3, min.features = 200)
NS_P12_1 <- CreateSeuratObject(counts = NS_P12_1, project = "NS_P12_1", min.cells = 3, min.features = 200)
NS_P12_2 <- CreateSeuratObject(counts = NS_P12_2, project = "NS_P12_2", min.cells = 3, min.features = 200)
# ... repeat for others

###-----------------------------
### 2. Merge and Quality Control (QC)
###-----------------------------
# 合并所有Seurat对象并创建用于下游分析的新Seurat对象
metadata.sub <- c("orig.ident","nCount_RNA","nFeature_RNA")
scRNA <- CreateSeuratObject(merged_seurat@assays$RNA@counts, meta.data = merged_seurat@meta.data[,metadata.sub])
scRNA1 <- CreateSeuratObject(scRNA1@assays$RNA@counts, meta.data = scRNA1@meta.data[,metadata.sub])


# 计算质控指标
# 计算细胞中线粒体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
scRNA$percent.mt1 <- round(scRNA$percent.mt, 3)
# 计算细胞中核糖体基因比例
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^Rp[ls]")
# 计算红细胞比例
HB.genes <- c("Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbe1", "Hbg1", "Hbg2", "Hbm", "Hbq1", "Hbz")
HB.genes <- CaseMatch(HB.genes, rownames(scRNA))
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 


# 保存 QC summary（建议作为Supplementary Table）
qc_summary <- scRNA@meta.data %>% select(orig.ident, nFeature_RNA, nCount_RNA, percent.mt, percent.rb, percent.HB)
write.csv(summary(qc_summary), "QC/QC_metric_summary.csv")

# 可视化过滤前指标
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
qc_plots <- lapply(qc_features, function(f) {
  VlnPlot(scRNA, group.by = "orig.ident", features = f, pt.size = 0, raster = FALSE) + 
    NoLegend() + ggtitle(f)
})
ggsave("QC/vlnplot_qc.pdf", plot = wrap_plots(qc_plots, nrow = 3), width = 16, height = 8)

# 阈值可视化（示意）
VlnPlot(scRNA, features = "nFeature_RNA") + geom_hline(yintercept = c(500, 3000))
VlnPlot(scRNA, features = "nCount_RNA") + geom_hline(yintercept = 22000)
VlnPlot(scRNA, features = "percent.mt") + geom_hline(yintercept = 10)

# 过滤低质量细胞
scRNA <- subset(scRNA, subset = 
                  nFeature_RNA > 500 & nFeature_RNA < 3000 &
                  nCount_RNA < 22000 &
                  percent.mt < 10)

# 可视化过滤后结果
qc_plots_post <- lapply(qc_features, function(f) {
  VlnPlot(scRNA, group.by = "orig.ident", features = f, pt.size = 0.01, raster = FALSE) + 
    NoLegend() + ggtitle(paste0(f, " (Post-QC)"))
})
ggsave("QC/vlnplot_after_qc.pdf", plot = wrap_plots(qc_plots_post, nrow = 3), width = 16, height = 8)

# 补充散点图（nFeature vs nCount）
pdf("QC/feature_scatter.pdf", width = 6, height = 5)
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


#####batch
scRNA@meta.data <- scRNA@meta.data %>%
  mutate(batch = case_when(
    orig.ident %in% c("NS", "Lps") ~ "batch1",
    orig.ident %in% c("Lps_5_1", "Lps_5_2","NS_1","NS_2") ~ "batch2",
    orig.ident %in% c("LPS_10_LPS_24h","NS_P12_1","NS_P12_2","LPS_5_LPS_P12_1","LPS_5_LPS_P12_2") ~ "batch3",
    orig.ident %in% c("NS_P7_1","NS_P7_2","LPS_5_P7_1","LPS_5_P7_2") ~ "batch4",
    TRUE ~ orig.ident))


#---------  3. 降维聚类分群         -------------------------------
### Load required packages
library(Seurat)      # v4.2.1
library(dplyr)       # v1.1.1
library(harmony)     # v0.1.1

#  ---------- 3.1: With Batch Integration (Harmony)---------------------------

# Define metadata to keep
metadata.sub <- c("orig.ident", "S.Score", "G2M.Score", "Phase", "group", "batch")

# Create Seurat object
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts,
                            meta.data = scRNA@meta.data[, metadata.sub])

# SCTransform normalization
scRNA <- SCTransform(scRNA)

# Run PCA
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)

# Run Harmony integration based on batch variable
scRNA <- RunHarmony(scRNA,
                    group.by.vars = "batch",
                    assay.use = "SCT",
                    max.iter.harmony = 20)

# Visualize elbow plot to determine PCs
ElbowPlot(scRNA, ndims = 50)

# Select PCs for downstream analysis
pc.num <- 1:30

# UMAP based on Harmony-reduced dimensions
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = pc.num)

# Clustering using different resolutions
for (i in c(0.1, 0.3, 0.5, 0.7, 1)) {
  scRNA <- scRNA %>%
    FindNeighbors(reduction = "harmony", dims = pc.num) %>%
    FindClusters(resolution = i)}

# Save object
saveRDS(scRNA, file = "scRNA_HarmonyIntegrated.rds")


###-------------------------------------------------------------------------
#-------3.2: Without Batch Integration -----------------------------------
metadata.sub <- c("orig.ident", "S.Score", "G2M.Score", "Phase", "group", "batch")
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts,
                            meta.data = scRNA@meta.data[, metadata.sub])

# SCTransform normalization
scRNA <- SCTransform(scRNA)

# Run PCA
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)

# Elbow plot to inspect PC variance
ElbowPlot(scRNA, ndims = 50)

# Select PCs
pc.num <- 1:30

# Run UMAP on PCA space (no batch correction)
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = pc.num)

# Clustering using a range of resolutions
for (i in c(0.1, 0.2, 0.3, 0.5, 0.7, 1)) {
  scRNA <- scRNA %>%
    FindNeighbors(reduction = "pca", dims = pc.num) %>%
    FindClusters(resolution = i)}

# Save object
saveRDS(scRNA, file = "scRNA_NoBatchCorrection.rds")

-------------------------------------------------------------------------------

### -------4--Marker Gene Identification and DEG Analysis --------------


### Load required packages
library(Seurat)
library(future)        # v1.28.0

### Set parallel processing parameters
options(future.globals.maxSize = 2000000 * 1024^3)  # Set memory limit
plan(multisession, workers = 10)

###---------------------------------------------------------------------------
###------- 4.1: Calculate Marker Genes for All Clusters---------------------

Idents(scRNA) <- "celltype"
# Run marker gene identification
ClusterMarker <- FindAllMarkers(
  object = scRNA,
  assay = "SCT",
  slot = "data")
# Save marker gene results
write.csv(ClusterMarker, "Marker_Genes_for_P3.csv", row.names = FALSE)



###------------- 4.2: DEG Between Specific Clusters-----------------------------


options(future.globals.maxSize = 20000 * 1024^3)
plan(multisession, workers = 20)
# Compare between cluster "0" and "6" under clustering result SCT_snn_res.0.5
deg <- FindMarkers(
  object = scRNA,
  ident.1 = "0",
  ident.2 = "6",
  group.by = "SCT_snn_res.0.5",
  logfc.threshold = 0,
  min.pct = 0,
  assay = "RNA",
  slot = "counts")

# Save DEG result
write.csv(deg, file = "Cluster_0_vs_6_DEG.csv", row.names = TRUE)


--------------------------------------------------------------------------------

###亚群命名#   5. cluster annotation     ###########

metadata <- scRNA@meta.data
metadata$celltype <- recode(metadata$SCT_snn_res.0.5,
                            `0` = "NDM",`1` = "Mg1",`2` = "Mg2")
scRNA@meta.data <- metadata
scRNA$celltype<- factor(scRNA$celltype,levels = 
                          c("NDM","Mg1","Mg2"))

--------------------------------------------------------------------------------
###  ------6.Visualization ----------------------------------------------------

###-----6.1. Cluster Tree Visualization-----------------------------

library(clustree)   # Version 0.5.0

# Plot clustering relationships based on different resolutions
clustree(scRNA@meta.data, prefix = "SCT_snn_res.") +
  guides(edge_colour = FALSE, edge_alpha = FALSE) +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "blue", high = "red") +
  theme(legend.position = "bottom")
ggsave("ClusterTree.tiff", width = 8, height = 8, dpi = 600)

# ---------------6.2. UMAP by Cell Type and Group-----------------------------

library(SCP)   # Version 0.4.2
CellDimPlot(
  srt = scRNA,
  group.by = "celltype",
  reduction = "umap",
  theme_use = "theme_blank",
  label = TRUE,
  label_insitu = TRUE)
ggsave("umap_celltype.tiff", width = 6, height = 4, dpi = 600)
ggsave("umap_celltype.pdf", width = 4, height = 4, dpi = 600)

# UMAP split by group (e.g., NS_P7 vs LPS_P7)
CellDimPlot(
  srt = scRNA,
  group.by = "celltype",
  reduction = "umap",
  theme_use = "theme_blank",
  label = FALSE,
  label_insitu = FALSE,
  show_stat = FALSE,
  split.by = "group")
ggsave("umap_group.tiff", width = 15, height = 12, dpi = 600)

#---------------6.3. Feature Plot on UMAP  -----------------------------

# Plot gene expression (e.g., Gpnmb) on UMAP
FeatureDimPlot(
  object = scRNA,
  features = "Gpnmb",
  reduction = "UMAP",
  cells.highlight = TRUE,
  theme_use = "theme_blank",
  show_stat = FALSE,
  legend.position = "none")
ggsave("umap/Gpnmb.tiff", width = 2, height = 2, dpi = 600)


###----------------6.4. Cell Type Ratio per Group-----------------------------

# Cell proportion trend plot across groups
CellStatPlot(
  srt = scRNA,
  stat.by = "celltype",
  group.by = "group",
  plot_type = "trend")
ggsave("cellratio.pdf", width = 4, height = 3)
ggsave("cellratio.tiff", width = 4, height = 3, dpi = 600)


###----------------6.5. gene boxplot-----------------------------

FeatureStatPlot(srt = scRNA, stat.by = c("Pgk1","Pgam1","Pkm", "Ldha" ,"Aif1","P2ry12"),
                fill.by = "group", plot_type = "box", 
                group.by = "celltype",  bg.by = "celltype",stack = TRUE, flip = F) 
ggsave("gene.tiff", width= 8, height = 4,dpi = 600)



# ---------------------6.6 Gene Expression Heatmap ----------------------------

library(readxl)
library(Seurat)
library(ComplexHeatmap)##2.14.0
library(circlize)

### Load gene-cluster annotation
###-------------------------
gene_data <- read_excel("DATA/p3 all  cell.xlsx", sheet = 1)
colnames(gene_data) <- c("cluster", "gene")

### Extract expression matrix

exp_matrix <- GetAssayData(scRNA, assay = "SCT", slot = "data")

# 仅保留在表达矩阵中存在的基因
valid_genes <- gene_data$gene[gene_data$gene %in% rownames(exp_matrix)]

# 构建表达矩阵
h_data <- as.matrix(exp_matrix[valid_genes, ])

# 获取列注释（细胞类型）
celltypes <- scRNA$celltype[colnames(h_data)]

### 色板与注释设置
# 自定义颜色（示例色板）
color_palette <- c(
  '#E5D2DD', '#53A85F', '#F3B1A0', "#FFDD44", "#9467BD",
  "#E377C2", "#D62728", "#8C564B", "#2CA02C", '#23452F')
names(color_palette) <- unique(celltypes)

# 列注释（细胞类型）
top_anno <- columnAnnotation(
  celltype = factor(celltypes, levels = unique(celltypes)),
  col = list(celltype = color_palette),
  show_annotation_name = FALSE)

# 行注释（基因所属 cluster）
row_clusters <- gene_data$cluster[gene_data$gene %in% rownames(h_data)]
row_anno <- rowAnnotation(
  cluster = factor(row_clusters, levels = unique(row_clusters)),
  col = list(cluster = color_palette),
  show_annotation_name = FALSE)

###-------------------------
### 数据归一化函数（min-max scaling）

min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
scaled_data <- t(apply(h_data, 1, min_max_scale))

###-------------------------
### 绘制热图

ht <- Heatmap(
  scaled_data,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = factor(celltypes, levels = unique(celltypes)),
  row_split = factor(row_clusters, levels = unique(row_clusters)),
  top_annotation = top_anno,
  left_annotation = row_anno,
  col = colorRampPalette(c("#4978b3", "white", "#FF3333"))(100),
  heatmap_legend_param = list(
    at = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2),
    title = "Expression",
    title_position = "leftcenter-rot"
  ),
  border = TRUE,
  use_raster = TRUE,
  column_gap = unit(1, "mm"),
  row_gap = unit(1, "mm")
)

###-------------------------
### 保存为 PDF
pdf("FIGURE/GENE_EXPRESSION.pdf", width = 6.5, height = 6)
draw(ht)
dev.off()



### ----------6.7 DotPlot: Marker Gene Expression across Cell Types-------------

library(Seurat)
library(ggplot2)
library(readxl)

### Step 1. Load gene list
markers <- read_excel("gene.xlsx", sheet = 1)
markers <- as.character(markers$gene)

### Step 2. Set identity

Idents(scRNA) <- "celltype"

### Step 3. Generate DotPlot & extract data-----
p <- DotPlot(scRNA, features = markers)
dot_data <- p$data  # Extracted for custom ggplot

### Step 4. Plot using ggplot2---
p <- ggplot(dot_data, aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = NA))) + 
  scale_fill_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1', '#bef0b0', '#fdf4af',
                                   '#f9b64b', '#ec840e', '#ca443d', '#a51a49')) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank(),
    axis.text.y = element_text(color = 'black', size = 12),
    axis.text.x = element_text(color = 'black', size = 12, angle = 90, hjust = 1, vjust = 0.5))

### Step 5. Save to PDF
ggsave("FIGURE/DOTPLOT_NOT_SLICED.pdf", plot = p, width = 6, height = 2)


------------------------------------------------------------------------------
##7 --------------- GO Enrichment Analysis & Bubble Plot-------------


### Load Required Libraries
library(clusterProfiler)  ###4.2.2
library(dplyr)
library(org.Mm.eg.db)###3.14.0
library(GOplot)##1.0.2
library(ggplot2)
library(cowplot)
library(aplot)  # for insert_left
library(openxlsx)  # for reading Excel files

### Create output directory
dir.create("go_enrichment", showWarnings = FALSE)

### Load DEG list from Excel
DEG <- read.xlsx("module_gene.xlsx")

### Filter significant DEGs (upregulated & adjusted p < 0.01)
markers <- DEG %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0, p_val_adj < 0.01) %>%
  ungroup()

write.csv(markers, "go_enrichment/GO_gene.csv", row.names = FALSE)

### Map gene symbols to ENTREZ IDs
gid <- bitr(unique(markers$gene), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
markers <- full_join(markers, gid, by = c('gene' = 'SYMBOL'))

### GO enrichment analysis (Biological Process) by cluster
ego <- compareCluster(
  ENTREZID ~ cluster,
  data = markers,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP"
)

### Convert to readable format
ego_readable <- setReadable(ego, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
ego_df <- ego_readable@compareClusterResult
write.csv(ego_df, file = "go_enrichment/Enrichment_GO.csv", row.names = FALSE)

### Prepare data for GO bubble plot
path <- ego_df

### Add p-value group annotations
path$adjust_group <- cut(
  path$p.adjust,
  breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  labels = c("<0.0001", "<0.001", "<0.01", "<0.05", "<0.1", ">0.1")
)

### Ensure factor levels are ordered
path$adjust_group <- factor(path$adjust_group, levels = c("<0.0001", "<0.001", "<0.01", "<0.05", "<0.1", ">0.1"))
path$cluster <- factor(path$cluster, levels = c("Mg1", "Mg2", "Mg3", "Mg4", "Mg5", "Mg6", "Inf_1", "2a", "2b", "2c", "2d", "2e"))

### Main bubble plot
bubble_plot <- ggplot(path, aes(cluster, Description)) +
  geom_vline(aes(xintercept = cluster), color = "#D3D3D3") +
  geom_hline(aes(yintercept = Description), color = "#E8E8E8") +
  geom_point(aes(color = adjust_group, size = Count), shape = 19, stroke = 2) +
  scale_color_manual(values = c(
    "<0.0001" = "#67000D",
    "<0.001"  = "#EF3B2C",
    "<0.01"   = "#FB6A4A",
    "<0.05"   = "#FC9272",
    "<0.1"    = "#FEE0D2",
    ">0.1"    = "#FFF5F0"
  )) +
  cowplot::theme_cowplot() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5),
    legend.direction = "horizontal",
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL) +
  guides(
    color = guide_legend(override.aes = list(size = 4), order = 3),
    size  = guide_legend(override.aes = list(size = c(6, 7, 8)), order = 1)
  ) +
  scale_y_discrete(position = "right")

### Left annotation panel (optional)
# 'p' should be a dataframe with columns: Description, Annotation
left_anno <- p %>%
  as.data.frame() %>%
  mutate(p = "") %>%
  ggplot(aes(p, Description, fill = Annotation)) +
  geom_tile() +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  xlab(NULL) + ylab(NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
    legend.direction = "horizontal",
    legend.position = "top"
  ) +
  labs(fill = "Annotations") +
  scale_fill_manual(values = c("#8B4513", "#6B8E23", "#cc6d33", "#6c8c8b", "#4682B4"))

### Combine bubble plot and left annotation (if available)
combined_plot <- bubble_plot %>% insert_left(left_anno, width = 0.08)

### Save final plot
ggsave("FIGURE/GO_Bubble.pdf", plot = combined_plot, width = 18, height = 8)


------------------------------------------------------------------------------
### 8--UCell Pathway Score & Visualization -------------------------------------

library(UCell)  ## 2.1.1
library(ggrastr)## 1.01

# 读取基因集
markers <- readRDS("DATA/step0-Microglia_15_Pathways_geneset0.rds")

# 只保留 dendrite development 通路
features <- list("dendrite development" = markers[["dendrite development"]])

# UCell 打分
marker_score <- AddModuleScore_UCell(scRNA, features = features)

# 通路打分列名自动会是 "dendrite.development_UCell"，我们统一下名字
colnames(marker_score@meta.data)[grep("UCell$", colnames(marker_score@meta.data))] <- 
  gsub("\\.", " ", colnames(marker_score@meta.data)[grep("UCell$", colnames(marker_score@meta.data))])

# 提取分数数据
i <- "dendrite development_UCell"
data <- FetchData(marker_score, vars = c("group", i))
data$cellid <- data$group

# 提取细胞数
cell_number <- as.data.frame(table(data$cellid))
cell_number$number <- paste0("n=", cell_number$Freq)
colnames(cell_number)[1] <- "cellid"

# 计算每组最大值 + 0.05 作为注释位置
cell_number$y <- sapply(cell_number$cellid, function(g) {
  max(data[data$cellid == g, i]) + 0.05})

# 保证分组顺序
data$cellid <- factor(data$cellid, levels = c("NS_P3", "NS_P7", "NS_P12"))

# 组间比较设置
comparisons <- list(c("NS_P3", "NS_P7"), c("NS_P7", "NS_P12"), c("NS_P3", "NS_P12"))

# 颜色
colors <- c('#507BA8', '#F38D37', '#F1CE60')

# 绘图
p <- ggplot(data, aes(x = cellid, y = data[[i]], fill = cellid, color = cellid)) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(size = 6, face = "plain", color = "black"),
    axis.text.y = element_text(size = 6, face = "plain", color = "black"),
    plot.title = element_blank(),
    axis.title.y = element_text(color = 'black', size = 8, face = "bold", vjust = 0.5),
    legend.position = 'none') +
  labs(x = NULL, y = NULL, title = "dendrite development") +
  geom_jitter_rast(col = "#00000033", pch = 19, cex = 0.8, stroke = 0.01, 
                   position = position_jitter(0.15), alpha = 0.3) +
  scale_fill_manual(values = colors) +
  geom_boxplot(color = 'black', outlier.shape = NA, alpha = 0.8, size = 0.5, width = 0.4) +
  stat_compare_means(comparisons = comparisons, method = "t.test", 
                     label = "p.signif", size = 4, vjust = -0.5) +
  geom_text(data = cell_number, aes(x = cellid, y = y, label = number),
            inherit.aes = FALSE, size = 2.5)

# 保存图像
ggsave("FIGURE/dendrite_development_score.pdf", p,
       width = 4.5, height = 3.5, units = "cm", dpi = 600)


-------------------------9 monocle--------------------------
  --------------------------10 scFEA--------------------------
