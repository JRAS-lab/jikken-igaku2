#### ライブラリの読み込み ####
library(Seurat)

#### データの読み込み ####
# gene expression matrixの場所を指定して読み込む
pbmc.data <- Read10X(data.dir = "./gene_expression_matrix/")

# Seuratのデータ形式に変換
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)


#### low-quality cellのフィルタリング ####
# MT-から始まるミトコンドリアRNAを"percent.mt"列としてデータに追加
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# ミトコンドリアRNAの割合をバイオリンプロットで確認
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# ミトコンドリアRNAが5%以下、遺伝子数が200~4500の間の細胞のみ抽出
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)

#### データの正規化（normalization) ####
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#### クラスタリング ####
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label=TRUE)

#### PCAによる次元削減 ####
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#### マーカー遺伝子の分布を表示 ####
FeaturePlot(pbmc, features = c("IL7R", "CD4", "CCR7", 
                               "S100A4", "CD14", "CD8A",
                               "CD4", "FCGR3A", "NCAM1", 
                               "GNLY","MS4A1", "CD79A", 
                               "PPBP","FCER1A"), pt.size =0.0000001)



#### Cell Typeのアノテーション ####
pbmc <- RenameIdents(pbmc, 
                     `0` = "CD4+ T", `1` = "CD14+ monocyte", `2` = "NK",
                     `3` = "CD4+ T", `4` = "CD8+ T", `5` = "CD8+ T",
                     `6` = "B", `7` = "CD16+ monocyte", `8` = "CD4+ T",
                     `9` = "DC", `10` = "CD14+ monocyte", `11` = "Platelet", 
                     `12` = "CD4+ T", `13` = "Unknown")

DimPlot(pbmc, reduction = "umap", label=TRUE, repel=TRUE)
