#### ライブラリの読み込み ####
library(Seurat)

pbmc <- readRDS("./cell_type_anotation.rds")

DEG <- FindMarkers(pbmc, ident.1="CD4+ T", ident.2="NK")

head(DEG)

write.csv(DEG, "DEG1_CD4+T_vs_NK.csv")
