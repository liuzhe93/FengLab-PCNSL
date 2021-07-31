# set workpath and load R packages
# conda activate Renv
setwd("/ifs1/Grp8/liuzhe/scRNA/")
rm(list = ls())
#remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
#remotes::install_github("jlmelville/uwot")
#install.packages('Seurat')
library("Seurat")
#install.packages('tidyverse')
library("tidyverse")
library("Matrix")
library("scales")
library("cowplot")
library("RCurl")
options(stringsAsFactors = F)
library("scRNAseq")
library("scater")
# devtools::install_github(repo = "satijalab/seurat", ref = "loom")
library("loomR")
library("patchwork")
library("scuttle")
library("dplyr")
library("tibble")
library("HCAData")
library("SingleR")
library("org.Hs.eg.db")
library("clusterProfiler")
library("vroom")
library("celldex")
library("dittoSeq")
set.seed(2020)
library("monocle3")
suppressPackageStartupMessages(library("SingleCellExperiment"))
library("ggplot2"); theme_set(theme_bw())
library("DuoClustering2018")
require(scry)
library("scran")
library("DoubletFinder")
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
# load the scRNA-seq data
# create each individual Seurat object for every sample
# create Seurat object cancer and normal
seurat_data<-Read10X(data.dir=paste0("0_counts/","cancer"))
#seurat_obj<-CreateSeuratObject(counts=seurat_data,project="cancer")
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="cancer")
assign("cancer",seurat_obj)
cancer<-NormalizeData(cancer)
cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 3000)
cancer <- ScaleData(cancer)
cancer <- RunPCA(cancer)
cancer <- RunUMAP(cancer, reduction = "pca", dims = 1:21)
cancer <- RunTSNE(cancer, reduction = "pca", dims = 1:21)
## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
cancer <- SCTransform(cancer)
cancer <- RunPCA(cancer)
cancer <- RunUMAP(cancer, reduction = "pca", dims = 1:21)
cancer <- RunTSNE(cancer, reduction = "pca", dims = 1:21)
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_brain <- paramSweep_v3(cancer, PCs = 1:21, sct = FALSE)
sweep.stats_brain <- summarizeSweep(sweep.res.list_brain, GT = FALSE)
head(sweep.stats_brain)
bcmvn_brain <- find.pK(sweep.stats_brain)
mpK<-as.numeric(as.vector(bcmvn_brain$pK[which.max(bcmvn_brain$BCmetric)]))
DoubletRate = ncol(cancer)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
#DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
cancer <- FindNeighbors(cancer, reduction = "pca", dims = 1:21)
cancer <- FindClusters(cancer, resolution = 0.6)
levels(cancer)
cancer <- RunUMAP(cancer, reduction = "pca", dims = 1:21)
cancer <- RunTSNE(cancer, reduction = "pca", dims = 1:21)
annotations <- cancer@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(DoubletRate*nrow(cancer@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
cancer <- doubletFinder_v3(cancer, PCs = 1:21, pN = 0.05, pK = 5e-04, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
cancer <- doubletFinder_v3(cancer, PCs = 1:21, pN = 0.05, pK = 5e-04, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.05_5e-04_9717", sct = FALSE)
pdf("test_Double.pdf", width = 10, height = 10)
DimPlot(cancer, reduction = "tsne", group.by = "DF.classifications_0.05_5e-04_8535")
dev.off()
cancer@meta.data$singledouble<-cancer@meta.data$'DF.classifications_0.05_5e-04_8535'
cancer.singlet <- subset(cancer, subset = singledouble == "Singlet")
cancer$'DF.classifications_0.05_5e-04_8535'<-Idents(cancer)
pdf("test_Double.pdf", width = 10, height = 10)
DimPlot(cancer, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8, group.by = 'DF.classifications_0.05_5e-04_8535')
dev.off()
cancer_singlet<-cancer.singlet
pdf("1_preoperation/figures/pre-operation/cancer.counts.vs.features.pdf")
plot(x=cancer_singlet@meta.data$nCount_RNA,y=cancer_singlet@meta.data$nFeature_RNA)
dev.off()
# check the metadata in the new Seurat objects
head(cancer_singlet@meta.data)
tail(cancer_singlet@meta.data)
# Create .RData object to load at any time
save(cancer_singlet, file="1_preoperation/data/cancer.combined.RData")
cancer_singlet$log10GenesPerUMI <- log10(cancer_singlet$nFeature_RNA) / log10(cancer_singlet$nCount_RNA)
cancer_singlet$mitoRatio <- PercentageFeatureSet(object = cancer_singlet, pattern = "^MT-")
cancer_singlet$mitoRatio <- cancer_singlet@meta.data$mitoRatio / 100
cancermetadata <- cancer_singlet@meta.data
cancermetadata$cells <- rownames(cancermetadata)
cancermetadata <- cancermetadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
cancer_singlet
cancer_singlet@meta.data <- cancermetadata
counts <- GetAssayData(object = cancer_singlet, slot = "counts")
cancer_singlet <- CreateSeuratObject(counts, meta.data = cancer@meta.data)
cancer_singlet$label <- "cancer"
cancer_norm <- NormalizeData(cancer_singlet, normalization.method = "LogNormalize", scale.factor = 10000)
cancer_norm <- FindVariableFeatures(cancer_norm, selection.method = "vst", nfeatures = 3000)
pdf("1_preoperation/figures/pre-operation/cancer_Visualize_QC.pdf", width = 12, height = 6)
VlnPlot(cancer_singlet, features = c("nFeature_SCT", "nCount_SCT", "mitoRatio"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(cancer_singlet, feature1 = "nCount_SCT", feature2 = "mitoRatio")
plot2 <- FeatureScatter(cancer_singlet, feature1 = "nCount_SCT", feature2 = "nFeature_SCT")
pdf("1_preoperation/figures/pre-operation/cancer_FeatureScatter.pdf", width = 12, height = 6)
CombinePlots(plots = list(plot1, plot2))
dev.off()
top30 <- head(VariableFeatures(cancer_norm), 30)
pdf("1_preoperation/figures/pre-operation/cancer_VariableFeatures.pdf", width = 12, height = 6)
plot1 <- VariableFeaturePlot(cancer_norm)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
# filter cancer
filtered_cancer <- subset(x = cancer_singlet, 
                          subset= (nUMI < 20000) & 
                            (nGene > 200) &
                            (nGene < 2000) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.1))
filtered_cancer
cancer_singlet[["percent_HBA1"]] <- PercentageFeatureSet(filtered_cancer, features = "HBA1")
cancer_singlet[["percent_HBB"]] <- PercentageFeatureSet(filtered_cancer, features = "HBB")
cancer_singlet <- subset(x = cancer_singlet, subset= percent_HBA1 < 0.001 | percent_HBB < 0.001)
counts <- GetAssayData(object = filtered_cancer, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_cancer <- CreateSeuratObject(filtered_counts, meta.data = cancer_singlet@meta.data)
filtered_cancer$label <- "cancer"
save(filtered_cancer, file="1_preoperation/data/filtered_cancer.RData")
filtered_cancer_norm<-NormalizeData(filtered_cancer)

setwd("/ifs1/Grp8/liuzhe/scRNA/")
seurat_data<-Read10X(data.dir=paste0("0_counts/normal/","HFA567_total.filtered_gene_matrices"))
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="normal")
assign("normal1",seurat_obj)
normal1<-NormalizeData(normal1)
seurat_data<-Read10X(data.dir=paste0("0_counts/normal/","HFA570_total.filtered_gene_matrices"))
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="normal")
assign("normal2",seurat_obj)
normal2<-NormalizeData(normal2)
seurat_data<-Read10X(data.dir=paste0("0_counts/normal/","HFA571_total.filtered_gene_matrices"))
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="normal")
assign("normal3",seurat_obj)
normal3<-NormalizeData(normal3)
normal.normalized.combined <- merge(normal1, y = c(normal2, normal3), add.cell.ids = c("N1", "N2", "N3"), project = "normal", merge.data = TRUE)
normal<-normal.normalized.combined
pdf("1_preoperation/figures/pre-operation/normal.counts.vs.features.pdf")
plot(x=normal@meta.data$nCount_RNA,y=normal@meta.data$nFeature_RNA)
dev.off()
# check the metadata in the new Seurat objects
head(normal@meta.data)
tail(normal@meta.data)
# Create .RData object to load at any time
save(normal, file="1_preoperation/data/normal.combined.RData")
normal$log10GenesPerUMI <- log10(normal$nFeature_RNA) / log10(normal$nCount_RNA)
normal$mitoRatio <- PercentageFeatureSet(object = normal, pattern = "^MT-")
normal$mitoRatio <- normal@meta.data$mitoRatio / 100
normalmetadata <- normal@meta.data
normalmetadata$cells <- rownames(normalmetadata)
normalmetadata <- normalmetadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
normal
normal@meta.data <- normalmetadata
counts <- GetAssayData(object = normal, slot = "counts")
normal <- CreateSeuratObject(counts, meta.data = normal@meta.data)
normal$label <- "normal"
normal_norm <- NormalizeData(normal, normalization.method = "LogNormalize", scale.factor = 10000)
normal_norm <- FindVariableFeatures(normal_norm, selection.method = "vst", nfeatures = 3000)
pdf("1_preoperation/figures/pre-operation/normal_Visualize_QC.pdf", width = 12, height = 6)
VlnPlot(normal, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(normal, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(normal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("1_preoperation/figures/pre-operation/normal_FeatureScatter.pdf", width = 12, height = 6)
CombinePlots(plots = list(plot1, plot2))
dev.off()
top30 <- head(VariableFeatures(normal_norm), 30)
pdf("1_preoperation/figures/pre-operation/normal_VariableFeatures.pdf", width = 12, height = 6)
plot1 <- VariableFeaturePlot(normal_norm)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()
# filter normal
filtered_normal <- subset(x = normal, 
                          subset= (nUMI < 20000) & 
                            (nGene > 100) &
                            (nGene < 3000) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.1))
counts <- GetAssayData(object = filtered_normal, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_normal <- CreateSeuratObject(filtered_counts, meta.data = normal@meta.data)
filtered_normal$label <- "normal"
save(filtered_normal, file="1_preoperation/data/filtered_normal.RData")
filtered_normal_norm<-NormalizeData(filtered_normal)
WBY.anchors <- FindIntegrationAnchors(object.list = list(filtered_cancer_norm, filtered_normal_norm), dims = 1:30)
save(WBY.anchors, file="1_preoperation/data/integrated.anchors_seurat20210604.RData")
WBY.combined <- IntegrateData(anchorset = WBY.anchors, dims = 1:30)
WBY.combined <- FindVariableFeatures(WBY.combined, selection.method = "vst", nfeatures = 3000)
save(WBY.combined, file="1_preoperation/data/integrated.combined_seurat20210510.RData")

DefaultAssay(WBY.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
WBY.combined <- ScaleData(WBY.combined, verbose = FALSE, vars.to.regress = c("nUMI", "mitoRatio"))
#Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA. Therefore, the default in ScaleData() is only to perform scaling on the previously identified variable features (2,000 by default). 
WBY.combined <- ScaleData(WBY.combined)
save(WBY.combined, file="1_preoperation/data/WBY.combined_scaled.RData")
WBY.combined <- RunPCA(WBY.combined, npcs = 30, verbose = FALSE)
print(WBY.combined[["pca"]], dims = 1:5, nfeatures = 5)
pdf("1_preoperation/data/VizDimLoadings.pdf")
VizDimLoadings(WBY.combined, dims = 1:2, reduction = "pca")
dev.off()
pdf("1_preoperation/data/DimPlot.pdf")
DimPlot(WBY.combined, reduction = "pca")
dev.off()
pdf("1_preoperation/data/DimHeatmap.pc1.pdf")
DimHeatmap(WBY.combined, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf("1_preoperation/data/DimHeatmap.all.pdf")
DimHeatmap(WBY.combined, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
WBY.combined <- JackStraw(WBY.combined, num.replicate = 100, dims = 30)
WBY.combined <- ScoreJackStraw(WBY.combined, dims = 1:30)
pdf("1_preoperation/data/Determine_dimensionality.pdf", width = 24, height = 18)
p1 <- JackStrawPlot(WBY.combined, dims = 1:30)
p2 <- ElbowPlot(WBY.combined,ndims = 30)
plot_grid(p1, p2)
dev.off()
save(WBY.combined, file="1_preoperation/data/integrated.combined_beforepcs.RData")
# determine the resolution
library(Seurat)
library(clustree)
WBY.combined <- FindNeighbors(WBY.combined, reduction = "pca", dims = 1:21)
WBY.combined <- FindClusters(WBY.combined, resolution = 0.5)
WBY.combined <- RunUMAP(WBY.combined, reduction = "pca", dims = 1:21)
WBY.combined <- RunTSNE(WBY.combined, reduction = "pca", dims = 1:21)
save(WBY.combined, file="1_preoperation/data/integrated.combined_beforepcs.RData")

#sce <- as.SingleCellExperiment(WBY.combined)
#sce <- FindNeighbors(sce, dims = 1:21)
#save(sce, file="1_preoperation/data/integrated.combined.sce_beforepcs.RData")
#sce <- FindClusters(
#  object = sce,
#  resolution = c(seq(.1,1.6,.2))
#)
#pdf("1_preoperation/data/clusters.pdf",width=30,height=15)
#clustree(sce@meta.data, prefix = "integrated_snn_res.")
#colnames(sce@meta.data)
#dev.off()

#DefaultAssay(WBY.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
#WBY.combined<-FindVariableFeatures(WBY.combined, selection.method = "vst", nfeatures = 3000)
#WBY.combined <- ScaleData(WBY.combined)
#save(WBY.combined, file="1_preoperation/data/WBY.combined_scaled.RData")
#WBY.combined <- RunPCA(WBY.combined, npcs = 30, verbose = FALSE)
#print(WBY.combined[["pca"]], dims = 1:5, nfeatures = 5)
#pdf("1_preoperation/figures/pre-operation/VizDimLoadings.pdf")
#VizDimLoadings(WBY.combined, dims = 1:2, reduction = "pca")
#dev.off()
#pdf("1_preoperation/figures/pre-operation/DimPlot.pdf")
#DimPlot(WBY.combined, reduction = "pca")
#dev.off()
#pdf("1_preoperation/figures/pre-operation/DimHeatmap.pc1.pdf")
#DimHeatmap(WBY.combined, dims = 1, cells = 500, balanced = TRUE)
#dev.off()
#pdf("1_preoperation/figures/pre-operation/DimHeatmap.all.pdf")
#DimHeatmap(WBY.combined, dims = 1:30, cells = 500, balanced = TRUE)
#dev.off()
#WBY.combined <- JackStraw(WBY.combined, num.replicate = 100, dims = 30)
#WBY.combined <- ScoreJackStraw(WBY.combined, dims = 1:30)
#pdf("1_preoperation/figures/pre-operation/Determine_dimensionality.pdf", width = 24, height = 18)
#p1 <- JackStrawPlot(WBY.combined, dims = 1:30)
#p2 <- ElbowPlot(WBY.combined,ndims = 30)
#plot_grid(p1, p2)
#dev.off()
#save(WBY.combined, file="1_preoperation/data/integrated.combined_beforepcs.RData")
# determine the resolution
r=0.7
WBY.combined <- FindClusters(WBY.combined, resolution = r)
levels(WBY.combined)
save(WBY.combined,file="2_annotation/seurat/WBY.combined.res0.7.RData")
WBY.pca21.markers <- FindAllMarkers(object = WBY.combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
save(WBY.pca21.markers,file="2_annotation/seurat/WBY.pca21.markers.res0.7.RData")
top30<-WBY.pca21.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_annotation/seurat/top30.pca21.markers.csv",sep=",",quote=F)
head(Idents(WBY.combined), 5)
write.table(WBY.pca21.markers,"2_annotation/seurat/WBY.pca21.markers.csv",sep=",",quote=F)

#library(tidyverse)
#umap_tx = WBY.combined.pca25@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(label = WBY.combined.pca25@meta.data$label) %>% cbind(subpop = WBY.combined.pca25@meta.data$integrated_snn_res.0.5)  %>% cbind(indiv = WBY.combined.pca25@meta.data$seq_folder)
#write.csv(umap_tx,file="2_scRNA-seq/annotation/Seurat/umap_tx.csv",quote=F)
#tsne_tx = WBY.combined.pca25@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(label = WBY.combined.pca25@meta.data$label) %>% cbind(subpop = WBY.combined.pca25@meta.data$integrated_snn_res.0.5)  %>% cbind(indiv = WBY.combined.pca25@meta.data$seq_folder)
#write.csv(tsne_tx,file="2_scRNA-seq/annotation/Seurat/tsne_tx.csv",quote=F)
#umap_tx_ind = WBY.combined.pca25@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(tx = WBY.combined.pca25@meta.data$seq_folder)
#write.csv(umap_tx_ind,file="annotation/Seurat/umap_tx_ind.csv",quote=F)
#ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2, color=tx)) + geom_point() + 
#scale_color_manual(values=c("group1_untreated" = "darkblue", 
#                            "group1_treated" = "darkred"))
#tsne_tx = WBY.combined.pca25@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(tx = WBY.combined.pca25@meta.data$label)
#tsne_tx_ind = WBY.combined.pca25@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(tx = WBY.combined.pca25@meta.data$seq_folder)

# Visualization
pdf(paste0("2_annotation/seurat/umap.pca21.res",r,".splitbyLabel.pdf"),width=20,height=10)
DimPlot(WBY.combined, reduction = "umap", label = TRUE, pt.size=1,label.size = 8, split.by = 'label', group.by = 'integrated_snn_res.0.7')
dev.off()
pdf(paste0("2_annotation/seurat/umap.pca21.res",r,".pdf"),width=10,height=10)
DimPlot(WBY.combined, reduction = "umap", label = TRUE, pt.size=1,label.size = 8, group.by = 'integrated_snn_res.0.7')
dev.off()
pdf(paste0("2_annotation/seurat/tsne.pca21.res",r,".splitbyLabel.pdf"),width=20,height=10)
DimPlot(WBY.combined, reduction = "tsne", label = TRUE, pt.size=0.1,label.size = 8, split.by = 'label', group.by = 'integrated_snn_res.0.7')
dev.off()
pdf(paste0("2_annotation/seurat/tsne.pca21.res",r,".pdf"),width=10,height=10)
DimPlot(WBY.combined, reduction = "tsne", label = TRUE, pt.size=0.1,label.size = 8, group.by = 'integrated_snn_res.0.7')
dev.off()

prop.table(table(Idents(WBY.combined), WBY.combined$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(WBY.combined), WBY.combined$label)))
prop.table(table(Idents(WBY.combined)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(WBY.combined))))
write.csv(x = allsampleprop.each,file = '2_annotation/seurat/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = '2_annotation/seurat/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(WBY.combined))
pro.total <- table(Idents(WBY.combined),WBY.combined$label)
table(Idents(WBY.combined),WBY.combined$label)
pro.each <- table(Idents(WBY.combined),WBY.combined$label)
write.csv(x =pro.total,file = '2_annotation/seurat/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = '2_annotation/seurat/anno.pro.each.csv',quote = T,row.names = T)
save(WBY.combined,file="2_annotation/seurat/WBY.combined.pca21.res0.7.17clusters.aftercluster.autoSeurat.nolabel.RData")
save(WBY.combined,file="2_annotation/seurat/WBY.combined.pca21.aftercluster.autoSeurat.nolabel.RData")
# annotation
new.cluster.ids<-c("B_cell", "T_cell", "EN", "IN", "IN", "RG", "Macrophage", "EN", "INP", "Dendritic_cell", "ENP", "Astrocyte", "Oligodendrocyte", "Miningeal_cell", "OPC", "EN")
names(new.cluster.ids) <- levels(WBY.combined)
WBY.combined <- RenameIdents(WBY.combined, new.cluster.ids)
WBY.combined$celltype<-Idents(WBY.combined)
save(WBY.combined,file="/ifs1/Grp8/liuzhe/scRNA/2_annotation/seurat/WBY.combined.pca21.res0.7.afteranno.RData")
WBY.combined.markers <- FindAllMarkers(object = WBY.combined, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(WBY.combined.markers,"2_annotation/seurat/WBY.pca21.markers.csv",sep=",",quote=F)
top30<-WBY.combined.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_annotation/seurat/top30.pca21.markers.csv",sep=",",quote=F)
save(WBY.combined.markers,file="/ifs1/Grp8/liuzhe/scRNA/2_annotation/seurat/WBY.combined.pca21.res0.7.marker.afteranno.RData")
# calculate the percentage and count the numbers
prop.table(table(Idents(WBY.combined), WBY.combined$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(WBY.combined), WBY.combined$label)))
prop.table(table(Idents(WBY.combined)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(WBY.combined))))
write.csv(x = allsampleprop.each,file = '2_annotation/seurat/anno.allsample.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = '2_annotation/seurat/anno.allsample.total.prop.csv',quote = T,row.names = T)
table(Idents(WBY.combined))
pro.total <- table(Idents(WBY.combined),WBY.combined$label)
table(Idents(WBY.combined),WBY.combined$label)
pro.each <- table(Idents(WBY.combined),WBY.combined$label)
write.csv(x =pro.total,file = '2_annotation/seurat/anno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = '2_annotation/seurat/anno.pro.each.csv',quote = T,row.names = T)
# visualization
WBY.combined$celltype <- Idents(WBY.combined)
pdf(paste0("2_annotation/seurat/umap.pca21.res",r,".splitbyLabel.pdf"),width=20,height=10)
DimPlot(WBY.combined, reduction = "umap", label = TRUE, pt.size=1.5,label.size = 8, split.by = 'label', group.by = "celltype")
dev.off()
pdf(paste0("2_annotation/seurat/umap.pca21.res",r,".pdf"),width=10,height=10)
DimPlot(WBY.combined, reduction = "umap", label = TRUE, pt.size=1.5,label.size = 8, group.by = 'celltype')
dev.off()
pdf(paste0("2_annotation/seurat/tsne.pca21.res",r,".splitbyLabel.pdf"),width=20,height=10)
DimPlot(WBY.combined, reduction = "tsne", label = TRUE, pt.size=1,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
pdf(paste0("2_annotation/seurat/tsne.pca21.res",r,".pdf"),width=10,height=10)
DimPlot(WBY.combined, reduction = "tsne", label = TRUE, pt.size=1,label.size = 8, group.by = 'label')
dev.off()

delta.genes <- c("CD79A","CD37","CCL5","CRIP1","ENC1","EGR1","HBA2","CXCR4","MSMO1","RPS26","LYZ","C1QA","HIST1H4C","TOP2A","CST3","CPVL","EOMES","CORO1C","TTYH1","FGFBP3","PLP1","MBP","TIMP1","IGFBP7","APOD","S100B")
pdf("2_annotation/anno/dittoDotPlot.pdf",width=15,height=8)
dittoDotPlot(WBY.combined, vars = delta.genes, group.by = "celltype",scale = FALSE)
dev.off()
pdf("2_annotation/anno/heatmap.top30.pdf",width=24,height=18)
DoHeatmap(WBY.combined,features=top30$gene,cells = 1:500, size = 4, angle = 90, disp.min=-2, disp.max=2) + scale_fill_gradientn(colours=c("blue","white","red"))
dev.off()
library(psych)
library(pheatmap)
AverageExp<-AverageExpression(WBY.combined,features=unique(top30$gene))
typeof(AverageExp)
head(AverageExp$RNA)
DefaultAssay(WBY.combined) <- "integrated"
pdf("2_annotation/anno/averageExptop30.clusters.pdf")
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r,cluster_row = FALSE,cluster_col = FALSE)
dev.off()
pdf("2_annotation/anno/heatmap.top30.pdf",width=24,height=18)
DoHeatmap(WBY.combined,features=top30$gene,cells = 1:500, size = 4, angle = 90, disp.min=-2, disp.max=2) + scale_fill_gradientn(colours=c("blue","white","red"))
dev.off()
DefaultAssay(WBY.combined) <- "RNA"
features.plot <- c("CD37","CD79A","EEF1B2","CCL5","CRIP1","TRAC","ENC1","SLA","NRP1","PLS3","PDE4DIP","MEG3","MSMO1","FDFT1","TSPAN13","LYZ","FTL","FTH1","TOP2A","UBE2C","ZWINT","CST3","SNX3","VIM","EOMES","CORO1C","ADGRG1","PON2","CLU","GFAP","PLP1","CRYAB","CLDN11","MGP","C1S","OLIG1","PLLP","CMTM5")

pdf("2_annotation/anno/markergenes.dotplot.pdf",width = 10, height = 8)
DotPlot(object = WBY.combined, features = features.plot,  cols = c("lightgrey", "red"))
dev.off()

pdf("2_annotation/anno/markergenes.dotplot.pdf",width = 14, height = 8)
DotPlot(object = WBY.combined, features=features.plot,dot.scale = 10,cols = c("gray90", "red")) + RotatedAxis()
dev.off()


pdf("2_annotation/anno/markergenes.Bcell.pdf", width = 16, height = 8)
p1 <- FeaturePlot(WBY.combined, features = c("CD37", "CD79A"), pt.size = 1.5, combine = FALSE, reduction="tsne" )
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(2, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()
pdf("2_annotation/anno/markergenes.MPC.pdf", width = 16, height = 8)
p1 <- FeaturePlot(WBY.combined, features = c("TUBA1B", "HMGB2"), pt.size = 1.5, combine = FALSE, reduction="tsne" )
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(3, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()





WBY.combined_cancer<-subset(x=WBY.combined,subset = label == "cancer")
cellcom <- subset(WBY.combined_cancer, subset = (celltype == "B_cell" | celltype == "T_cell" | celltype == "Macrophage" | celltype == "Dendritic_cell"))
write.table(as.matrix(cellcom@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(cellcom@meta.data), cellcom@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  ????????в?????NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)





write.table(as.matrix(WBY.combined@assays$RNA@data), '3_cellphone/data/cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(WBY.combined@meta.data), WBY.combined@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  ????????в?????NA
write.table(meta_data, '3_cellphone/data/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)


require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(WBY.combined))
names(geneList)<-row.names(WBY.combined)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(WBY.combined.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(WBY.combined.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="2_annotation/anno/gene_erichment_results.RData")


write.csv(gene_erichment_results[["B_cell"]][["topGO"]],"2_annotation/anno/GO_B_cell.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["T_cell"]][["topGO"]],"2_annotation/anno/GO_T_cell.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["EN"]][["topGO"]],"2_annotation/anno/GO_EN.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["IN"]][["topGO"]],"2_annotation/anno/GO_IN.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["RG"]][["topGO"]],"2_annotation/anno/GO_RG.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["Macrophage"]][["topGO"]],"2_annotation/anno/GO_Macrophage.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["INP"]][["topGO"]],"2_annotation/anno/GO_INP.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["Dendritic_cell"]][["topGO"]],"2_annotation/anno/GO_Dendritic_cell.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["ENP"]][["topGO"]],"2_annotation/anno/GO_ENP.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["Astrocyte"]][["topGO"]],"2_annotation/anno/GO_Astrocyte.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["Oligodendrocyte"]][["topGO"]],"2_annotation/anno/GO_Oligodendrocyte.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["Miningeal_cell"]][["topGO"]],"2_annotation/anno/GO_Minigeal_cell.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["OPC"]][["topGO"]],"2_annotation/anno/GO_OPC.GO.csv",quote=F,row.names=F)



#MPC
gene_erichment_results[["Multilymphoid_progenitor_cell"]][["topGO"]][1:5,]
pdf("2_annotation/anno/clusterAnalysis_MPC.1.pdf",width=8,height=10)
library(enrichplot)
dotplot(gene_erichment_results[["Multilymphoid_progenitor_cell"]][["DGN"]], showCategory=30) 
dev.off()

#Bcell
gene_erichment_results[["B_cell"]][["topGO"]][1:5,]
pdf("2_annotation/anno/clusterAnalysis_B_cell.1.pdf",width=8,height=10)
library(enrichplot)
dotplot(gene_erichment_results[["B_cell"]][["DGN"]], showCategory=30) 
dev.off()


WBY.combined$celltype<-Idents(WBY.combined)
WBY_cancer<-subset(x=WBY.combined,subset = label == "cancer")
B_cell.subset <- subset(WBY_cancer@meta.data, celltype=="B_cell")
scRNAsub.tme <- subset(WBY_cancer, cells=row.names(B_cell.subset))
scRNAsub.tme <- FindVariableFeatures(scRNAsub.tme, selection.method = "vst", nfeatures = 2000)
scale.genes.tme <-  rownames(scRNAsub.tme)
scRNAsub.tme <- ScaleData(scRNAsub.tme, features = scale.genes.tme)
scRNAsub.tme <- RunPCA(scRNAsub.tme, features = VariableFeatures(scRNAsub.tme))
pdf("2_annotation/subcluster/Determine.tme.pcnumber.pdf")
ElbowPlot(scRNAsub.tme, ndims=20, reduction="pca")
dev.off()
pc.num=1:8
##???????
scRNAsub.tme <- FindNeighbors(scRNAsub.tme, dims = pc.num) 
scRNAsub.tme <- FindClusters(scRNAsub.tme, resolution = 0.6)
table(scRNAsub.tme@meta.data$seurat_clusters)
metadata <- scRNAsub.tme@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'2_annotation/subcluster/tme.cell_cluster.csv',row.names = F)
##????????
#tSNE
scRNAsub.tme = RunTSNE(scRNAsub.tme, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub.tme, 'tsne')
write.csv(embed_tsne,'2_annotation/subcluster/tme.embed_tsne.csv')
pdf("2_annotation/subcluster/tsne_Bcell.pdf")
DimPlot(scRNAsub.tme, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8)
dev.off()
diff.wilcox = FindAllMarkers(scRNAsub.tme)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "2_annotation/subcluster/tme.diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "2_annotation/subcluster/tme.top30_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub.tme,file="2_annotation/subcluster/scRNAsub.tme.RData")
new.cluster.ids<-c("B_cell-1", "B_cell-2", "B_cell-3", "B_cell-3", "B_cell-1", "B_cell-2", "Plasma_cell")
names(new.cluster.ids) <- levels(scRNAsub.tme)
scRNAsub.tme <- RenameIdents(scRNAsub.tme, new.cluster.ids)
scRNAsub.tme <- RunUMAP(scRNAsub.tme, dims = 1:8)
scRNAsub.tme$celltype<-Idents(scRNAsub.tme)
save(scRNAsub.tme,file="2_annotation/subcluster/scRNAsub.tme.afteranno.RData")
scRNAsub.tme.markers <- FindAllMarkers(object = scRNAsub.tme, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.tme.markers,"2_annotation/subcluster/scRNAsub.tme.markers.anno.csv",sep=",",quote=F)
save(scRNAsub.tme.markers,file="2_annotation/subcluster/scRNAsub.tme.markers.afteranno.RData")
pdf("2_annotation/subcluster/tsne.bcell.integrate.pdf", width = 10, height = 10)
DimPlot(scRNAsub.tme, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
top30<-scRNAsub.tme.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_annotation/subcluster/scRNAsub.tme.top30.anno.csv",sep=",",quote=F)

write.table(as.matrix(scRNAsub.tme@assays$RNA@data), 'cellphonedb_Bcell_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(scRNAsub.tme@meta.data), scRNAsub.tme@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  ????????в?????NA
write.table(meta_data, 'cellphonedb_Bcell_meta.txt', sep='\t', quote=F, row.names=F)



library("SingleR")
#refdata <- MonacoImmuneData()
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
#immu.se <- DatabaseImmuneCellExpressionData() 
testdata <- GetAssayData(scRNAsub.tme, slot="data")
clusters <- scRNAsub.tme@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, 
                    ref = list(HP = hpca.se , BP = bpe.se),
                    labels = list(hpca.se$label.main , bpe.se$label.main), 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts",de.method="wilcox")           
table(cellpred$labels)
pdf("2_annotation/subcluster/test.pdf", width=18 ,height=9)
plotScoreHeatmap(cellpred)
dev.off()
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"2_annotation/subcluster/tme.celltype_singleR.csv",row.names = F)
scRNAsub.tme@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAsub.tme@meta.data[which(scRNAsub.tme@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub.tme, group.by="celltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsub.tme, group.by="celltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("2_annotation/subcluster/tme.tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("2_annotation/subcluster/tme.UMAP_celltype.pdf", p2, width=10 ,height=6)
ggsave("2_annotation/subcluster/tme.celltype.pdf", p3, width=10 ,height=5)
ggsave("2_annotation/subcluster/tme.celltype.png", p3, width=10 ,height=5)



require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(scRNAsub.tme))
names(geneList)<-row.names(scRNAsub.tme)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(all.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(all.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="2_annotation/subcluster/Bcell_gene_erichment_results.RData")
library(enrichplot)

require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(scRNAsub.tme))
names(geneList)<-row.names(scRNAsub.tme)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(scRNAsub.tme.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(scRNAsub.tme.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="2_annotation/anno/Bcellsub.gene_erichment_results.RData")
gene_erichment_results[["B_cell-1"]][["topGO"]][1:5,]

write.csv(gene_erichment_results[["B_cell-1"]][["topGO"]],"2_annotation/subcluster/Bcellsub.B_cell-1.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["Plasma_cell"]][["topGO"]],"2_annotation/subcluster/Bcellsub.Plasma_cell.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["B_cell-2"]][["topGO"]],"2_annotation/subcluster/Bcellsub.B_cell-2.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["B_cell-3"]][["topGO"]],"2_annotation/subcluster/Bcellsub.B_cell-3.GO.csv",quote=F,row.names=F)


    
names(new.cluster.ids) <- levels(scRNAsub.tme)
scRNAsub.tme <- RenameIdents(scRNAsub.tme, new.cluster.ids)
scRNAsub.tme <- RunUMAP(scRNAsub.tme, dims = 1:8)
save(scRNAsub.tme,file="data/scRNAsub.tme.afteranno.RData")
scRNAsub.tme.markers <- FindAllMarkers(object = scRNAsub.tme, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.tme.markers,"subcluster/scRNAsub.tme.markers.anno.csv",sep=",",quote=F)
save(scRNAsub.tme.markers,file="subcluster/scRNAsub.tme.markers.afteranno.RData")
pdf("subcluster/scRNAsub.tme.integrate.pdf", width = 36, height = 18)
p1 <- DimPlot(scRNAsub.tme, reduction = "umap", group.by = "label")
p2<-DimPlot(scRNAsub.tme, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
prop.table(table(Idents(scRNAsub.tme), scRNAsub.tme$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(scRNAsub.tme), scRNAsub.tme$label)))
prop.table(table(Idents(scRNAsub.tme)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(scRNAsub.tme))))
write.csv(x = allsampleprop.each,file = '2_annotation/subcluster/Bcellprop.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = '2_annotation/subcluster/Bcellprop.prop.csv',quote = T,row.names = T)
table(Idents(scRNAsub.tme))
pro.total <- table(Idents(scRNAsub.tme),scRNAsub.tme$label)
table(Idents(scRNAsub.tme),scRNAsub.tme$label)
pro.each <- table(Idents(scRNAsub.tme),scRNAsub.tme$label)
write.csv(x =pro.total,file = '2_annotation/subcluster/Bcellanno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = '2_annotation/subcluster/Bcellanno.pro.each.csv',quote = T,row.names = T)





pdf("results/pca16res0.3/umap.pca16.res0.3.integrate.pdf", width = 36, height = 18)
p1 <- DimPlot(glioma.combined.pca16, reduction = "umap", group.by = "label")
p2<-DimPlot(glioma.combined.pca16, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

tme.markers <- FindAllMarkers(object = scRNAsub.tme, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
save(WBY.pca21.markers,file="2_annotation/seurat/WBY.pca21.markers.res0.7.RData")
top30<-WBY.pca21.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)

plot1 = DimPlot(scRNAsub.tme, reduction = "tsne") 
ggsave("2_annotation/subcluster/tme.tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("2_annotation/subcluster/tme.tSNE.png", plot = plot1, width = 8, height = 7)




features.plot <- c("CD37", "CD79A", "TUBA1B", "HMGB2")
pdf("2_annotation/anno/test.RidgePlot.pdf",width = 20, height = 20)
RidgePlot(object = WBY.combined, features = features.plot,ncol=3)
dev.off()



pdf("2_annotation/anno/markergenes.MPC.pdf")
FeaturePlot(WBY.combined, features = "TUBA1B", blend.threshold = 1,reduction="tsne")
dev.off()
pdf("2_annotation/anno/markergenes.Bcell.pdf")
FeaturePlot(WBY.combined, features = "CD79A", reduction="tsne")
dev.off()












##############################################################################################################################
# auto annotation by SingleR training dataset
WBY.combined@meta.data$cell.type <- Idents(WBY.combined)
test <- as.SingleCellExperiment(WBY.combined)
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
#immu.se <- DatabaseImmuneCellExpressionData()
sceESC<-LaMannoBrainData('human-es')
sceEmidBrain<-LaMannoBrainData('human-embryo')
sceIPSC<-LaMannoBrainData('human-ips')
sceESC<-sceESC[,!is.na(sceESC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceESC <- logNormCounts(sceESC)
sceEmidBrain<-sceEmidBrain[,!is.na(sceEmidBrain$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceEmidBrain <- logNormCounts(sceEmidBrain)
sceIPSC<-sceIPSC[,!is.na(sceIPSC$Cell_type)]
# SingleR() expects reference datasets to be normalized and log-transformed.
sceIPSC <- logNormCounts(sceIPSC)
Anno <- SingleR(test = test,
                ref = list(HP = hpca.se , BP = bpe.se, ESC=sceESC, midbrain=sceEmidBrain, iPSC=sceIPSC),
                labels = list(hpca.se$label.main , bpe.se$label.main, sceESC$Cell_type, sceEmidBrain$Cell_type, sceIPSC$Cell_type),
                method = "cluster",
                cluster = test$cell.type)
table(Anno$pruned.labels)
Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
#??细??注????息???????拥?Seurat??????去
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(WBY.combined)
WBY.combined <- RenameIdents(WBY.combined, new.cluster.ids)
head(Idents(WBY.combined), 5)
levels(Idents(WBY.combined))
WBY.combined$celltype<-Idents(WBY.combined)
# Run non-linear dimensional reduction (UMAP/tSNE)
WBY.combined <- RunUMAP(WBY.combined, reduction = "pca", dims = 1:21)
# Visualization
WBY.combined <- RunTSNE(WBY.combined, dims = 1:21)
pdf("/ifs1/Grp8/liuzhe/scRNA/2_annotation/singleR/umap.singleRanno.pdf",width=10,height=10)
DimPlot(WBY.combined, reduction = "umap",group.by ="celltype",pt.size=2)
dev.off()
pdf("/ifs1/Grp8/liuzhe/scRNA/2_annotation/singleR/tsne.singleRanno.pdf",width=10,height=10)
DimPlot(WBY.combined, reduction = "tsne",group.by ="celltype",pt.size=2)
dev.off()
pdf("/ifs1/Grp8/liuzhe/scRNA/2_annotation/singleR/tsne.singleRanno.split.pdf",width=20,height=10)
DimPlot(WBY.combined, reduction = "tsne", label = TRUE, pt.size=0.1,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()

###DEGs analysis

pdf("2_annotation/anno/markergenes.macrophage.pdf", width = 16, height = 8)
p1 <- FeaturePlot(WBY.combined, features = c("LYZ", "NPC2"), pt.size = 1.5, combine = FALSE, reduction="tsne" )
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(2, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()

pdf("2_annotation/anno/markergenes.Meningeal_cell.pdf", width = 16, height = 8)
p1 <- FeaturePlot(WBY.combined, features = c("TIMP1", "IGFBP7"), pt.size = 1.5, combine = FALSE, reduction="tsne" )
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(2, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()

pdf("2_annotation/anno/markergenes.Tcell.pdf", width = 16, height = 8)
p1 <- FeaturePlot(WBY.combined, features = c("LTB", "CCL5"), pt.size = 1.5, combine = FALSE, reduction="tsne" )
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(2, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()

pdf("2_annotation/anno/markergenes.CD8PlusDendritic_cell.pdf", width = 16, height = 8)
p1 <- FeaturePlot(WBY.combined, features = c("CST3", "CPVL"), pt.size = 1.5, combine = FALSE, reduction="tsne" )
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(2, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()

WBY.combined$celltype<-Idents(WBY.combined)
pdf("2_annotation/anno/markergenes.macrophage.vioplot.pdf",width=15,height=10)
VlnPlot(WBY.combined, features = c("CST3","NPC2"), group.by = 'celltype', pt.size = 0) 
dev.off()





macrophage.cells <- subset(WBY.combined, idents = "Macrophage")
Idents(macrophage.cells) <- "label"
avg.macrophage.cells <- as.data.frame(log1p(AverageExpression(macrophage.cells, verbose = FALSE)$RNA))
avg.macrophage.cells$gene <- rownames(avg.macrophage.cells)
write.table(avg.macrophage.cells,"2_annotation/DEGs/avg.macrophage.cells.csv",sep=",",quote=F)
WBY.combined$celltype.label <- paste(Idents(WBY.combined), WBY.combined$label, sep = "_")
WBY.combined$celltype <- Idents(WBY.combined)
Idents(WBY.combined) <- "celltype.label"
scRNAsub.macrophage.data <- subset(WBY.combined, subset = (celltype.label == "Macrophage_cancer" | celltype.label == "Macrophage_normal"))
write.table(GetAssayData(scRNAsub.macrophage.data),"macrophage.data.csv",sep=",",quote=F)


immune.response <- FindMarkers(WBY.combined, ident.1 = "Macrophage_cancer", ident.2 = "Macrophage_normal", verbose = FALSE)
head(immune.response, n = 15)
write.table(immune.response,"2_annotation/DEGs/macrophage.degs.csv",sep=",")

library(devtools)
install_github("jingshuw/SAVERX")



pdf("2_annotation/DEGs/macrophage.figure1.pdf",width = 20, height = 18)
FeaturePlot(WBY.combined, features = c("SPP1","EGR2","CXCL8","CCL3","IL1B","CX3CR1","IL1B","CMC1","PTGDS","ANXA6","DUSP4","CTNNB1","PTTG1"), reduction = "tsne", split.by = "label", max.cutoff = 3, 
            cols = c("grey", "red"))
dev.off()
pdf("2_annotation/DEGs/macrophage.figure2.pdf",width = 10, height = 18)
plots <- VlnPlot(WBY.combined, features = c("TLR4","TLR10","BTG2","NFKBIA","NFKBIZ"), split.by = "label", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("2_annotation/DEGs/macrophage.figure3.pdf")
genes.to.label = c("SPP1","EGR2","CXCL8","CCL3","IL1B","CX3CR1","IL1B","CMC1","PTGDS","ANXA6","DUSP4","CTNNB1","PTTG1")
p1 <- ggplot(avg.macrophage.cells, aes(cancer, normal)) + geom_point() + ggtitle("macrophage")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
dev.off()


WBY.combined$celltype<-Idents(WBY.combined)
pdf("2_annotation/anno/markergenes.macrophage.vioplot.pdf",width=15,height=10)
VlnPlot(WBY.combined, features = c("CST3","NPC2"), group.by = 'celltype', pt.size = 0) 
dev.off()

macrophage.cells <- subset(WBY.combined, idents = "Macrophage")
Idents(macrophage.cells) <- "label"
avg.macrophage.cells <- as.data.frame(log1p(AverageExpression(macrophage.cells, verbose = FALSE)$RNA))
avg.macrophage.cells$gene <- rownames(avg.macrophage.cells)
write.table(avg.macrophage.cells,"2_annotation/DEGs/avg.macrophage.cells.csv",sep=",",quote=F)
WBY.combined$celltype.label <- paste(Idents(WBY.combined), WBY.combined$label, sep = "_")
WBY.combined$celltype <- Idents(WBY.combined)
Idents(WBY.combined) <- "celltype.label"
immune.response <- FindMarkers(WBY.combined, ident.1 = "Macrophage_cancer", ident.2 = "Macrophage_normal", verbose = FALSE)
head(immune.response, n = 15)
write.table(immune.response,"2_annotation/DEGs/macrophage.degs.csv",sep=",")
degs_macrophage<-immune.response
degs_filtered<-degs_macrophage[(degs_macrophage$avg_log2FC>log2(1.5)|degs_macrophage$avg_log2FC<=-log2(1.5) ) & (degs_macrophage$p_val_adj<0.05),]
features_degs<-degs_filtered[order(degs_filtered$avg_log2FC),]
pdf("2_annotation/DEGs/macrophage.heatmap.pdf",width = 20, height = 18)
DoHeatmap(macrophage.cells, features = row.names(features_degs)) + NoLegend()
dev.off()
pdf("2_annotation/DEGs/macrophage.heatmap.pdf",width = 20, height = 18)
DoHeatmap(macrophage.cells, features = row.names(features_degs))
dev.off()



pdf("2_annotation/DEGs/macrophage.figure1.pdf",width = 20, height = 18)
FeaturePlot(WBY.combined, features = c("TLR4","TLR10","BTG2","NFKBIA","NFKBIZ"), reduction = "tsne", split.by = "label", max.cutoff = 3, 
            cols = c("grey", "red"))
dev.off()
pdf("2_annotation/DEGs/macrophage.figure2.pdf",width = 10, height = 18)
plots <- VlnPlot(WBY.combined, features = c("TLR4","TLR10","BTG2","NFKBIA","NFKBIZ"), split.by = "label", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("2_annotation/DEGs/macrophage.figure3.pdf")
genes.to.label = c("TLR4","TLR10","BTG2","NFKBIA","NFKBIZ")
p1 <- ggplot(avg.macrophage.cells, aes(cancer, normal)) + geom_point() + ggtitle("macrophage")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
dev.off()

### T cell

WBY.combined$celltype<-Idents(WBY.combined)
WBY_cancer<-subset(x=WBY.combined,subset = label == "cancer")
T_cell.subset <- subset(WBY_cancer@meta.data, celltype=="T_cell")
scRNAsub.T_cell <- subset(WBY_cancer, cells=row.names(T_cell.subset))
scRNAsub.T_cell <- FindVariableFeatures(scRNAsub.T_cell, selection.method = "vst", nfeatures = 2000)
scale.genes.T_cell <-  rownames(scRNAsub.T_cell)
scRNAsub.T_cell <- ScaleData(scRNAsub.T_cell, features = scale.genes.T_cell)
scRNAsub.T_cell <- RunPCA(scRNAsub.T_cell, features = VariableFeatures(scRNAsub.T_cell))
pdf("2_annotation/subcluster/Determine.T_cell.pcnumber.pdf")
ElbowPlot(scRNAsub.T_cell, ndims=20, reduction="pca")
dev.off()
pc.num=1:9
##???????
scRNAsub.T_cell <- FindNeighbors(scRNAsub.T_cell, dims = pc.num) 
scRNAsub.T_cell <- FindClusters(scRNAsub.T_cell, resolution = 0.8)
table(scRNAsub.T_cell@meta.data$seurat_clusters)
metadata <- scRNAsub.T_cell@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'2_annotation/subcluster/T_cell.cell_cluster.csv',row.names = F)
##????????
#tSNE
scRNAsub.T_cell = RunTSNE(scRNAsub.T_cell, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub.T_cell, 'tsne')
write.csv(embed_tsne,'2_annotation/subcluster/T_cell.embed_tsne.csv')
pdf("2_annotation/subcluster/tsne_T_cell.pdf")
DimPlot(scRNAsub.T_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8) + NoLegend()
dev.off()
diff.wilcox = FindAllMarkers(scRNAsub.T_cell)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "2_annotation/subcluster/T_cell.diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "2_annotation/subcluster/T_cell.top30_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub.T_cell,file="2_annotation/subcluster/scRNAsub.T_cell.RData")

library(SingleR)
#refdata <- MonacoImmuneData()
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
#immu.se <- DatabaseImmuneCellExpressionData() 
testdata <- GetAssayData(scRNAsub.T_cell, slot="data")
clusters <- scRNAsub.T_cell@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, 
                    ref = list(HP = hpca.se , BP = bpe.se),
                    labels = list(hpca.se$label.main , bpe.se$label.main), 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts",de.method="wilcox")           
table(cellpred$labels)
pdf("2_annotation/subcluster/test.pdf", width=18 ,height=9)
plotScoreHeatmap(cellpred)
dev.off()
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"2_annotation/subcluster/T_cell.celltype_singleR.csv",row.names = F)
scRNAsub.T_cell@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAsub.T_cell@meta.data[which(scRNAsub.T_cell@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub.T_cell, group.by="celltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsub.T_cell, group.by="celltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("2_annotation/subcluster/T_cell.tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("2_annotation/subcluster/T_cell.UMAP_celltype.pdf", p2, width=10 ,height=6)
ggsave("2_annotation/subcluster/T_cell.celltype.pdf", p3, width=10 ,height=5)
ggsave("2_annotation/subcluster/T_cell.celltype.png", p3, width=10 ,height=5)

diff.wilcox = FindAllMarkers(scRNAsub.T_cell)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(scRNAsub.T_cell))
names(geneList)<-row.names(scRNAsub.T_cell)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(all.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(all.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="2_annotation/subcluster/T_cell_gene_erichment_results.RData")

write.csv(gene_erichment_results[["T_helper"]][["topGO"]],"2_annotation/subcluster/Tcellsub.T_helper.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["NKT"]][["topGO"]],"2_annotation/subcluster/Tcellsub.NKT.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["Multilymphoid_progenitor_cell"]][["topGO"]],"2_annotation/subcluster/Tcellsub.MPC.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["T_cell"]][["topGO"]],"2_annotation/subcluster/Tcellsub.T_cell.csv",quote=F,row.names=F)

load("2_annotation/subcluster/scRNAsub.T_cell.RData")
new.cluster.ids<-c("T_helper", "NKT", "Multilymphoid_progenitor_cell", "NKT", "Multilymphoid_progenitor_cell","NA", "T_cell", "NA")
names(new.cluster.ids) <- levels(scRNAsub.T_cell)
scRNAsub.T_cell <- RenameIdents(scRNAsub.T_cell, new.cluster.ids)
scRNAsub.T_cell$celltype<-Idents(scRNAsub.T_cell )
scRNAsub.T_cell_rmNA <- subset(scRNAsub.T_cell, subset = (celltype == "T_helper" | celltype == "NKT" | celltype == "Multilymphoid_progenitor_cell" | celltype == "T_cell"))
scRNAsub.T_cell<-scRNAsub.T_cell_rmNA
scRNAsub.T_cell <- RunUMAP(scRNAsub.T_cell, dims = 1:8)
scRNAsub.T_cell$celltype<-Idents(scRNAsub.T_cell)
save(scRNAsub.T_cell,file="2_annotation/subcluster/scRNAsub.T_cell.afteranno.RData")
scRNAsub.T_cell.markers <- FindAllMarkers(object = scRNAsub.T_cell, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.T_cell.markers,"2_annotation/subcluster/scRNAsub.T_cell.markers.anno.csv",sep=",",quote=F)
save(scRNAsub.T_cell.markers,file="2_annotation/subcluster/scRNAsub.T_cell.markers.afteranno.RData")
pdf("2_annotation/subcluster/tsne.Tcell.integrate.pdf", width = 10, height = 10)
DimPlot(scRNAsub.T_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
top30<-scRNAsub.T_cell.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_annotation/subcluster/scRNAsub.T_cell.top30.anno.csv",sep=",",quote=F)

features.plot <- c("LTB","NEAT1","IL7R","DUSP1","CCL5","JUNB","UBE2C","STMN1","BIRC5","NCL","HNRNPA3","DUT")
pdf("2_annotation/subcluster/Tcellsub.dittoDotPlot.pdf",width = 9, height = 4)
DotPlot(object = scRNAsub.T_cell, features=features.plot,dot.scale = 12,cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()
prop.table(table(Idents(scRNAsub.T_cell), scRNAsub.T_cell$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(scRNAsub.T_cell), scRNAsub.T_cell$label)))
prop.table(table(Idents(scRNAsub.T_cell)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(scRNAsub.T_cell))))
write.csv(x = allsampleprop.each,file = '2_annotation/subcluster/Tcellprop.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = '2_annotation/subcluster/Tcellprop.prop.csv',quote = T,row.names = T)
table(Idents(scRNAsub.T_cell))
pro.total <- table(Idents(scRNAsub.T_cell),scRNAsub.T_cell$label)
table(Idents(scRNAsub.T_cell),scRNAsub.T_cell$label)
pro.each <- table(Idents(scRNAsub.T_cell),scRNAsub.T_cell$label)
write.csv(x =pro.total,file = '2_annotation/subcluster/Tcellanno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = '2_annotation/subcluster/Tcellanno.pro.each.csv',quote = T,row.names = T)

## Dendritic_cell

WBY.combined$celltype<-Idents(WBY.combined)
WBY_cancer<-subset(x=WBY.combined,subset = label == "cancer")
CD8PlusDendritic_cell.subset <- subset(WBY_cancer@meta.data, celltype=="Dendritic_cell")
scRNAsub.CD8PlusDendritic_cell <- subset(WBY_cancer, cells=row.names(CD8PlusDendritic_cell.subset))
scRNAsub.CD8PlusDendritic_cell <- FindVariableFeatures(scRNAsub.CD8PlusDendritic_cell, selection.method = "vst", nfeatures = 2000)
scale.genes.CD8PlusDendritic_cell <-  rownames(scRNAsub.CD8PlusDendritic_cell)
scRNAsub.CD8PlusDendritic_cell <- ScaleData(scRNAsub.CD8PlusDendritic_cell, features = scale.genes.CD8PlusDendritic_cell)
scRNAsub.CD8PlusDendritic_cell <- RunPCA(scRNAsub.CD8PlusDendritic_cell, features = VariableFeatures(scRNAsub.CD8PlusDendritic_cell))
pdf("2_annotation/subcluster/Determine.CD8PlusDendritic_cell.pcnumber.pdf")
ElbowPlot(scRNAsub.CD8PlusDendritic_cell, ndims=20, reduction="pca")
dev.off()
pc.num=1:18
##???????
scRNAsub.CD8PlusDendritic_cell <- FindNeighbors(scRNAsub.CD8PlusDendritic_cell, dims = pc.num) 
scRNAsub.CD8PlusDendritic_cell <- FindClusters(scRNAsub.CD8PlusDendritic_cell, resolution = 1)
table(scRNAsub.CD8PlusDendritic_cell@meta.data$seurat_clusters)
metadata <- scRNAsub.CD8PlusDendritic_cell@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'2_annotation/subcluster/CD8PlusDendritic_cell.cell_cluster.csv',row.names = F)
##????????
#tSNE
scRNAsub.CD8PlusDendritic_cell = RunTSNE(scRNAsub.CD8PlusDendritic_cell, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub.CD8PlusDendritic_cell, 'tsne')
write.csv(embed_tsne,'2_annotation/subcluster/CD8PlusDendritic_cell.embed_tsne.csv')
pdf("2_annotation/subcluster/tsne_CD8PlusDendritic_cell.pdf")
DimPlot(scRNAsub.CD8PlusDendritic_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8)
dev.off()
diff.wilcox = FindAllMarkers(scRNAsub.CD8PlusDendritic_cell)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "2_annotation/subcluster/CD8PlusDendritic_cell.diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "2_annotation/subcluster/CD8PlusDendritic_cell.top30_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub.CD8PlusDendritic_cell,file="2_annotation/subcluster/scRNAsub.CD8PlusDendritic_cell.RData")

library(SingleR)
#refdata <- MonacoImmuneData()
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
#immu.se <- DatabaseImmuneCellExpressionData() 
testdata <- GetAssayData(scRNAsub.CD8PlusDendritic_cell, slot="data")
clusters <- scRNAsub.CD8PlusDendritic_cell@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, 
                    ref = list(HP = hpca.se , BP = bpe.se),
                    labels = list(hpca.se$label.main , bpe.se$label.main), 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts",de.method="wilcox")           
table(cellpred$labels)
pdf("2_annotation/subcluster/test.pdf", width=18 ,height=9)
plotScoreHeatmap(cellpred)
dev.off()
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"2_annotation/subcluster/CD8PlusDendritic_cell.celltype_singleR.csv",row.names = F)
scRNAsub.CD8PlusDendritic_cell@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAsub.CD8PlusDendritic_cell@meta.data[which(scRNAsub.CD8PlusDendritic_cell@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNAsub.CD8PlusDendritic_cell, group.by="celltype", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNAsub.CD8PlusDendritic_cell, group.by="celltype", label=T, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave("2_annotation/subcluster/CD8PlusDendritic_cell.tSNE_celltype.pdf", p1, width=7 ,height=6)
ggsave("2_annotation/subcluster/CD8PlusDendritic_cell.UMAP_celltype.pdf", p2, width=10 ,height=6)
ggsave("2_annotation/subcluster/CD8PlusDendritic_cell.celltype.pdf", p3, width=10 ,height=5)
ggsave("2_annotation/subcluster/CD8PlusDendritic_cell.celltype.png", p3, width=10 ,height=5)

diff.wilcox = FindAllMarkers(scRNAsub.CD8PlusDendritic_cell)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(top30, "2_annotation/subcluster/CD8PlusDendritic_cell.top30_diff_genes_wilcox.csv", row.names = F)
prop.table(table(Idents(scRNAsub.CD8PlusDendritic_cell), scRNAsub.CD8PlusDendritic_cell$label))
allsampleprop.each <-as.data.frame(prop.table(table(Idents(scRNAsub.CD8PlusDendritic_cell), scRNAsub.CD8PlusDendritic_cell$label)))
prop.table(table(Idents(scRNAsub.CD8PlusDendritic_cell)))
allsampleprop.total <-as.data.frame(prop.table(table(Idents(scRNAsub.CD8PlusDendritic_cell))))
write.csv(x = allsampleprop.each,file = '2_annotation/subcluster/Denprop.each.prop.csv',quote = T,row.names = T)
write.csv(x = allsampleprop.total,file = '2_annotation/subcluster/Denprop.prop.csv',quote = T,row.names = T)
table(Idents(scRNAsub.CD8PlusDendritic_cell))
pro.total <- table(Idents(scRNAsub.CD8PlusDendritic_cell),scRNAsub.CD8PlusDendritic_cell$label)
table(Idents(scRNAsub.CD8PlusDendritic_cell),scRNAsub.CD8PlusDendritic_cell$label)
pro.each <- table(Idents(scRNAsub.CD8PlusDendritic_cell),scRNAsub.CD8PlusDendritic_cell$label)
write.csv(x =pro.total,file = '2_annotation/subcluster/Denanno.pro.total.csv',quote = T,row.names = T)
write.csv(x =pro.each,file = '2_annotation/subcluster/Denanno.pro.each.csv',quote = T,row.names = T)
require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(scRNAsub.CD8PlusDendritic_cell))
names(geneList)<-row.names(scRNAsub.CD8PlusDendritic_cell)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(all.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(all.markers,cluster==c1)$gene
  gene_erichment_results[[c1]]=list()
  testgeneList=geneList
  testgeneList[which(newwallgenes %in% testgenes)]= 1
  #gene_erichment_results=list()
  tab1=c()
  for(ont in c("BP","MF")){
    sampleGOdata<-suppressMessages(new("topGOdata",description="Simple session",ontology=ont,allGenes=as.factor(testgeneList),
                                       nodeSize=10,annot=annFUN.org,mapping="org.Hs.eg.db",ID="entrez"))
    resultTopGO.elim<-suppressMessages(runTest(sampleGOdata,algorithm="elim",statistic="Fisher"))
    
    resultTopGO.classic<-suppressMessages(runTest(sampleGOdata,algorithm="classic",statistic="Fisher"))
    tab1<-rbind(tab1,GenTable(sampleGOdata,Fisher.elim=resultTopGO.elim,Fisher.classic=resultTopGO.classic,orderBy="Fisher.elim",
                              topNodes=200))
  }
  gene_erichment_results[[c1]][["topGO"]]=tab1
  x<-suppressMessages(enrichDO(gene=names(testgeneList)[testgeneList==1],ont="DO",pvalueCutoff=1,pAdjustMethod="BH",universe=names(testgeneList),
                               minGSSize=5,maxGSSize=500,readable=T))
  gene_erichment_results[[c1]][["DO"]]=x
  dgn<-suppressMessages(enrichDGN(names(testgeneList)[testgeneList==1]))
  gene_erichment_results[[c1]][["DGN"]]=dgn
}
save(gene_erichment_results,file="2_annotation/subcluster/CD8PlusDendritic_cell_gene_erichment_results.RData")

write.csv(gene_erichment_results[["cDC"]][["topGO"]],"2_annotation/subcluster/CD8PlusDendritic_cellsub.cDC.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["pDC"]][["topGO"]],"2_annotation/subcluster/CD8PlusDendritic_cellsub.pDC.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["mDC"]][["topGO"]],"2_annotation/subcluster/CD8PlusDendritic_cellsub.mDC.GO.csv",quote=F,row.names=F)

library(enrichplot)
#CD8PlusDendritic_cell_0
gene_erichment_results[["Dendritic_cell"]][["topGO"]][1:5,]
pdf("2_annotation/subcluster/CD8PlusDendritic_cell_Dendritic_cell.pdf",width=8,height=10)
dotplot(gene_erichment_results[["Dendritic_cell"]][["DGN"]], showCategory=30) 
dev.off()
#CD8PlusDendritic_cell_1
gene_erichment_results[["Monocyte_derived_dendritic_cell"]][["topGO"]][1:5,]
pdf("2_annotation/subcluster/CD8PlusDendritic_cell_Monocyte_derived_dendritic_cell.pdf",width=8,height=10)
dotplot(gene_erichment_results[["Monocyte_derived_dendritic_cell"]][["DGN"]], showCategory=30) 
dev.off()
#CD8PlusDendritic_cell_2
gene_erichment_results[["NPC"]][["topGO"]][1:5,]
pdf("2_annotation/subcluster/CD8PlusDendritic_cell_NPC.pdf",width=8,height=10)
dotplot(gene_erichment_results[["NPC"]][["DGN"]], showCategory=30) 
dev.off()
#CD8PlusDendritic_cell_3
gene_erichment_results[["Lymphoid_dendritic_cells"]][["topGO"]][1:5,]
pdf("2_annotation/subcluster/CD8PlusDendritic_cell_Lymphoid_dendritic_cells.pdf",width=8,height=10)
dotplot(gene_erichment_results[["Lymphoid_dendritic_cells"]][["DGN"]], showCategory=30) 
dev.off()




load("2_annotation/subcluster/scRNAsub.CD8PlusDendritic_cell.RData")
new.cluster.ids<-c("cDC", "pDC", "mDC", "NA","NA")
names(new.cluster.ids) <- levels(scRNAsub.CD8PlusDendritic_cell)
scRNAsub.CD8PlusDendritic_cell <- RenameIdents(scRNAsub.CD8PlusDendritic_cell, new.cluster.ids)
scRNAsub.CD8PlusDendritic_cell$celltype<-Idents(scRNAsub.CD8PlusDendritic_cell )
scRNAsub.CD8PlusDendritic_cell_rmNA <- subset(scRNAsub.CD8PlusDendritic_cell, subset = (celltype == "cDC" | celltype == "pDC" | celltype == "mDC"))
scRNAsub.CD8PlusDendritic_cell<-scRNAsub.CD8PlusDendritic_cell_rmNA
scRNAsub.CD8PlusDendritic_cell <- RunUMAP(scRNAsub.CD8PlusDendritic_cell, dims = 1:18)
scRNAsub.CD8PlusDendritic_cell$celltype<-Idents(scRNAsub.CD8PlusDendritic_cell)
save(scRNAsub.CD8PlusDendritic_cell,file="2_annotation/subcluster/scRNAsub.CD8PlusDendritic_cell.afteranno.RData")
scRNAsub.CD8PlusDendritic_cell.markers <- FindAllMarkers(object = scRNAsub.CD8PlusDendritic_cell, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.CD8PlusDendritic_cell.markers,"2_annotation/subcluster/scRNAsub.CD8PlusDendritic_cell.markers.anno.csv",sep=",",quote=F)
save(scRNAsub.CD8PlusDendritic_cell.markers,file="2_annotation/subcluster/scRNAsub.CD8PlusDendritic_cell.markers.afteranno.RData")
pdf("2_annotation/subcluster/tsne.CD8PlusDendritic_cell.integrate.pdf", width = 10, height = 10)
DimPlot(scRNAsub.CD8PlusDendritic_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
pdf("2_annotation/subcluster/tsne.CD8PlusDendritic_cell.integrate.nolegend.pdf", width = 10, height = 10)
DimPlot(scRNAsub.CD8PlusDendritic_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8, split.by = 'label', group.by = 'celltype')+NoLegend()
dev.off()
top30<-scRNAsub.CD8PlusDendritic_cell.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_annotation/subcluster/scRNAsub.CD8PlusDendritic_cell.top30.anno.csv",sep=",",quote=F)

features.plot <- c("HIST1H4C","STMN1","HMGB2","IGLC2","GAS5","IL32","CST3","HLA-DPA1","HLA-DQB1")
pdf("2_annotation/subcluster/CD8PlusDendritic_cell.dittoDotPlot.pdf",width = 10, height = 5)
DotPlot(object = scRNAsub.CD8PlusDendritic_cell, features = features.plot,dot.scale = 12,cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()


Meningeal_cell.cells <- subset(WBY.combined, idents = "Meningeal_cell")
Idents(Meningeal_cell.cells) <- "label"
avg.Meningeal_cell.cells <- as.data.frame(log1p(AverageExpression(Meningeal_cell.cells, verbose = FALSE)$RNA))
avg.Meningeal_cell.cells$gene <- rownames(avg.Meningeal_cell.cells)
write.table(avg.Meningeal_cell.cells,"2_annotation/DEGs/avg.Meningeal_cell.cells.csv",sep=",",quote=F)
WBY.combined$celltype.label <- paste(Idents(WBY.combined), WBY.combined$label, sep = "_")
WBY.combined$celltype <- Idents(WBY.combined)
Idents(WBY.combined) <- "celltype.label"
Meningeal_cell.response <- FindMarkers(WBY.combined, ident.1 = "Meningeal_cell_cancer", ident.2 = "Meningeal_cell_normal", verbose = FALSE)
head(Meningeal_cell.response, n = 15)
write.table(Meningeal_cell.response,"2_annotation/DEGs/Meningeal_cell.degs.csv",sep=",")
degs_meningeal<-Meningeal_cell.response
degs_filtered<-degs_meningeal[(degs_meningeal$avg_log2FC>log2(1.5)|degs_meningeal$avg_log2FC<=-log2(1.5) ) & (degs_meningeal$p_val_adj<0.05),]
features_degs<-degs_filtered[order(degs_filtered$avg_log2FC),]
pdf("2_annotation/DEGs/meningeal.heatmap.pdf",width = 20, height = 18)
DoHeatmap(Meningeal_cell.cells, features = row.names(features_degs)) + NoLegend()
dev.off()
degs_meningeal<-Meningeal_cell.response
degs_filtered<-degs_meningeal[(degs_meningeal$avg_log2FC>log2(1.5)|degs_meningeal$avg_log2FC<=-log2(1.5) ) & (degs_meningeal$p_val_adj<0.05),]
features_degs<-degs_filtered[order(degs_filtered$avg_log2FC),]
pdf("2_annotation/DEGs/meningeal.heatmap.label.pdf",width = 20, height = 18)
DoHeatmap(Meningeal_cell.cells, features = row.names(features_degs))
dev.off()

pdf("2_annotation/DEGs/Meningeal_cell.figure1.pdf",width = 20, height = 18)
FeaturePlot(WBY.combined, features = c("IL6ST","IRF1","PRDM1","IGFBP1","IGFBP2"), reduction = "tsne", split.by = "label", max.cutoff = 3, 
            cols = c("grey", "red"))
dev.off()
pdf("2_annotation/DEGs/Meningeal_cell.figure2.pdf",width = 10, height = 18)
plots <- VlnPlot(WBY.combined, features = c("CARD16","IRF1","PRDM1","IGFBP2","IGFBP7"), split.by = "label", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
dev.off()
pdf("2_annotation/DEGs/Meningeal_cell.figure3.pdf")
genes.to.label = c("CARD16","IRF1","PRDM1","IGFBP2","IGFBP7")
p1 <- ggplot(avg.Meningeal_cell.cells, aes(cancer, normal)) + geom_point() + ggtitle("Meningeal_cell")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1
dev.off()

load("2_annotation/subcluster/scRNAsub.tme.afteranno.RData")
features.plot <- c("NCL","NEAT1","LTB","TCL1A","CD79B","CD79A","MT-ND3","CST7","CXCL13","IGHG3","IGHG4","IGHG1")
pdf("2_annotation/subcluster/Bcellsub.dittoDotPlot.pdf",width = 10, height = 5)
DotPlot(object = scRNAsub.tme, features=features.plot,dot.scale = 12,cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()


load("2_annotation/subcluster/scRNAsub.tme.afteranno.RData")
features.plot <- c("B2M","CCL5","TUBB","LDHA","HIST1H4C","ARL6IP1","CD69","TSC22D3","CD44","UCP2")
pdf("2_annotation/subcluster/Bcellsub.dittoDotPlot.pdf",width = 10, height = 5)
DotPlot(object = scRNAsub.tme, features=features.plot,dot.scale = 12,cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()


load("2_annotation/subcluster/scRNAsub.MPC.afteranno.RData")
features.plot <- c("LTB","SAT1","CD79A","JUNB","CD69","HIST1H4C","UBE2C","TUBB4B","CKS1B","TYMS")
pdf("2_annotation/subcluster/MPCsub.dittoDotPlot.pdf",width = 10, height = 3)
DotPlot(object = scRNAsub.MPC, features=features.plot,dot.scale = 12,cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()


##############################################################################################################################





# 加载需要的R包
library(Seurat)
library(monocle3)
data <- as(as.matrix(scRNAsub.tme@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub.tme@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
monocle_cds@phenoData@data[["cell_type"]] <- as.character(scRNAsub.tme@active.ident)
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
#Filtering low-quality cells
monocle_cds <- detectGenes(monocle_cds, min_expr = 3 )
disp_table <- dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
monocle_cds <- reduceDimension(
 monocle_cds,
 max_components = 2,
 method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
pdf("plot_cell_trajectory.pdf", width = 24, height = 18)
plot_cell_trajectory(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
dev.off()
pdf("plot_cell_trajectory.pdf", width = 24, height = 18)
plot_cell_trajectory(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "cell_type")
dev.off()

