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
# load the scRNA-seq data
# create each individual Seurat object for every sample
# create Seurat object cancer
seurat_data<-Read10X(data.dir=paste0("0_counts/","cancer"))
#seurat_obj<-CreateSeuratObject(counts=seurat_data,project="cancer")
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="cancer")
assign("cancer",seurat_obj)
cancer<-NormalizeData(cancer)
write.table(GetAssayData(cancer),"exp_cancer.csv",sep=",",quote=F)
pdf("1_preoperation/figures/pre-operation/cancer.counts.vs.features.pdf")
plot(x=cancer@meta.data$nCount_RNA,y=cancer@meta.data$nFeature_RNA)
dev.off()
# check the metadata in the new Seurat objects
head(cancer@meta.data)
tail(cancer@meta.data)
# Create .RData object to load at any time
save(cancer, file="1_preoperation/data/cancer.combined.RData")
cancer$log10GenesPerUMI <- log10(cancer$nFeature_RNA) / log10(cancer$nCount_RNA)
cancer$mitoRatio <- PercentageFeatureSet(object = cancer, pattern = "^MT-")
cancer$mitoRatio <- cancer@meta.data$mitoRatio / 100
cancermetadata <- cancer@meta.data
cancermetadata$cells <- rownames(cancermetadata)
cancermetadata <- cancermetadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
cancer
cancer@meta.data <- cancermetadata
counts <- GetAssayData(object = cancer, slot = "counts")
cancer <- CreateSeuratObject(counts, meta.data = cancer@meta.data)
cancer$label <- "cancer"
cancer_norm <- NormalizeData(cancer, normalization.method = "LogNormalize", scale.factor = 10000)
cancer_norm <- FindVariableFeatures(cancer_norm, selection.method = "vst", nfeatures = 3000)
pdf("1_preoperation/figures/pre-operation/cancer_Visualize_QC.pdf", width = 12, height = 6)
VlnPlot(cancer, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(cancer, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(cancer, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
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
filtered_cancer <- subset(x = cancer, 
                          subset= (nUMI < 20000) & 
                            (nGene > 100) &
                            (nGene < 3500) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.40))
counts <- GetAssayData(object = filtered_cancer, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_cancer <- CreateSeuratObject(filtered_counts, meta.data = cancer@meta.data)
filtered_cancer$label <- "cancer"
save(filtered_cancer, file="1_preoperation/data/filtered_cancer.RData")
write.table(GetAssayData(filtered_cancer),"exp_filtered_cancer.csv",sep=",",quote=F)
filtered_cancer_norm<-NormalizeData(filtered_cancer)
# create Seurat object normal
setwd("/ifs1/Grp8/liuzhe/scRNA/")
seurat_data<-Read10X(data.dir=paste0("0_counts/normal/","HFA567_total.filtered_gene_matrices"))
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="normal")
assign("normal1",seurat_obj)
normal1<-NormalizeData(normal1)
write.table(GetAssayData(normal1),"exp_normal1.csv",sep=",",quote=F)
seurat_data<-Read10X(data.dir=paste0("0_counts/normal/","HFA570_total.filtered_gene_matrices"))
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="normal")
assign("normal2",seurat_obj)
normal2<-NormalizeData(normal2)
write.table(GetAssayData(normal2),"exp_normal2.csv",sep=",",quote=F)
seurat_data<-Read10X(data.dir=paste0("0_counts/normal/","HFA571_total.filtered_gene_matrices"))
seurat_obj<-CreateSeuratObject(counts=seurat_data,project="normal")
assign("normal3",seurat_obj)
normal3<-NormalizeData(normal3)
write.table(GetAssayData(normal3),"exp_normal3.csv",sep=",",quote=F)
normal.normalized.combined <- merge(normal1, y = c(normal2, normal3), add.cell.ids = c("N1", "N2", "N3"), project = "normal", merge.data = TRUE)
normal<-normal.normalized.combined
write.table(GetAssayData(normal),"exp_normal.csv",sep=",",quote=F)
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
                          subset= (nUMI < 25000) & 
                            (nGene > 100) &
                            (nGene < 6000) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))
counts <- GetAssayData(object = filtered_normal, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 20
filtered_counts <- counts[keep_genes, ]
filtered_normal <- CreateSeuratObject(filtered_counts, meta.data = normal@meta.data)
filtered_normal$label <- "normal"
save(filtered_normal, file="1_preoperation/data/filtered_normal.RData")
write.table(GetAssayData(filtered_normal),"exp_filtered_normal.csv",sep=",",quote=F)
filtered_normal_norm<-NormalizeData(filtered_normal)
# integrate normal and cancer sample
WBY.anchors <- FindIntegrationAnchors(object.list = list(filtered_cancer_norm, filtered_normal_norm), dims = 1:30)
save(WBY.anchors, file="1_preoperation/data/integrated.anchors_seurat20210719.RData")
WBY.combined <- IntegrateData(anchorset = WBY.anchors, dims = 1:30)
WBY.combined <- FindVariableFeatures(WBY.combined, selection.method = "vst", nfeatures = 3000)
save(WBY.combined, file="1_preoperation/data/integrated.combined_seurat20210510.RData")
write.table(GetAssayData(WBY.combined),"exp_WBY.combined.csv",sep=",",quote=F)
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
# calculate the percentage of each auto annotation cluster
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
new.cluster.ids<-c("T_cell", "Multilymphoid_progenitor_cell", "T_cell", "Microglia", "B_cell", "Macrophage", "Astrocyte", "Neuron", "Interneuron", "T_cell", "CD8+Dendritic_cell", "Intermediate_progenitor_cell", "Astrocyte", "Oligodendrocyte", "Meningeal_cell", "Neural_progenitor_cell", "Oligodendrocyte_progenitor_cell")
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

delta.genes <- c("CCL5","RGS1","TUBA1B","HMGB2","CADM2","CSRP2","SLC25A5","CD79A","LYZ","NPC2","MSMO1","FGFBP3","MAP1B","NFIB","PLS3","SLAIN1","CST3","HLA-DPA1","EOMES","BAALC","PLP1","PPP1R14A","MGP","SPARCL1","ASPM","CENPF","OLIG2","OLIG1")
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
features.plot <- c("CCL5","RGS1","TUBA1B","HMGB2","CADM2","CSRP2","SLC25A5","CD79A","LYZ","NPC2","MSMO1",
                   "FGFBP3","MAP1B","NFIB","PLS3","SLAIN1","CST3","HLA-DPA1","EOMES","BAALC","PLP1","PPP1R14A",
                   "MGP","SPARCL1","ASPM","CENPF","OLIG2","OLIG1")
pdf("2_annotation/anno/markergenes.dotplot.pdf",width = 10, height = 8)
DotPlot(object = WBY.combined, features = features.plot,  cols = c("lightgrey", "red"))
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
# B-cell recluster
B_cell.subset <- subset(WBY_cancer@meta.data, celltype=="B_cell")
scRNAsub.B_cell <- subset(WBY_cancer, cells=row.names(B_cell.subset))
scRNAsub.B_cell <- FindVariableFeatures(scRNAsub.B_cell, selection.method = "vst", nfeatures = 2000)
scale.genes.tme <-  rownames(scRNAsub.B_cell)
scRNAsub.B_cell <- ScaleData(scRNAsub.B_cell, features = scale.genes.tme)
scRNAsub.B_cell <- RunPCA(scRNAsub.B_cell, features = VariableFeatures(scRNAsub.B_cell))
pdf("2_annotation/subcluster/Determine.tme.pcnumber.pdf")
ElbowPlot(scRNAsub.B_cell, ndims=20, reduction="pca")
dev.off()
pc.num=1:10
scRNAsub.B_cell <- FindNeighbors(scRNAsub.B_cell, dims = pc.num) 
scRNAsub.B_cell <- FindClusters(scRNAsub.B_cell, resolution = 0.6)
table(scRNAsub.B_cell@meta.data$seurat_clusters)
metadata <- scRNAsub.B_cell@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'2_annotation/subcluster/tme.cell_cluster.csv',row.names = F)
#tSNE
scRNAsub.B_cell = RunTSNE(scRNAsub.B_cell, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub.B_cell, 'tsne')
write.csv(embed_tsne,'2_annotation/subcluster/tme.embed_tsne.csv')
pdf("2_annotation/subcluster/tsne_Bcell.pdf")
DimPlot(scRNAsub.B_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8)
dev.off()
diff.wilcox = FindAllMarkers(scRNAsub.B_cell)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "2_annotation/subcluster/tme.diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "2_annotation/subcluster/tme.top30_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub.B_cell,file="2_annotation/subcluster/scRNAsub.B_cell.RData")
new.cluster.ids<-c("B_cell-1", "Pre-B_cell", "MKI67+_progenitor_cell", "B_cell-2", "B_cell-1", "B_cell-3")
names(new.cluster.ids) <- levels(scRNAsub.B_cell)
scRNAsub.B_cell <- RenameIdents(scRNAsub.B_cell, new.cluster.ids)
scRNAsub.B_cell <- RunUMAP(scRNAsub.B_cell, dims = 1:10)
scRNAsub.B_cell$celltype<-Idents(scRNAsub.B_cell)
save(scRNAsub.B_cell,file="2_annotation/subcluster/scRNAsub.B_cell.afteranno.RData")
scRNAsub.B_cell.markers <- FindAllMarkers(object = scRNAsub.B_cell, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.B_cell.markers,"2_annotation/subcluster/scRNAsub.B_cell.markers.anno.csv",sep=",",quote=F)
save(scRNAsub.B_cell.markers,file="2_annotation/subcluster/scRNAsub.B_cell.markers.afteranno.RData")
pdf("2_annotation/subcluster/tsne.bcell.integrate.pdf", width = 10, height = 10)
DimPlot(scRNAsub.B_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
top30<-scRNAsub.B_cell.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_annotation/subcluster/scRNAsub.B_cell.top30.anno.csv",sep=",",quote=F)
require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(scRNAsub.B_cell))
names(geneList)<-row.names(scRNAsub.B_cell)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(scRNAsub.B_cell.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(scRNAsub.B_cell.markers,cluster==c1)$gene
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
write.csv(gene_erichment_results[["Pre-B_cell"]][["topGO"]],"2_annotation/subcluster/Bcellsub.Pre-B_cell.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["MKI67+_progenitor_cell"]][["topGO"]],"2_annotation/subcluster/Bcellsub.MKI67+_progenitor_cell.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["B_cell-2"]][["topGO"]],"2_annotation/subcluster/Bcellsub.B_cell-2.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["B_cell-3"]][["topGO"]],"2_annotation/subcluster/Bcellsub.B_cell-3.GO.csv",quote=F,row.names=F)

# MPC recluster
WBY.combined$celltype<-Idents(WBY.combined)
WBY_cancer<-subset(x=WBY.combined,subset = label == "cancer")
MPC.subset <- subset(WBY_cancer@meta.data, celltype=="Multilymphoid_progenitor_cell")
scRNAsub.MPC <- subset(WBY_cancer, cells=row.names(MPC.subset))
scRNAsub.MPC <- FindVariableFeatures(scRNAsub.MPC, selection.method = "vst", nfeatures = 2000)
scale.genes.MPC <-  rownames(scRNAsub.MPC)
scRNAsub.MPC <- ScaleData(scRNAsub.MPC, features = scale.genes.MPC)
scRNAsub.MPC <- RunPCA(scRNAsub.MPC, features = VariableFeatures(scRNAsub.MPC))
pdf("2_annotation/subcluster/Determine.MPC.pcnumber.pdf")
ElbowPlot(scRNAsub.MPC, ndims=20, reduction="pca")
dev.off()
pc.num=1:12
scRNAsub.MPC <- FindNeighbors(scRNAsub.MPC, dims = pc.num) 
scRNAsub.MPC <- FindClusters(scRNAsub.MPC, resolution = 1.2)
table(scRNAsub.MPC@meta.data$seurat_clusters)
metadata <- scRNAsub.MPC@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'2_annotation/subcluster/MPC.cell_cluster.csv',row.names = F)
#tSNE
scRNAsub.MPC = RunTSNE(scRNAsub.MPC, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub.MPC, 'tsne')
write.csv(embed_tsne,'2_annotation/subcluster/MPC.embed_tsne.csv')
pdf("2_annotation/subcluster/tsne_MPC.pdf")
DimPlot(scRNAsub.MPC, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8)
dev.off()
diff.wilcox = FindAllMarkers(scRNAsub.MPC)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "2_annotation/subcluster/MPC.diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "2_annotation/subcluster/MPC.top30_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub.MPC,file="2_annotation/subcluster/scRNAsub.MPC.RData")
new.cluster.ids<-c("MPC-1", "MPC-1", "MPC-1", "MPC-2", "MPC-1", "MPC-2", "MPC-2", "MPC-2", "MPC-1", "MPC-2")
names(new.cluster.ids) <- levels(scRNAsub.MPC)
scRNAsub.MPC <- RenameIdents(scRNAsub.MPC, new.cluster.ids)
scRNAsub.MPC <- RunUMAP(scRNAsub.MPC, dims = 1:12)
scRNAsub.MPC$celltype<-Idents(scRNAsub.MPC)
save(scRNAsub.MPC,file="2_annotation/subcluster/scRNAsub.MPC.afteranno.RData")
scRNAsub.MPC.markers <- FindAllMarkers(object = scRNAsub.MPC, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(scRNAsub.MPC.markers,"2_annotation/subcluster/scRNAsub.MPC.markers.anno.csv",sep=",",quote=F)
save(scRNAsub.MPC.markers,file="2_annotation/subcluster/scRNAsub.MPC.markers.afteranno.RData")
pdf("2_annotation/subcluster/tsne.MPC.integrate.pdf", width = 10, height = 10)
DimPlot(scRNAsub.MPC, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8, split.by = 'label', group.by = 'celltype')
dev.off()
top30<-scRNAsub.MPC.markers %>% group_by(cluster) %>% top_n(n=30,wt=avg_log2FC)
write.table(top30,"2_annotation/subcluster/scRNAsub.MPC.top30.anno.csv",sep=",",quote=F)

require(org.Hs.eg.db)
library(topGO)
library(DOSE)
#devtools::install_github("eliocamp/ggnewscale")
library("ggnewscale")
x=as.list(org.Hs.egALIAS2EG)
geneList<-rep(0,nrow(scRNAsub.MPC))
names(geneList)<-row.names(scRNAsub.MPC)
geneList<-geneList[intersect(names(geneList),names(x))]
newwallgenes=names(geneList)
for (ii in 1:length(geneList)){
  names(geneList)[ii]<-x[[names(geneList)[ii]]][1]
  
}
gene_erichment_results=list()
for (c1 in as.character(unique(levels(scRNAsub.MPC.markers$cluster)))){
  print(paste0("RUN ", c1))
  testgenes<-subset(scRNAsub.MPC.markers,cluster==c1)$gene
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
save(gene_erichment_results,file="2_annotation/anno/MPCsub.gene_erichment_results.RData")
write.csv(gene_erichment_results[["MPC-1"]][["topGO"]],"2_annotation/subcluster/MPCsub.MPC-1.GO.csv",quote=F,row.names=F)
write.csv(gene_erichment_results[["MPC-2"]][["topGO"]],"2_annotation/subcluster/MPCsub.MPC-2.GO.csv",quote=F,row.names=F)

# marker gene distribution
pdf("2_annotation/anno/markergenes.MPC.pdf")
FeaturePlot(WBY.combined, features = "TUBA1B", blend.threshold = 1,reduction="tsne")
dev.off()
pdf("2_annotation/anno/markergenes.Bcell.pdf")
FeaturePlot(WBY.combined, features = "CD79A", reduction="tsne")
dev.off()
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


###DEGs analysis for macrophage and Meningeal_cell between normal and cancer 
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
macrophage.response <- FindMarkers(WBY.combined, ident.1 = "Macrophage_cancer", ident.2 = "Macrophage_normal", verbose = FALSE)
head(macrophage.response, n = 15)
write.table(macrophage.response,"2_annotation/DEGs/macrophage.degs.csv",sep=",")
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
WBY.combined$celltype<-Idents(WBY.combined)
pdf("2_annotation/anno/markergenes.macrophage.vioplot.pdf",width=15,height=10)
VlnPlot(WBY.combined, features = c("CST3","NPC2"), group.by = 'celltype', pt.size = 0) 
dev.off()
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
### T cell recluster
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
pc.num=1:8
scRNAsub.T_cell <- FindNeighbors(scRNAsub.T_cell, dims = pc.num) 
scRNAsub.T_cell <- FindClusters(scRNAsub.T_cell, resolution = 0.7)
table(scRNAsub.T_cell@meta.data$seurat_clusters)
metadata <- scRNAsub.T_cell@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'2_annotation/subcluster/T_cell.cell_cluster.csv',row.names = F)
#tSNE
scRNAsub.T_cell = RunTSNE(scRNAsub.T_cell, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub.T_cell, 'tsne')
write.csv(embed_tsne,'2_annotation/subcluster/T_cell.embed_tsne.csv')
pdf("2_annotation/subcluster/tsne_T_cell.pdf")
DimPlot(scRNAsub.T_cell, reduction = "tsne", label = TRUE, pt.size=1.5,label.size = 8)
dev.off()
diff.wilcox = FindAllMarkers(scRNAsub.T_cell)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(all.markers, "2_annotation/subcluster/T_cell.diff_genes_wilcox.csv", row.names = F)
write.csv(top30, "2_annotation/subcluster/T_cell.top30_diff_genes_wilcox.csv", row.names = F)
save(scRNAsub.T_cell,file="2_annotation/subcluster/scRNAsub.T_cell.RData")
load("2_annotation/subcluster/scRNAsub.T_cell.RData")
new.cluster.ids<-c("NKT", "NA", "T_cell", "T_cell-1", "NA","T_cell", "T_cell", "T_helper_cell")
names(new.cluster.ids) <- levels(scRNAsub.T_cell)
scRNAsub.T_cell <- RenameIdents(scRNAsub.T_cell, new.cluster.ids)
scRNAsub.T_cell$celltype<-Idents(scRNAsub.T_cell )
scRNAsub.T_cell_rmNA <- subset(scRNAsub.T_cell, subset = (celltype == "NKT" | celltype == "T_cell" | celltype == "T_cell-1" | celltype == "T_helper_cell"))
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
features.plot <- c("NCL","HMGB2","RPLP0","RPS5","MT-ND5","LYZ","IL7R","KLF6")
pdf("2_annotation/anno/Tcellsub.dittoDotPlot.pdf",width = 10, height = 8)
DotPlot(object = scRNAsub.T_cell, features = features.plot,  cols = c("lightgrey", "red"))
dev.off()
## CD8+Dendritic_cell recluster
WBY.combined$celltype<-Idents(WBY.combined)
WBY_cancer<-subset(x=WBY.combined,subset = label == "cancer")
CD8PlusDendritic_cell.subset <- subset(WBY_cancer@meta.data, celltype=="CD8+Dendritic_cell")
scRNAsub.CD8PlusDendritic_cell <- subset(WBY_cancer, cells=row.names(CD8PlusDendritic_cell.subset))
scRNAsub.CD8PlusDendritic_cell <- FindVariableFeatures(scRNAsub.CD8PlusDendritic_cell, selection.method = "vst", nfeatures = 2000)
scale.genes.CD8PlusDendritic_cell <-  rownames(scRNAsub.CD8PlusDendritic_cell)
scRNAsub.CD8PlusDendritic_cell <- ScaleData(scRNAsub.CD8PlusDendritic_cell, features = scale.genes.CD8PlusDendritic_cell)
scRNAsub.CD8PlusDendritic_cell <- RunPCA(scRNAsub.CD8PlusDendritic_cell, features = VariableFeatures(scRNAsub.CD8PlusDendritic_cell))
pdf("2_annotation/subcluster/Determine.CD8PlusDendritic_cell.pcnumber.pdf")
ElbowPlot(scRNAsub.CD8PlusDendritic_cell, ndims=20, reduction="pca")
dev.off()
pc.num=1:13
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
load("2_annotation/subcluster/scRNAsub.CD8PlusDendritic_cell.RData")
new.cluster.ids<-c("Dendritic_cell", "Monocyte_derived_dendritic_cell", "Monocyte_derived_dendritic_cell", "NPC", "Lymphoid_dendritic_cells","Dendritic_cell", "Dendritic_cell")
names(new.cluster.ids) <- levels(scRNAsub.CD8PlusDendritic_cell)
scRNAsub.CD8PlusDendritic_cell <- RenameIdents(scRNAsub.CD8PlusDendritic_cell, new.cluster.ids)
scRNAsub.CD8PlusDendritic_cell$celltype<-Idents(scRNAsub.CD8PlusDendritic_cell )
scRNAsub.CD8PlusDendritic_cell_rmNA <- subset(scRNAsub.CD8PlusDendritic_cell, subset = (celltype == "Dendritic_cell" | celltype == "Monocyte_derived_dendritic_cell" | celltype == "Lymphoid_dendritic_cells"))
scRNAsub.CD8PlusDendritic_cell<-scRNAsub.CD8PlusDendritic_cell_rmNA
scRNAsub.CD8PlusDendritic_cell <- RunUMAP(scRNAsub.CD8PlusDendritic_cell, dims = 1:13)
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
features.plot <- c("HIST1H4C","TUBB","RPL22L1","RGCC","FOS","IRF8","RGS1","CCL4L2","CCL5")
pdf("2_annotation/anno/CD8PlusDendritic_cell.dittoDotPlot.pdf",width = 10, height = 8)
DotPlot(object = scRNAsub.CD8PlusDendritic_cell, features = features.plot,  cols = c("lightgrey", "red"))
dev.off()
# Meningeal_cell recluster
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