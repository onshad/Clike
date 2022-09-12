### Seurat V3
library(BoutrosLab.plotting.general)
library(dplyr)
library(monocle)
packageVersion("monocle")
library(canprot)
source("../functions.seurat.R")
library(Seurat)
packageVersion("Seurat")
# [1] ‘4.0.1’
# BiocManager::available()

setwd("/ ... ")

# saveRDS(TKO, "TKO.moreCells.rds")
# TKO=readRDS("TKO.moreCells.rds")

# saveRDS(TKOe, "TKOe.moreCells.rds")
# TKOe=readRDS("TKOe.moreCells.rds")

# saveRDS(TKOe, "TKOe.rds")
# TKOe=readRDS("TKOe.rds")



### read samples
### process from raw counts
# data.directory="/ ... /counts/mouse"
data.directory="/ ... /counts/mouse"
matrix.directory=""
dir.raw=list.dirs(path = data.directory, full.names = TRUE, recursive = F)
dir.raw=paste(dir.raw, matrix.directory, sep="")
# # the sample names' starting position 
pos.1=regexpr(data.directory, dir.raw, perl=TRUE)
pos.1=attr(pos.1, "match.length")+2
# # the sample names' ending position
# pos.2=regexpr(matrix.directory, dir.raw, perl=TRUE)-1
# pos.2=nchar(dir.raw)
pos.2=regexpr(".matrix", dir.raw, perl=TRUE)
sample.names=substring(dir.raw, pos.1, pos.2-1)
print(sample.names)
# [1] "SKO"   "SKOCa" "TKO"
# "pten", "pten.castration", "pten.p53.rb"
inte.list=list()# the list prepared for further dataset integration
for(i in 1:length(dir.raw))
{
   temp=Read10X(data.dir = dir.raw[i])
   inte.list[i]=CreateSeuratObject(counts = temp, project = sample.names[i], min.cells = 0, min.features = 100)
}
for (i in 1:length(inte.list))
{
  inte.list[[i]][["percent.mt"]] <- PercentageFeatureSet(inte.list[[i]], pattern = "^mt-")
}
# inte.list.bak=inte.list
# inte.list=inte.list.bak
for(i in 1:length(inte.list))
{
  png(paste("temp.featureScatter_",sample.names[i],"_quality.png", sep=""), width = 15*150, height = 4*150 )
  plot1 <- FeatureScatter(inte.list[[i]], feature1 = "nFeature_RNA", feature2 = "percent.mt")
  plot1=plot1+geom_vline(xintercept=500)+geom_vline(xintercept=1000)+geom_hline(yintercept=20)
  plot2 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2=plot2+geom_vline(xintercept=3000)+geom_vline(xintercept=1500)+geom_hline(yintercept=20)
  plot3 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2 + plot3)
  dev.off()
}
for(i in 1:length(inte.list))
{
  # temp<- subset(inte.list[[i]], subset = nFeature_RNA > 1000 & nCount_RNA > 6000 & percent.mt < 3 ) # & top6k==TRUE)  previous version
  temp<- subset(inte.list[[i]], subset = nFeature_RNA > 500 & nCount_RNA > 1500 & percent.mt < 20 ) # & top6k==TRUE)
  print(sample.names[i])
  print(dim(temp))
}
### nFeature_RNA>500 and percent.mt<60 and NOT - restricted on top 6000 cells with highest nFeature_RNA
for(i in 1:length(inte.list))
{
  ### NOT - only consider cells from top 6000 cells
  # inte.list[[i]]@meta.data$top6k=rep(TRUE, ncol(inte.list[[i]]))
  # temp=sort(inte.list[[i]]@meta.data$nFeature_RNA, decreasing=TRUE)[6000]
  # inte.list[[i]]@meta.data$top6k[inte.list[[i]]@meta.data$nFeature_RNA<=temp]=FALSE
  # print(table(inte.list[[i]]@meta.data$top6k))
  # inte.list[[i]]<- subset(inte.list[[i]], subset = nFeature_RNA > 500 & percent.mt < 60 ) # & top6k==TRUE) # --->>> this is last version
  inte.list[[i]]<- subset(inte.list[[i]], subset = nFeature_RNA > 500 & nCount_RNA > 1500 & percent.mt < 20 ) # & top6k==TRUE)
  # inte.list[[i]]@meta.data$top6k=NULL
}
# # filter out cells with > 2 x median nFeature_RNA
# for(i in 1:length(inte.list))
# {
#   inte.list[[i]]@meta.data$top6k=NULL
#   temp=median(inte.list[[i]]@meta.data$nFeature_RNA)*2
#   inte.list[[i]]@meta.data$double=rep(FALSE, ncol(inte.list[[i]]))
#   inte.list[[i]]@meta.data$double[inte.list[[i]]@meta.data$nFeature_RNA>temp]=TRUE
#   print(sample.names[i])
#   print(table(inte.list[[i]]@meta.data$double))
#   print(temp)
#   ### NO need to, the nFeature_RNA~nCount_RNA plot show continuous distribution
#   inte.list[[i]]<- subset(inte.list[[i]], subset = (top6k==FALSE) )
#   inte.list[[i]]@meta.data$double=NULL
# }
# NOT TO individually normalize each batch
# for(i in 1:length(inte.list))
# {
#   inte.list[[i]]=SCTransform(inte.list[[i]], vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE)
# }
# TKO=merge(x = inte.list[[1]], y = inte.list[2:length(inte.list)] )
TKO=inte.list[[3]]
dim(TKO)
# [1] 31253  9624
TKO@meta.data$fig.sample=as.character(TKO@meta.data$orig.ident)
# TKO@meta.data$fig.sample=factor(TKO@meta.data$fig.sample, levels=c("control.sample1", "control.sample2", "smoke1m.sample1")
# pat2pheno=data.frame(sample=sort(unique(as.character(TKO@meta.data$fig.sample))))

# pat2pheno$sample
# #  [1] "control.sample1"    "control.sample2"    "offspring.sample1"
# #  [4] "offspring.sample2"  "offspring.sample3"  "offspring2.sample1"
# #  [7] "quit1p5m.sample1"   "quit3m.sample1"     "smoke1m.sample1"
# # [10] "smoke1m.sample2"    "smoke3m.sample1"    "smoke3m.sample2"
# # [13] "smoke3m.sample3"    "smoke6m.sample1"    "smoke6m.sample2"
# pat2pheno$pheno=c("control", "control", "offspring", "offspring", "offspring", "offspring2", "quit1p5m", "quit3m", "smoke1m", "smoke1m", "smoke3m", "smoke3m", "smoke3m", "smoke6m", "smoke6m")
# pat2pheno

# TKO@meta.data$pheno=pat2pheno$pheno[match(as.character(TKO@meta.data$fig.sample), pat2pheno$sample)]
# # TKO@meta.data$pheno=as.character(TKO@meta.data$pheno)
# table(TKO@meta.data$pheno)
# # check ingredients in meta.data
# str(TKO@meta.data)

TKO<- SCTransform(TKO, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes=FALSE, method = "glmGamPoi")
DefaultAssay(TKO) <- "SCT"
str(TKO@assays$SCT@var.features)
# chr [1:3000] "Cxcl2" "Col1a1" "Ccl3" "Il1b" "Col1a2" "Sparc" "S100a9" ...

TKO <- RunPCA(TKO, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(TKO,  ndims = 100)
dev.off()
TKO <- RunUMAP(TKO, dims = 1:25, verbose = FALSE)
TKO <- FindNeighbors(TKO, dims = 1:25, verbose = FALSE)
TKO <- FindClusters(TKO, verbose = FALSE, resolution = 0.2)
# DefaultAssay(TKO)="RNA"
# TKO <- NormalizeData(TKO, normalization.method = "LogNormalize", scale.factor = 10000)
# TKO <- ScaleData(TKO, features = rownames(TKO), assay="RNA" )
# DefaultAssay(TKO)="SCT"

pdf("dimplot.pdf")
DimPlot(TKO, label = TRUE) # + NoLegend()
dev.off()

TKO@meta.data$fig.cluster=TKO@meta.data$seurat_clusters
table(TKO@meta.data$fig.cluster)
#    0    1    2    3    4    5    6    7    8    9   10   11
# 1843 1782 1582 1114 1114  560  489  453  302  174  133   78

TKO <- SetIdent(TKO, value = "fig.cluster")
png("temp_cluster_dimplot.png", width=5*100, height=5*100)
# pdf("temp_cluster_dimplot.pdf")
DimPlot(TKO, label = FALSE, group.by="fig.cluster", pt.size=0.01)# + NoLegend()
dev.off()
png("temp_patient_sample.png", width=5*100, height=5*100)
DimPlot(TKO, label = FALSE, group.by="fig.sample", pt.size=0.01)# + NoLegend()
dev.off()
TKO.markers <- FindAllMarkers(TKO, min.pct=0.2, max.cells.per.ident=500, assay="SCT")
# saveRDS(TKO.markers, "TKO.onceSCT.markers.as.fig.cluster.res02.moreCells.rds")
# TKO.markers=readRDS("TKO.onceSCT.markers.as.fig.cluster.res02.moreCells.rds")


### check celltype markers' expression to make sure every cluster's celltype
markers=TKO.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[markers$avg_log2FC>log2(1.5) & markers$pct.1>0.2, ]
### celltypeAnno() need to be run in directory "processing"
source("../functions.seurat.R")
setwd("..")
temp.celltype=celltypeAnno(TKO, markers, mouse=TRUE)
setwd("./ ... ")
# [1] "Decided cell type:"
#  [1] "Macrophage: 1.5 Macrophage_cell201801: 1.43"
#  [2] "Luminal cell: 0.67 Epithelial cell: 0.61"
#  [3] "Luminal cell: 0.66"
#  [4] "Fibroblast_cell201801: 2"
#  [5] "Luminal cell: 0.56 Epithelial cell: 0.56"
#  [6] "Neutrophil: 1.37"
#  [7] "Luminal cell: 0.27 Epithelial cell: 0.25"
#  [8] "T cell: 1.6"
#  [9] "Basal cell_profRen: 2.35"
# [10] "Endothelial cell: 4.88"
# [11] "B cell Lymph node: 2.68 B/Plasma cell_cell201801: 2.47"
# [12] "Myofibroblast: 3.41"
table(names(temp.celltype)==rownames(TKO@meta.data))
temp.celltype[temp.celltype=="B cell Lymph node"]="B/Plasma cell"
temp.celltype[temp.celltype=="Luminal cell"]="Epithelial cell"
temp.celltype[temp.celltype=="Basal cell"]="Epithelial cell"
table(names(temp.celltype)==rownames(TKO@meta.data))
table(temp.celltype)
#    B/Plasma cell Endothelial cell  Epithelial cell       Fibroblast
#              133              174             5269             1114
#       Macrophage    Myofibroblast       Neutrophil           T cell
#             1843               78              560              453
TKO@meta.data$celltype=temp.celltype
table(is.na(TKO@meta.data$celltype))
# FALSE
#  9624
# saveRDS(TKO, "TKO.moreCells.rds")
# TKO=readRDS("TKO.moreCells.rds")

png("temp_cluster_dimplot.png", width=5*100, height=5*100)
# pdf("temp_cluster_dimplot.pdf")
DimPlot(TKO, label = TRUE, group.by="fig.cluster", pt.size=0.01)# + NoLegend()
dev.off()

pdf("TKO.dimplot.celltype.pdf", width=5, height=4.5)
# pdf("dimplot.celltype.pdf")
DimPlot(TKO, label = TRUE, group.by="celltype", pt.size=0.5, raster=TRUE)# + NoLegend()
dev.off()

# png("dimplot.phase.png", width=5*100, height=5*100)
# DimPlot(TKO, label = TRUE, group.by="Phase", pt.size=0.01)# + NoLegend()
# dev.off()

png("dimplot.sample.png", width=5*100, height=5*100)
DimPlot(TKO, label = FALSE, group.by="fig.sample", pt.size=0.01)# + NoLegend()
dev.off()

png("dimplot.celltype.split.by.sample.png",width=30*100, height=5*100)
DimPlot(TKO, label = TRUE, split.by = 'fig.sample', group.by="celltype", pt.size=0.05, label.size=8)# + NoLegend()
dev.off()


### write markers of clusters as a table
markers=TKO.markers
markers=markers[order(markers$cluster, markers$avg_log2FC, decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.TKO.res02.cluster","csv",sep="."), row.name=T)


# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
require(nichenetr)
check.genes=convert_human_to_mouse_symbols(check.genes)
check.genes[c(1, 12, 13)]=c("Klkb1", "nFeature_RNA", "nFeature_SCT")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # PCa markers
# check.genes=c("FOLH1", "SCHLAP1", "AMACR","CCL2", "CD74", "LY6D", "KLK3", "AR", "ACPP")
# check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

DefaultAssay(TKO)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(TKO, features = check.genes)
dev.off()
DefaultAssay(TKO)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(TKO, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


### plot heatmap to validate cell type
# heatmap top 10 DEGs for every cluster
markers=TKO.markers
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(2), ]
# celltype.markers=readRDS("../mouse.celltype.markers.rds")
# temp.union=Reduce(union, celltype.markers)
# markers=markers[markers$gene%in%temp.union, ]
str(markers)
# 'data.frame':   411 obs. of  7 variables
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
pdf(paste("heatmap.TKO.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=nrow(top10)/8)
temp.heatmap<-DoHeatmap( subset(TKO, downsample=1000), features = top10$gene, group.by = "fig.cluster", slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### further analysis focused on epithelial cells
table(TKO@meta.data$celltype)
   # B/Plasma cell Endothelial cell  Epithelial cell       Fibroblast
   #           133              174             5269             1114
   #    Macrophage    Myofibroblast       Neutrophil           T cell
   #          1843               78              560              453
epitypes=c("Epithelial cell", "Luminal cell", "Basal cell")
TKOe=subset(TKO, subset=(celltype%in%epitypes))
# str(TKOe)
TKOe@commands=list()


dim(TKOe)
# [1] 19875  5269
DefaultAssay(TKOe)="SCT"
str(TKOe@assays$SCT@var.features)
# chr [1:3000] "Cxcl2" "Col1a1" "Ccl3" "Il1b" "Col1a2" "Sparc" "S100a9" ...
# TKOe <- FindVariableFeatures(TKOe, assay="SCT", nfeatures = 3000)
# str(TKOe@assays$SCT@var.features)
# top10 <- head(VariableFeatures(TKOe, assay="SCT"), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(TKOe, assay = "SCT")
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# toplot=plot1 + plot2
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()

TKOe=SCTransform(TKOe, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
TKOe <- RunPCA(TKOe, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(TKOe,  ndims = 100)
dev.off()

### clustering basd on PCA
{
  DefaultAssay(TKOe)="SCT"
  TKOe <- RunUMAP(TKOe, dims = 1:25, verbose = FALSE)
  TKOe <- FindNeighbors(TKOe, dims = 1:25, verbose = FALSE)

  DefaultAssay(TKOe)="SCT" 
  TKOe <- FindClusters(TKOe, verbose = FALSE, resolution = 0.2)
  table(TKOe@meta.data$seurat_clusters)
  TKOe.markers <- FindAllMarkers(TKOe, assay="SCT")
  # saveRDS(TKOe.markers, "TKOe.markers.as.fig.cluster.byPCA.moreCells.rds")
  # TKOe.markers =readRDS("TKOe.markers.as.fig.cluster.byPCA.moreCells.rds")

}

# ### clustering basd on Harmony
# {
#   require(harmony)
#   DefaultAssay(TKOe)
#   # [1] "SCT"
#   Sys.time()
#   TKOe <- RunHarmony(TKOe, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
#   dev.off()
#   Sys.time()

#   # TKO <- RunUMAP(TKO, dims = 1:25, verbose = FALSE)
#   TKOe <- RunUMAP(TKOe, reduction = "harmony", dims = 1:25, verbose = FALSE)
#   TKOe <- FindNeighbors(TKOe, reduction = "harmony", dims = 1:25, verbose = FALSE)
#   TKOe <- FindClusters(TKOe, verbose = FALSE, resolution = 0.2)
#   table(TKOe@meta.data$seurat_clusters)
#   #    0    1    2    3    4    5    6    7    8    9   10
#   # 4069 3853 2325 1659 1446 1289  499  229   88   76   43
#   TKOe.markers <- FindAllMarkers(TKOe, assay="SCT")
#   # saveRDS(TKOe.markers, "TKOe.markers.as.fig.cluster.byHarmony.rds")
#   # TKOe.markers =readRDS("TKOe.markers.as.fig.cluster.byHarmony.rds")
# }

markers=TKOe.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.TKOe.cluster.moreCells","csv",sep="."), row.name=T)

TKOe@meta.data$fig.cluster=as.character(TKOe@meta.data$seurat_clusters)
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(TKOe, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
png("dimplot.sample.png", height=5*100, width=5*100)
DimPlot(TKOe, label = FALSE, group.by="fig.sample")
dev.off()


# cluster 6 is T cell(CD3E+), remove it
# cluster 7 is B cell(MS4A1+, CD79A+), remove it
# remove clusters 8, 9 and 10 with too few cells (<100 cells)
# TKOe=subset(TKOe, cells = colnames(TKOe)[!TKOe@meta.data$fig.cluster%in%c("6", "7", "8", "9", "10")])
table(TKOe@meta.data$fig.cluster)
#    0    1    2    3    4    5    6    7
# 1896 1084 1015  679  268  204   93   30


# check known epithelial cell types
check.genes=c("KLK3", "KRT8", "EPCAM", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
check.genes=convert_human_to_mouse_symbols(check.genes)
check.genes[c(1, 12, 13)]=c("Klkb1", "nFeature_RNA", "nFeature_SCT")

# Clike
check.genes=c("IL1A", "ARL14", "FABP4", "CYP1A1",   "UPK1A", "CD55", "TGFBR3", "FGFR3", "KRT8")
check.genes=convert_human_to_mouse_symbols(check.genes)


TKOe <- SetIdent(TKOe, value = "fig.cluster")
DefaultAssay(TKOe)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(TKOe, features = check.genes, pt.size=1)
dev.off()
DefaultAssay(TKOe)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(TKOe, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()



### re-cluster TKOe_cluster0 cells
{
  TKOeC4=subset(TKOe, subset=(seurat_clusters%in%c("4")))
  dim(TKOeC4)
  # [1] 18586   267
  DefaultAssay(TKOeC4)="SCT"
  str(TKOeC4@assays$SCT@var.features)
  # chr [1:3000] "Cbr2" "Pigr" "Ly6a" "Krt15" "Nupr1" "Gsto1" "Lcn2" "Cxcl10" ...
  # TKOeC4 <- FindVariableFeatures(TKOeC4, assay="SCT", nfeatures = 3000)
  # str(TKOeC4@assays$SCT@var.features)
  # top10 <- head(VariableFeatures(TKOeC4, assay="SCT"), 10)
  # # plot variable features with and without labels
  # plot1 <- VariableFeaturePlot(TKOeC4, assay = "SCT")
  # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # toplot=plot1 + plot2
  # png("temp.png", width=10*100, height=5*100)
  # print(toplot)
  # dev.off()

  TKOeC4=SCTransform(TKOeC4, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  TKOeC4 <- RunPCA(TKOeC4, assay="SCT", verbose = FALSE, npcs = 100)
  ElbowPlot(TKOeC4,  ndims = 100)
  dev.off()

  ### clustering basd on PCA
  {
    DefaultAssay(TKOeC4)="SCT"
    # TKOeC4 <- RunUMAP(TKOeC4, dims = 1:25, verbose = FALSE)
    TKOeC4 <- FindNeighbors(TKOeC4, dims = 1:25, verbose = FALSE)

    DefaultAssay(TKOeC4)="SCT" 
    TKOeC4 <- FindClusters(TKOeC4, verbose = FALSE, resolution = 0.5)
    table(TKOeC4@meta.data$seurat_clusters)
    #   0   1
    # 199  68
    TKOeC4.markers <- FindAllMarkers(TKOeC4, assay="SCT")
    # saveRDS(TKOeC4.markers, "TKOeC4.markers.as.fig.cluster.byPCA.moreCells.rds")
    # TKOeC4.markers =readRDS("TKOeC4.markers.as.fig.cluster.byPCA.moreCells.rds")
  }

  ### clustering basd on Harmony
  # {
  #   require(harmony)
  #   DefaultAssay(TKOeC4)
  #   # [1] "SCT"
  #   Sys.time()
  #   TKOeC4 <- RunHarmony(TKOeC4, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
  #   dev.off()
  #   Sys.time()

  #   # TKO <- RunUMAP(TKO, dims = 1:25, verbose = FALSE)
  #   TKOeC4 <- RunUMAP(TKOeC4, reduction = "harmony", dims = 1:25, verbose = FALSE)
  #   TKOeC4 <- FindNeighbors(TKOeC4, reduction = "harmony", dims = 1:25, verbose = FALSE)
  #   TKOeC4 <- FindClusters(TKOeC4, verbose = FALSE, resolution = 0.2)
  #   table(TKOeC4@meta.data$seurat_clusters)
  #    0    1    2    3    4    5    6    7
  # 1905 1074 1011  680  267  205   83   30
  #   # TKOeC4.markers <- FindAllMarkers(TKOeC4, assay="SCT")
  #   # saveRDS(TKOeC4.markers, "TKOeC4.markers.as.fig.cluster.byHarmony.rds")
  #   # TKOeC4.markers =readRDS("TKOeC4.markers.as.fig.cluster.byHarmony.rds")
  # }

  TKOeC4@meta.data$fig.cluster=as.character(TKOeC4@meta.data$seurat_clusters)
  png("dimplot.cluster.png", height=5*100, width=5*100)
  DimPlot(TKOeC4, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  dev.off()
  png("dimplot.sample.png", height=5*100, width=5*100)
  DimPlot(TKOeC4, label = TRUE, group.by="fig.sample")
  dev.off()
  # png("dimplot.gleason.png", height=5*100, width=5*100)
  # DimPlot(TKOeC4, label = TRUE, group.by="gleason")
  # dev.off()
  # png("dimplot.pheno.png", height=5*100, width=5*100)
  # DimPlot(TKOeC4, label = TRUE, group.by="pheno")
  # dev.off()


  # check known epithelial cell types
  check.genes=c("KLK3", "KRT8", "EPCAM", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
  check.genes=convert_human_to_mouse_symbols(check.genes)
  check.genes[c(1, 12, 13)]=c("Klkb1", "nFeature_RNA", "nFeature_SCT")

  # Clike
  check.genes=c("IL1A", "ARL14", "FABP4", "CYP1A1",   "UPK1A", "CD55", "TGFBR3", "FGFR3", "KRT8")
  library(nichenetr)
  check.genes=convert_human_to_mouse_symbols(check.genes)


  TKOeC4 <- SetIdent(TKOeC4, value = "fig.cluster")
  DefaultAssay(TKOeC4)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(TKOeC4, features = check.genes, pt.size=1)
  dev.off()
  DefaultAssay(TKOeC4)="SCT"
  png("temp.VlnPlot.png.png", width=10*100, height=8*100)
  temp=VlnPlot(TKOeC4, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()


  # ### name clusters
  current.cluster.ids <- c(0:1)
  new.cluster.ids <- c("4.0", "4.1")
  # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
  TKOeC4@meta.data$fig.cluster <- plyr::mapvalues(x = as.character(TKOeC4@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  pdf("TKOeC4.dimplot.fig.cluster.pdf", height=4.5, width=5)
  DimPlot(TKOeC4, label = TRUE, group.by="fig.cluster", raster=TRUE)
  dev.off()

  c0type=TKOeC4@meta.data$fig.cluster
  names(c0type)=rownames(TKOeC4@meta.data)
  table(c0type)
  # 4.0 4.1
  # 192  76

  ### TKOe annotation
  # add refined basal/typeC types into epitype slot
  TKOe@meta.data$fig.cluster=as.character(TKOe@meta.data$seurat_clusters)
  TKOe@meta.data[names(c0type), "fig.cluster"]=c0type
  table(TKOe@meta.data$fig.cluster)
  #    0    1    2    3  4.0  4.1    5    6    7
  # 1905 1074 1011  680  199   68  205   83   30
  TKOe@meta.data$fig.cluster=factor(TKOe@meta.data$fig.cluster, levels=c("0", "1", "2", "3", "4.0", "4.1", "5", "6", "7"))
  table(is.na( TKOe@meta.data$seurat_clusters))
  # FALSE
  #  5255
  pdf("TKOe.dimplot.fig.cluster.pdf", height=4.5, width=5)
  DimPlot(TKOe, label = TRUE, group.by="fig.cluster", raster=TRUE)
  dev.off()

  TKOe <- SetIdent(TKOe, value = "fig.cluster")
  table(TKOe@ active.ident)
  #    0    1    2    3  4.0  4.1    5    6    7
  # 1905 1074 1011  680  199   68  205   83   30
  TKOe.markers <- FindAllMarkers(TKOe, assay="SCT", max.cells.per.ident=500)
  # saveRDS(TKOe.markers, "TKOe.markers.as.fig.cluster.byHarmony.rds")
  # TKOe.markers =readRDS("TKOe.markers.as.fig.cluster.byHarmony.rds")

  # indi.plot(TKOeC4, tocheck="fig.patient")
  # indi.plot(TKOeC4, tocheck="fig.cluster")

  # saveRDS(TKOeC4, "TKOeC4.rds")
  # TKOeC4=readRDS("TKOeC4.rds")

}


TKOe@meta.data$celltype="NE"
TKOe@meta.data$celltype[TKOe@meta.data$fig.cluster=="4.1"]="Clike"
TKOe@meta.data$celltype[TKOe@meta.data$fig.cluster=="4.0"]="Intermediate"
TKOe@meta.data$celltype=factor(TKOe@meta.data$celltype, levels=c("Clike", "Intermediate", "NE"))
# png("dimplot.gleason.png", height=5*100, width=5*100)
# DimPlot(TKOe, label = TRUE, group.by="gleason")
# dev.off()
# png("dimplot.pheno.png", height=5*100, width=5*100)
# DimPlot(TKOe, label = TRUE, group.by="pheno")
# dev.off()


# saveRDS(TKOe, "TKOe.moreCells.rds")
# TKOe=readRDS("TKOe.moreCells.rds")


pdf("TKOe.dimplot.celltype.pdf", height=5, width=6)
DimPlot(TKOe, label = TRUE, group.by="celltype", label.size=5, raster=TRUE)
dev.off()


# check known epithelial cell types
check.genes=c("KLK3", "KRT8", "EPCAM", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
check.genes=convert_human_to_mouse_symbols(check.genes)
check.genes[c(1, 12, 13)]=c("Klkb1", "nFeature_RNA", "nFeature_SCT")

# Clike
check.genes=c("IL1A", "ARL14", "FABP4", "CYP1A1",   "UPK1A", "CD55", "TGFBR3", "FGFR3", "KRT8")
check.genes=convert_human_to_mouse_symbols(check.genes)

check.genes=c("Chga", "Nse",  "Il1a", "Arl14", "Cd55", "Upk1a", "nFeature_RNA", "Krt4", "Psca")
check.genes=c("Rmst", "Ascl1", "Neurod1", "Chga")
check.genes=c( "Il1rn", "Il1a", "Psca", "Tacstd2", "Krt4", "Ar", "Krt5", "Krt19", "Chga", "Syp", "Ncam1", "Epcam", "Ceacam1", "nFeature_RNA" )


TKOe <- SetIdent(TKOe, value = "celltype") # "fig.cluster" , "fig.celltype", "celltype"
DefaultAssay(TKOe)="SCT"
png("temp.featureplot.png", height=10*100, width=10*1100)
FeaturePlot(TKOe, features = check.genes, pt.size=1)
dev.off()
DefaultAssay(TKOe)="SCT"
pdf("TKOe.VlnPlot.pdf", width=5, height=15)
temp=VlnPlot(TKOe, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 2, colour = "black", shape = 95)
}
print(temp)
dev.off()


check.genes=c( "Il1a", "Arl14", "Fabp4", "Snx31", "Cyp1a1", "Fgfr3", "Tgfbr3", "Ncam1", "Ceacam1", "Krt19")
check.genes=c("CD55", "IL1RN", "IL1A", "UPK1A", "UPK1B", "PIM1", "PMAIP1")

str(epi.tumor.clike.markers) # created by part11_mergeEpi.V3_excludeIncidente.R
check.genes=epi.tumor.clike.markers
library(nichenetr)
check.genes=convert_human_to_mouse_symbols(check.genes)
check.genes
#   S100A4     IL1A    IL1RN    UPK1B     <NA>     AREG     HPGD   SPINK1
# "S100a4"   "Il1a"  "Il1rn"  "Upk1b"       NA   "Areg"   "Hpgd" "Spink1"
#    FOXQ1   TRIM31     PIM1    FABP4     PSCA     UPK2     <NA>     SNCG
#  "Foxq1" "Trim31"   "Pim1"  "Fabp4"   "Psca"   "Upk2"       NA   "Sncg"
#     <NA>    DHRS2     GPX2     BMP2    UPK1A    UPK3A    FBLN1     CSTB
#       NA  "Dhrs2"   "Gpx2"   "Bmp2"  "Upk1a"  "Upk3a"  "Fbln1"   "Cstb"
check.genes[5]="Gm40318"
check.genes[15]="Akr1c21"

TKOe<- SetIdent(TKOe, value = "celltype") # "fig.cluster" , "fig.celltype"
temp=TKOe
temp=ScaleData(temp, features = rownames(temp@assays$RNA@counts), assay="RNA")
temp=subset(temp, downsample=300)
# temp=SCTransform(temp, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
# The following features were omitted as they were not found in the scale.data slot for the SCT assay: Akr1c21, Padi3
# levels(temp)=c("Clike", "TypeC", "Basal", "Luminal")
pdf(paste("heatmap.TKOe.Clike.pdf", sep=""),  width=5, height=5)#height=length(tocheck)/3
temp.heatmap<-DoHeatmap( temp, features = check.genes, slot="scale.data", assay="RNA", angle = 0 , size=2)#+ scale_fill_gradientn(colors = c("green", "red"))
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


# ### name clusters
# current.cluster.ids <- c(0:4)
# new.cluster.ids <- c("tC20", "tB21", "tL22", "tC23", "tC24")
# # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
# TKOe@meta.data$epitype <- plyr::mapvalues(x = as.character(TKOe@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
# png("dimplot.epitype.png", height=5*100, width=5*100)
# DimPlot(TKOe, label = TRUE, group.by="epitype", label.size=8)
# dev.off()

# btype=TKOe@meta.data$epitype
# names(btype)=rownames(TKOe@meta.data)
# table(btype)

# indi.plot(TKOe, tocheck="fig.patient")
# indi.plot(TKOe, tocheck="fig.cluster")

# saveRDS(TKOe, "TKOe.rds")
# TKOe=readRDS("TKOe.rds")



### celltype-specific DEG identification
TKOe@meta.data$celltype=rep("NE", nrow(TKOe@meta.data))
TKOe@meta.data$celltype[TKOe@meta.data$fig.cluster=="4.1"]="Clike"
TKOe@meta.data$celltype[TKOe@meta.data$fig.cluster=="4.0"]="Intermediate"
table(TKOe@meta.data$celltype)
# Clike Intermediate           NE
#    68          199         4988
DefaultAssay(TKOe)="SCT"
### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(TKOe@meta.data$celltype)
celltype=sort(celltype)
TKOe = SetIdent(TKOe, value = "celltype")
for(i in celltype)
{
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        markerlist[[temp.name]]=FindMarkers(TKOe, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
      }
  }
}
# saveRDS(markerlist, "TKOe.celltype.pairwise.rds")
# markerlist=readRDS("TKOe.celltype.pairwise.rds")

### canonical process to find specific markers of Clike cells
{
  markerlist=readRDS("TKOe.celltype.pairwise.rds")
  sp.markerlist=list()
  celltype=unique(TKOe@meta.data$celltype)
  celltype=sort(celltype)
  for(i in celltype)
  {
    temp.genes=rownames(TKOe@assays$SCT@data)
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          cur.list=markerlist[[temp.name]]
          cur.list=cur.list[order(cur.list$avg_log2FC, decreasing=TRUE), ]
          # cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(2) & cur.list$pct.1>0.3 & cur.list$pct.2<0.1,]) # canonical version
          cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(1.5), ])
          temp.genes=intersect(temp.genes, cur.genes)
        }
    }
    sp.markerlist[[i]]=temp.genes
  }
  str(sp.markerlist)
  sp.markerlist[["Luminal"]]=c(sp.markerlist[["Luminal"]], "KLK3", "AR")

  temp=subset(TKOe, downsample=300)
  sp.markerlist=sp.markerlist[c("Clike", "Intermediate", "NE")]
  tocheck=Reduce(union, sp.markerlist)
  temp <- SetIdent(temp, value = "celltype")
  DefaultAssay(temp)="RNA"
  temp=NormalizeData(temp)
  temp= ScaleData(temp, features = rownames(temp), assay="RNA" )
  levels(temp)=c("Clike", "Intermediate", "NE")
  pdf(paste("heatmap.TKOe.typeC.pdf", sep=""),  width=10, height=length(tocheck)/3)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="RNA", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
  DefaultAssay(temp)="SCT"

  TKOe.clike.markers=sp.markerlist$Clike
  str(TKOe.clike.markers)
}


# best=c("CD55", "IL1RN", "IL1A", "UPK1A", "UPK1B", "PIM1", "PMAIP1")
tocheck=epi.tumor.clike.markers
library(nichenetr)
tocheck=convert_human_to_mouse_symbols(tocheck)
intersect(tocheck, TKOe.clike.markers)
# S100P --> Gm40318  # S100p
# AKR1C2 --> Akr1c21  # AI315367; 9430025F20Rik
"Gm40318"%in%rownames(temp@assays$RNA@counts)
"S100p"%in%rownames(temp@assays$RNA@counts)
"kr1c21"%in%rownames(temp@assays$RNA@counts)
"AI315367"%in%rownames(temp@assays$RNA@counts)
"9430025F20Rik"%in%rownames(temp@assays$RNA@counts)

### DEGs and heatmap in clusters
{
  TKOe.markers =readRDS("TKOe.markers.as.fig.cluster.byHarmony.rds")
  TKOe.markers$cluster=as.character(TKOe.markers$cluster)
  TKOe.markers$cluster=factor(TKOe.markers$cluster, levels=c("0", "1", "2", "3", "4.0", "4.1", "5", "6", "7"))
  TKOe@meta.data$fig.cluster=as.character(TKOe@meta.data$fig.cluster)
  TKOe@meta.data$fig.cluster=factor(TKOe@meta.data$fig.cluster, levels=c("0", "1", "2", "3", "4.0", "4.1", "5", "6", "7"))

  ### heatmap
  markers=TKOe.markers
  markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
  markers$foldchange=2^(markers$avg_log2FC)
  markers=markers[order(markers$cluster, markers$foldchange, decreasing=c(FALSE, TRUE)), ]

  write.csv(markers[,c("foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.TKOe.markers.as.fig.cluster","csv",sep="."), row.name=T)

  # & markers$pct.1>0.3 & markers$pct.2<0.1
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = (avg_log2FC) )
  temp=subset(TKOe, downsample=300)
  temp <- SetIdent(temp, value = "fig.cluster")
  tocheck=top10$gene
  pdf(paste("heatmap.TKOe.celltype.tocheck.pdf", sep=""),  width=(length(unique(tocheck))/4+5), height=length(tocheck)/8)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()


  tocheck=c( "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1" )
  tocheck=c( "Il1a", "Arl14", "Fabp4", "Snx31", "Cyp1a1", "Fgfr3", "Tgfbr3")

  prime<- SetIdent(prime, value = "celltype")
  temp=subset(prime, downsample=300)
  levels(temp)=c("Clike", "TypeC", "Basal", "Luminal")
  pdf(paste("heatmap.prime.Clike.pdf", sep=""),  width=5, height=5)#height=length(tocheck)/3
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)#+ scale_fill_gradientn(colors = c("green", "red"))
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

}


### DEGs and heatmap in clusters 4.0 and 4.1
{
  TKOe=SetIdent(TKOe, value="fig.cluster")
  markers.4=FindMarkers(TKOe, ident.1 = "4.1", ident.2 = "4.0", assay="SCT")
  # saveRDS(markers.4, "TKOe.markers.c41v40.as.fig.cluster.byHarmony.rds")
  markers.4 =readRDS("TKOe.markers.c41v40.as.fig.cluster.byHarmony.rds")

  ### heatmap
  markers=markers.4
  markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & markers$pct.1>0.1, ]
  markers$foldchange=2^(markers$avg_log2FC)
  markers=markers[order(markers$foldchange, decreasing=TRUE), ]
  write.csv(markers[,c("foldchange","pct.1", "pct.2", "p_val_adj")], file=paste("markers.TKOe.markers.c41v40.as.fig.cluster","csv",sep="."), row.name=T)

  # & markers$pct.1>0.3 & markers$pct.2<0.1
  up10=rownames(markers)[1:10]
  dn10=rownames(markers)[(nrow(markers)-9):nrow(markers)]
  top10=c(up10, dn10)
  temp=subset(TKOe, subset=fig.cluster%in%c("4.0", "4.1"), downsample=300)
  temp <- SetIdent(temp, value = "fig.cluster")
  tocheck=top10
  pdf(paste("heatmap.TKOe.celltype.tocheck.pdf", sep=""),  width=(length(unique(tocheck))/4+5), height=length(tocheck)/8)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
}



### enrichment
{
  ### qusage to find activated pathways
  require(qusage)
  packageVersion("qusage")
  # [1] '2.24.0'
  require(msigdbr)
  m_df = msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
  m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.reactome = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
  m_t2g.reactome = m_df.reactome %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.kegg = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
  m_t2g.kegg = m_df.kegg %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.hall = msigdbr(species = "Mus musculus", category = "H")
  m_t2g.hall = m_df.hall %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.regulate = msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:GTRD")
  m_t2g.regulate = m_df.regulate %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  # convert annotation from data.frame to list format
  kegg.list=list()
  kegg.name=unique(m_t2g.kegg$gs_name)
  for(i in kegg.name)
  {
    kegg.list[[i]]=m_t2g.kegg$gene_symbol[m_t2g.kegg$gs_name==i]
  }
  hall.list=list()
  hall.name=unique(m_t2g.hall$gs_name)
  for(i in hall.name)
  {
    hall.list[[i]]=m_t2g.hall$gene_symbol[m_t2g.hall$gs_name==i]
  }
    #seurat.all = a seurat object from scran
    #nm = name of the run
    #gs = gene set to test, a list
  # the "run_qusage_heatmap" function needs that @meta.data$for_temp_color is a numeric format vector with MIN=1 rather than MIN=0, whic represent the cluster.ID
  TKOe@meta.data$for_temp_color=as.numeric(TKOe@meta.data$fig.cluster)
  run_qusage_heatmap.seurat3(TKOe, nm = 'fig.cluster', kegg.list, my.seed=100)
  run_qusage_heatmap.seurat3(TKOe, nm = 'fig.cluster', hall.list, my.seed=100)

}


### plot heatmap for special categories of DEGs
setwd("..")
subtype_markers=subtypeMarker(TKOe, mouse=TRUE)
setwd("./ ... ")
TKOe.markers =readRDS("TKOe.markers.as.fig.cluster.byHarmony.rds")
TKOe.markers$cluster=as.character(TKOe.markers$cluster)
TKOe.markers$cluster=factor(TKOe.markers$cluster, levels=c("0", "1", "2", "3", "4.0", "4.1", "5", "6", "7"))

markers=TKOe.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.2, ]
# for TFs, change log2(1.5) to log(1.2)
# TKOe <- ScaleData(TKOe, features = rownames(TKOe), assay="RNA" )
num.cluster=length(unique(markers$cluster))
tocheck=c( "chemokine.ligand",
      "chemokine.receptor",
      "TF_fromLieBing",
      "active.check.point.ligand",
      "inhibi.check.point.ligand",
      "communication_ligand",
      "communication_receptor" )
for(j in tocheck)
{
  tempcheck=subtype_markers[[j]]
  temp.markers=markers[markers$gene%in%tempcheck, ]
  top10 <- temp.markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
  top10=top10[order(top10$cluster), ]
  if(nrow(top10)>0)
  {
    DefaultAssay(TKOe)="SCT"
    pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
    temp=subset(TKOe, downsample=300)
    temp <- SetIdent(temp, value = "fig.cluster")
    temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()
    DefaultAssay(TKOe)="SCT"
  }
}


### find Clike (cluster 4.1) specific markers 
{
  ### annotate Clike cells in all-cell-types data
  # TKO=readRDS("TKO.moreCells.rds")
  TKOeall=subset(TKO, fig.sample=="pten.p53.rb")
  # TKOe=readRDS("TKOe.rds")
  epic=paste("Epi.C", TKOe@meta.data$fig.cluster, sep="")
  names(epic)=rownames(TKOe@meta.data)
  table(names(epic)%in%colnames(TKOeall))
  TKOeall@meta.data$finer.celltype=TKOeall@meta.data$celltype
  TKOeall@meta.data[names(epic), "finer.celltype"]=epic
  # saveRDS(TKOeall, "TKOeall.rds")
  # TKOeall=readRDS("TKOeall.rds")
  
  png("temp_celltype_dimplot.png", width=5*100, height=5*100)
  # pdf("temp_cluster_dimplot.pdf")
  DimPlot(TKOeall, label = TRUE, group.by="finer.celltype", pt.size=0.01)# + NoLegend()
  dev.off()
  pdf("temp_celltype_dimplot.pdf", width=6, height=5)
  # pdf("temp_cluster_dimplot.pdf")
  DimPlot(TKOeall, label = FALSE, group.by="finer.celltype", pt.size=0.01)# + NoLegend()
  dev.off()
  png("temp_patient_sample.png", width=5*100, height=5*100)
  DimPlot(TKOeall, label = FALSE, group.by="fig.sample", pt.size=0.01)# + NoLegend()
  dev.off()



  check.genes=c("Rmst", "Ascl1", "Neurod1", "Chga")

  TKOeall@meta.data$finer.celltype=factor(TKOeall@meta.data$finer.celltype, levels=c(sort(unique(TKOeall@meta.data$finer.celltype))))
  TKOeall <- SetIdent(TKOeall, value = "finer.celltype")
  DefaultAssay(TKOeall)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(TKOeall, features = check.genes, pt.size=1)
  dev.off()
  DefaultAssay(TKOeall)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=8*100)
  temp=VlnPlot(TKOeall, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  pdf("temp.VlnPlot.pdf", width=10, height=8)
  temp=VlnPlot(TKOeall, features = check.genes, pt.size = 0, group.by="finer.celltype", ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()


  ### scatter plot
  # TKOeall.bak=TKOeall
  # TKOeall=TKOeall.bak
  celltypes=unique(TKOeall@meta.data$celltype)
  celltypes=celltypes[celltypes!="Epithelial cell"]
  TKOeall=subset(TKOeall, celltype%in%celltypes)
  TKOeall = SetIdent(TKOeall, value = "celltype")
  pdf(paste("temp.featureScatter.pdf", sep=""), width = 15, height =4)
  plot1 <- FeatureScatter(TKOeall, feature1 = "Neurod1", feature2 = "Ascl1")
  # plot1=plot1+geom_vline(xintercept=500)+geom_vline(xintercept=1000)+geom_hline(yintercept=60)+geom_hline(yintercept=40)+geom_hline(yintercept=30)
  plot2 <- FeatureScatter(TKOeall, feature1 = "Neurod1", feature2 = "Rmst")
  plot3 <- FeatureScatter(TKOeall, feature1 = "Rmst", feature2 = "Ascl1")
  print(plot1 + plot2 + plot3)
  dev.off()


  ### cross-comparison to identify specific markers 
  DefaultAssay(TKOeall)="SCT"
  ### pairwise customized cluster's DEGs
  markerlist=list()
  celltype=unique(TKOeall@meta.data$finer.celltype)
  celltype=sort(celltype)
  TKOeall = SetIdent(TKOeall, value = "finer.celltype")
  for(i in celltype)
  {
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          markerlist[[temp.name]]=FindMarkers(TKOeall, ident.1 = i, ident.2 = j, slot="data", assay="SCT", max.cells.per.ident=300)
        }
    }
  }
  # saveRDS(markerlist, "TKOeall.epiCandCelltype.pairwise.rds")


  # @data
  markerlist=readRDS("TKOeall.epiCandCelltype.pairwise.rds")
  # @scale.data
  # markerlist=readRDS("markerlist.TKOeall.pairwise.initiating.types.scaled.data.rds")
  ######
  celltypes=sort(unique(TKOeall@meta.data$finer.celltype))
  markerlist.filter=list()
  for(j in celltypes)
  {
    print("celltype:")
    print(j)
    print("------------>>>>>>>>>>>")
    sp.markers=rownames(TKOeall@assays$SCT@data)
    sp.markerlist="1"
    for(i in names(markerlist))
    {
      temp=strsplit(i, "_")[[1]][1]
      if(temp==j)
      {
        markers=markerlist[[i]]
        markers=markers[markers$p_val_adj<0.05, ]
        markers=markers[order(markers$avg_log2FC, decreasing = TRUE), ]
        markers=markers[markers$avg_log2FC>log2(1.5), ]
        # markers=markers[markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
        # markers=markers[markers$avg_log2FC>log2(2) , ]
        sp.markers=intersect(sp.markers, rownames(markers))
        temp.sp.markerlist=markers[sp.markers, "avg_log2FC"]
        names(temp.sp.markerlist)=sp.markers
        if(is.character(sp.markerlist))
        {
          sp.markerlist=temp.sp.markerlist
          # names(sp.markerlist)=names(temp.sp.markerlist)
        }else
        {
          temp.cbind=cbind(temp.sp.markerlist, sp.markerlist[sp.markers])
          sp.markerlist=apply(temp.cbind, 1, min)
        }
        print("pairwise:")
        print(i)
        print("length(sp.markers):")
        print(length(sp.markers))
      }
    }
    sp.markerlist=sp.markerlist[order(sp.markerlist, decreasing=TRUE)]
    markerlist.filter[[j]]=names(sp.markerlist)
  }
  str(markerlist.filter)
  write.table(markerlist.filter[["Epi.C4.1"]], file="Epi.C4.1", col.name=F, row.name=F, quote=F, sep='\t')
  for(i in names(markerlist.filter))
  {
    # markerlist.filter[[i]]=markerlist.filter[[i]][1:50]
    markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
  }
  markerlist.filter=markerlist.filter["Epi.C4.1"]
  markers=Reduce(union, markerlist.filter)
  str(markerlist.filter)
  str(markers)

  pdf(paste("temp.TKOe.Clike.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/7 )
  TKOeall = SetIdent(TKOeall, value = "finer.celltype")
  temp.heatmap<-DoHeatmap( subset(TKOeall, downsample=300), features = markers, group.by = "finer.celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

  # @data-based markers
  # write.table(markers, file="marker.union.TKOeall.pairwise.initiating.types.txt", col.name=F, row.name=F, quote=F, sep='\t')
  # @scale.data-based markers:
  write.table(markers, file="markerlist.TKOeall.pairwise.celltype.txt", col.name=F, row.name=F, quote=F, sep='\t')

  ### keep markers overlaped with human markers
  human.markers=read.table("markerlist.prim.pairwise.celltype.txt")$V1
  require(nichenetr)
  human.markers=convert_human_to_mouse_symbols(human.markers)
  human.markers[13]="Gm40318" #"S100P"
  markers=markers[markers%in%human.markers]
  pdf(paste("temp.TKOe.Clike.overlaped.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/7 )
  TKOeall = SetIdent(TKOeall, value = "finer.celltype")
  temp.heatmap<-DoHeatmap( subset(TKOeall, downsample=300), features = markers, group.by = "finer.celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

}


### find General Clike (cluster 4.0+4.1) specific markers 
{
  ### annotate Clike cells in all-cell-types data
  # TKO=readRDS("TKO.moreCells.rds")
  TKOeall=subset(TKO, fig.sample=="pten.p53.rb")
  # TKOe=readRDS("TKOe.rds")
  epic=paste("Epi.C", TKOe@meta.data$fig.cluster, sep="")
  epic[epic%in%c("Epi.C4.0", "Epi.C4.1")]="Epi.C4"
  names(epic)=rownames(TKOe@meta.data)
  table(names(epic)%in%colnames(TKOeall))
  TKOeall@meta.data$finer.celltype=TKOeall@meta.data$celltype
  TKOeall@meta.data[names(epic), "finer.celltype"]=epic
  # saveRDS(TKOeall, "TKOeall.rds")
  # TKOeall=readRDS("TKOeall.rds")
  
  png("temp_celltype_dimplot.png", width=5*100, height=5*100)
  # pdf("temp_cluster_dimplot.pdf")
  DimPlot(TKOeall, label = TRUE, group.by="finer.celltype", pt.size=0.01)# + NoLegend()
  dev.off()
  png("temp_patient_sample.png", width=5*100, height=5*100)
  DimPlot(TKOeall, label = FALSE, group.by="fig.sample", pt.size=0.01)# + NoLegend()
  dev.off()

  ### cross-comparison to identify specific markers 
  DefaultAssay(TKOeall)="SCT"
  ### pairwise customized cluster's DEGs
  markerlist=list()
  celltype=unique(TKOeall@meta.data$finer.celltype)
  celltype=sort(celltype)
  TKOeall = SetIdent(TKOeall, value = "finer.celltype")
  for(i in celltype)
  {
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          markerlist[[temp.name]]=FindMarkers(TKOeall, ident.1 = i, ident.2 = j, slot="data", assay="SCT", max.cells.per.ident=300)
        }
    }
  }
  # saveRDS(markerlist, "TKOeall.GeneralEpiCandCelltype.pairwise.rds")


  # @data
  markerlist=readRDS("TKOeall.GeneralEpiCandCelltype.pairwise.rds")
  # @scale.data
  # markerlist=readRDS("markerlist.TKOeall.pairwise.initiating.types.scaled.data.rds")
  ######
  ### only consider epithelial cluster to find surface marker
  {
    # TKOeall.bak=TKOeall
    # TKOeall=TKOeall.bak
    # TKOeall=subset(TKOeall, celltype=="Epithelial cell")
  }
  celltypes=sort(unique(TKOeall@meta.data$finer.celltype))
  markerlist.filter=list()
  for(j in celltypes)
  {
    print("celltype:")
    print(j)
    print("------------>>>>>>>>>>>")
    sp.markers=rownames(TKOeall@assays$SCT@data)
    sp.markerlist="1"
    for(i in names(markerlist))
    {
      temp=strsplit(i, "_")[[1]][1]
      temp2=strsplit(i, "_")[[1]][3]
      # if(temp==j)
      if((temp==j)&(temp2%in%celltypes))
      {
        markers=markerlist[[i]]
        markers=markers[markers$p_val_adj<0.05, ]
        markers=markers[order(markers$avg_log2FC, decreasing = TRUE), ]
        markers=markers[markers$avg_log2FC>log2(1.1), ]
        # markers=markers[markers$pct.1>0.3, ]
        # markers=markers[markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
        # markers=markers[markers$avg_log2FC>log2(2) , ]
        sp.markers=intersect(sp.markers, rownames(markers))
        temp.sp.markerlist=markers[sp.markers, "avg_log2FC"]
        names(temp.sp.markerlist)=sp.markers
        if(is.character(sp.markerlist))
        {
          sp.markerlist=temp.sp.markerlist
          # names(sp.markerlist)=names(temp.sp.markerlist)
        }else
        {
          temp.cbind=cbind(temp.sp.markerlist, sp.markerlist[sp.markers])
          sp.markerlist=apply(temp.cbind, 1, min)
        }
        print("pairwise:")
        print(i)
        print("length(sp.markers):")
        print(length(sp.markers))
      }
    }
    sp.markerlist=sp.markerlist[order(sp.markerlist, decreasing=TRUE)]
    markerlist.filter[[j]]=names(sp.markerlist)
  }
  str(markerlist.filter)
  write.table(markerlist.filter[["Epi.C4"]], file="Epi.C4.txt", col.name=F, row.name=F, quote=F, sep='\t')
  for(i in names(markerlist.filter))
  {
    # markerlist.filter[[i]]=markerlist.filter[[i]][1:50]
    markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
  }
  markerlist.filter=markerlist.filter["Epi.C4"]
  markers=Reduce(union, markerlist.filter)
  str(markerlist.filter)
  str(markers)

  ### only keep surface marker
  require(msigdbr)
  cc = msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC")
  cc = cc %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  cdgenes=sort(unique(cc$gene_symbol[cc$gs_name=="GO_CELL_SURFACE"]))
  table(markers%in%cdgenes)

  markers=markers[markers%in%cdgenes]

  pdf(paste("temp.TKOe.general.Clike.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/7+3 )
  TKOeall = SetIdent(TKOeall, value = "finer.celltype")
  temp.heatmap<-DoHeatmap( subset(TKOeall, downsample=300), features = markers, group.by = "finer.celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

  # @data-based markers
  # write.table(markers, file="marker.union.TKOeall.pairwise.initiating.types.txt", col.name=F, row.name=F, quote=F, sep='\t')
  # @scale.data-based markers:
  write.table(markers, file="markerlist.general.C4.TKOeall.pairwise.celltype.txt", col.name=F, row.name=F, quote=F, sep='\t')

  ### keep markers overlaped with human markers
  human.markers=read.table("markerlist.prim.pairwise.celltype.txt")$V1
  require(nichenetr)
  human.markers=convert_human_to_mouse_symbols(human.markers)
  # human.markers[13]="Gm40318" #"S100P"
  markers=markers[markers%in%human.markers]
  pdf(paste("temp.TKOe.Clike.overlaped.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/7+3 )
  TKOeall = SetIdent(TKOeall, value = "finer.celltype")
  temp.heatmap<-DoHeatmap( subset(TKOeall, downsample=300), features = markers, group.by = "finer.celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

}



### Clike/typeC/basal/luminal signature to identify cell type
{
  ### read markers from TKOeallary PCa data
  {
    # TKOe=readRDS("TKOe.rds")
    # @data
    markerlist=readRDS("TKOeall.celltype.pairwise.rds")
    # @scale.data
    # markerlist=readRDS("markerlist.TKOeall.pairwise.initiating.types.scaled.data.rds")
    ######
    celltypes=sort(unique(TKOe@meta.data$celltype))
    markerlist.filter=list()
    for(j in celltypes)
    {
      print("celltype:")
      print(j)
      print("------------>>>>>>>>>>>")
      sp.markers=rownames(TKOe@assays$SCT@data)
      sp.markerlist="1"
      for(i in names(markerlist))
      {
        temp=strsplit(i, "_")[[1]][1]
        if(temp==j)
        {
          markers=markerlist[[i]]
          markers=markers[markers$p_val_adj<0.05, ]
          markers=markers[order(markers$avg_log2FC, decreasing = TRUE), ]
          markers=markers[markers$avg_log2FC>log2(1.5), ]
          # markers=markers[markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
          # markers=markers[markers$avg_log2FC>log2(2) , ]
          sp.markers=intersect(sp.markers, rownames(markers))
          temp.sp.markerlist=markers[sp.markers, "avg_log2FC"]
          names(temp.sp.markerlist)=sp.markers
          if(is.character(sp.markerlist))
          {
            sp.markerlist=temp.sp.markerlist
            # names(sp.markerlist)=names(temp.sp.markerlist)
          }else
          {
            temp.cbind=cbind(temp.sp.markerlist, sp.markerlist[sp.markers])
            sp.markerlist=apply(temp.cbind, 1, min)
          }
          print("pairwise:")
          print(i)
          print("length(sp.markers):")
          print(length(sp.markers))
        }
      }
      sp.markerlist=sp.markerlist[order(sp.markerlist, decreasing=TRUE)]
      markerlist.filter[[j]]=names(sp.markerlist)
    }
    str(markerlist.filter)
    # $ Basal  : chr [1:7] "KRT5" "KRT15" "DSC3" "DST" ...
    # $ Clike  : chr [1:82] "SPINK1" "DHRS2" "UPK1A" "UPK2" ...
    # $ Luminal: chr [1:100] "KLK3" "NPY" "KLK2" "SMS" ...
    # $ TypeC  : chr [1:23] "SCGB1A1" "LTF" "OLFM4" "LCN2" ...
    for(i in names(markerlist.filter))
    {
      markerlist.filter[[i]]=markerlist.filter[[i]][1:7]
      markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
    }
    markers=Reduce(union, markerlist.filter)
    str(markerlist.filter)
    str(markers)
  }

  ### get detail of marker expression
  check.genes=markers
  check.genes=convert_human_to_mouse_symbols(check.genes)
  check.genes[17]="Klk4"
  TKOe <- SetIdent(TKOe, value = "fig.cluster")
  DefaultAssay(TKOe)="SCT"
  png("temp.featureplot.png", height=10*100, width=18*100)
  FeaturePlot(TKOe, features = check.genes, pt.size=1,ncol=7)
  dev.off()
  DefaultAssay(TKOe)="SCT"
  png("temp.VlnPlot.png.png", width=10*100, height=8*100)
  temp=VlnPlot(TKOe, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off() 

  check.genes=c("Gprc5a", "Tacstd2", "Epcam", "Cbr2", "Epcam", "Ceacam1")
  # check.genes=convert_human_to_mouse_symbols(check.genes)
  TKOeall <- SetIdent(TKOeall, value = "finer.celltype")
  DefaultAssay(TKOeall)="RNA"
  png("temp.featureplot.png", height=15*100, width=12*100)
  FeaturePlot(TKOeall, features = check.genes, pt.size=1,ncol=2)
  dev.off()

  check.genes=c("Gprc5a", "Tacstd2", "Epcam", "Cbr2", "Epcam", "Ceacam1")
  # check.genes=convert_human_to_mouse_symbols(check.genes)
  TKOe <- SetIdent(TKOe, value = "fig.celltype")
  DefaultAssay(TKOe)="RNA"
  png("temp.featureplot2.png", height=10*100, width=15*100)
  FeaturePlot(TKOe, features = check.genes, pt.size=1,ncol=3)
  dev.off()

  c41=subset(TKOe, fig.cluster=="4.1")
  c4=subset(TKOe, fig.cluster%in%c("4.0", "4.1"))
  table(c41@assays$RNA@counts["Cbr2", ]>0)
  # TRUE
  # 76
  table(c4@assays$RNA@counts["Gprc5a", ]>0)
  # FALSE  TRUE
  #  59   208
  table(TKOe@assays$RNA@counts["Cbr2", ]>0)
  # FALSE  TRUE
  #  3738  1517
  table(TKOe@assays$RNA@counts["Gprc5a", ]>0)
  # FALSE  TRUE
  #  4277   978
}


### monocle - all epithelial cells
{
  require(monocle)
  # TKOe=readRDS("TKOe.rds")
  temp=TKOe
  # temp=subset(TKOe, subset=(fig.zone=="Cent"))
  gene_metadata=data.frame(gene_short_name=rownames(temp@assays$SCT@counts))
  rownames(gene_metadata)=rownames(temp@assays$SCT@counts)
  TKOe.m <- newCellDataSet(  temp@assays$SCT@counts,
                            phenoData = new("AnnotatedDataFrame", temp@meta.data),
                            featureData = new("AnnotatedDataFrame", gene_metadata),
                            expressionFamily=negbinomial.size() )

  table(TKOe.m@phenoData@data$epitype)
  table(TKOe.m@phenoData@data$celltype)
  TKOe.m <- estimateSizeFactors(TKOe.m)
  TKOe.m <- estimateDispersions(TKOe.m)

  TKOe.m <- detectGenes(TKOe.m, min_expr = 0.1)
  ## only keep expressed genes
  # expressed_genes <- row.names(TKOe.m)[TKOe.m@featureData@data$num_cells_expressed>= 10]
  # TKOe.m <- TKOe.m[expressed_genes,]

  ### use all significant markers of clusters as ordering genes
  tumor.markers=readRDS("TKOe.markers.as.fig.cluster.byHarmony.rds")
  markers=tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[(markers$avg_log2FC)>log2(2), ]
  markers$foldChange=2^(markers$avg_log2FC)
  # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
  markers=markers[order(markers$foldChange, decreasing=TRUE), ]
  ordering.genes=unique(markers$gene)
  # ordering.genes=unique(rownames(markers))
  str(ordering.genes)
  # chr [1:731] "Svs5" "Svs4" "Svs2" "Pate14" "Svs6" "Spinkl" "Clu" "Pate4" ... # FC2
  #  chr [1:1485] "Svs5" "Svs4" "Svs2" "Pate14" "Svs6" "Spinkl" "Clu" ... # FC1.5

  TKOe.m <- setOrderingFilter(TKOe.m,  ordering.genes) # 
  pdf("plot_ordering_genes.pdf")
  plot_ordering_genes(TKOe.m)
  dev.off()

  TKOe.m <- reduceDimension(TKOe.m, max_components = 2,
      method = 'DDRTree')

  TKOe.m <- orderCells(TKOe.m)
  plot_cell_trajectory(TKOe.m, color_by = "fig.cluster")
  dev.off()

  pdf("plot_cell_trajectory_State_TKOe.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_TKOe.m.png")
  plot_cell_trajectory(TKOe.m, color_by = "State", cell_size=0.05, cell_name_size = 8)
  dev.off()

  # set a indicator for time, and sort again
  TKOe.m<-orderCells(TKOe.m, root_state = 2)
  table(rownames(TKOe.m@phenoData@data)==rownames(TKOe@meta.data))
  # TRUE
  # 5255
  TKOe.m@phenoData@data$fig.celltype=TKOe@meta.data$celltype


  pdf("plot_cell_trajectory_byPseudotime_TKOe.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_TKOe.m.png")
  plot_cell_trajectory(TKOe.m, color_by = "Pseudotime", cell_size=0.4)
  dev.off()

  pdf("plot_cell_trajectory_details_fig.cluster_TKOe.m.pdf", width=30, height=6)
  temp=plot_cell_trajectory(TKOe.m, color_by = "fig.cluster") +
      facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()

  pdf("plot_cell_trajectory_fig.cluster_TKOe.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_TKOe.m.png")
  plot_cell_trajectory(TKOe.m, color_by = "epitype", cell_size=0.4)
  dev.off()

  pdf("plot_cell_trajectory_celltype_TKOe.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_TKOe.m.png")
  plot_cell_trajectory(TKOe.m, color_by = "celltype", cell_size=0.05, cell_name_size = 8)
  dev.off()

  pdf("plot_cell_trajectory_fig.celltype_TKOe.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_TKOe.m.png")
  plot_cell_trajectory(TKOe.m, color_by = "fig.celltype", cell_size=0.05, cell_name_size = 8)
  dev.off()

  pdf("plot_cell_trajectory_fig.zone_TKOe.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_TKOe.m.png")
  plot_cell_trajectory(TKOe.m, color_by = "fig.zone", cell_size=0.05, cell_name_size = 8)
  dev.off()

  pdf("plot_cell_trajectory_pheno_TKOe.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_TKOe.m.png")
  plot_cell_trajectory(TKOe.m, color_by = "pheno", cell_size=0.1, cell_name_size = 8)
  dev.off()

  TKOe.m.time=TKOe.m@phenoData@data$Pseudotime
  names(TKOe.m.time)=rownames(TKOe.m@phenoData@data)
  TKOe@meta.data[names(temp.m.time), "pseudotime"]=temp.m.time

  # density plot of cell types' in pseudotime
  {
    # refer to www.r-graph-gallery.com/135-stacked-density-graph.html
    # refer to www.sthda.com/english/wiki/ggplot2-density-plot-quick-start-guide-r-software-and-data-visualization
    library(ggplot2)
    library(hrbrthemes)
    library(dplyr)
    library(tidyr)
    library(viridis)

    # give monocle's pseudotime to seurat object
    TKOe.m.time=TKOe.m@phenoData@data$Pseudotime
    names(TKOe.m.time)=rownames(TKOe.m@phenoData@data)
    TKOe@meta.data[names(TKOe.m.time), "pseudotime"]=TKOe.m.time

    temp.df=TKOe@meta.data[, c("fig.cluster", "pseudotime")]
    # get colors
    require("scales")
    uni.color=levels(TKOe@meta.data$fig.cluster)
    color=hue_pal()(length(uni.color))
    # color=color[length(oorder):1]# 
    # specify Clike's colour
    # color[length(oorder)]="#FF6600"
    pdf("f.3g.CellTypeDensityAlongPseudotime.pdf", width=10, height=5)
    # transparancy desnity plot. alpha: transparent degree
    # {
    #   p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=fig.celltype, fill=fig.celltype)) +
    #     geom_density(adjust=1.5, alpha=.4)
    # }
    # stacked density plot
    # {
      # p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=fig.celltype, fill=fig.celltype)) +
      # geom_density(adjust=1.5, position="fill")
    # }
    {
      p1<- ggplot(temp.df, aes(x=pseudotime, color=fig.cluster)) +
      geom_density() +
      labs(title="Cell type density curve",x="Pseudotime", y = "Density") +
      scale_color_manual(values=color)+
      theme_classic()
    }
    print(p1)
    dev.off()
  }


  # density plot of epitype in pseudotime
  {
    # refer to www.r-graph-gallery.com/135-stacked-density-graph.html
    # refer to www.sthda.com/english/wiki/ggplot2-density-plot-quick-start-guide-r-software-and-data-visualization
    library(ggplot2)
    library(hrbrthemes)
    library(dplyr)
    library(tidyr)
    library(viridis)

    # give monocle's pseudotime to seurat object
    TKOe.m.time=TKOe.m@phenoData@data$Pseudotime
    names(TKOe.m.time)=rownames(TKOe.m@phenoData@data)
    TKOe@meta.data[names(TKOe.m.time), "pseudotime"]=TKOe.m.time

    temp.df=TKOe@meta.data[, c("epitype", "pseudotime")]
    # get colors
    require("scales")
    uni.color=sort(unique(TKOe@meta.data$epitype))
    color=hue_pal()(length(uni.color))
    # color=color[length(oorder):1]# 
    # specify Clike's colour
    # color[length(oorder)]="#FF6600"
    pdf("f.3g.EpitypeDensityAlongPseudotime.pdf", width=10, height=5)
    # transparancy desnity plot. alpha: transparent degree
    # {
    #   p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=fig.celltype, fill=epitype)) +
    #     geom_density(adjust=1.5, alpha=.4)
    # }
    # stacked density plot
    # {
      # p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=fig.celltype, fill=epitype)) +
      # geom_density(adjust=1.5, position="fill")
    # }
    {
      p1<- ggplot(temp.df, aes(x=pseudotime, color=epitype)) +
      geom_density() +
      labs(title="Cell type density curve",x="Pseudotime", y = "Density") +
      scale_color_manual(values=color)+
      theme_classic()
    }
    print(p1)
    dev.off()
  }

  png("plot_cell_trajectory_details_fig.celltype_TKOe.m.png", width=30*100, height=6*100)
  temp=plot_cell_trajectory(TKOe.m, color_by = "fig.celltype") +
      facet_wrap(~celltype, nrow = 1)
  print(temp)
  dev.off()

  pdf("plot_cell_trajectory_details_fig.celltype_TKOe.m.pdf", width=30, height=6)
  temp=plot_cell_trajectory(TKOe.m, color_by = "fig.celltype") +
      facet_wrap(~fig.celltype, nrow = 1)
  print(temp)
  dev.off()


  png("plot_cell_trajectory_details_zone_TKOe.m.png", width=15*100, height=6*100)
  temp=plot_cell_trajectory(TKOe.m, color_by = "fig.zone") +
      facet_wrap(~fig.zone, nrow = 1)
  print(temp)
  dev.off()

  png("plot_cell_trajectory_details_epitype_TKOe.m.png", width=15*100, height=6*100)
  temp=plot_cell_trajectory(TKOe.m, color_by = "epitype") +
      facet_wrap(~epitype, nrow = 1)
  print(temp)
  dev.off()

  # eFig.D.3
  png("plot_cell_trajectory_details_cluster_TKOe.m.png", width=30*100, height=5*100)
  temp=plot_cell_trajectory(TKOe.m, color_by = "fig.cluster", cell_size = 0.8) +
      facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()
  pdf("plot_cell_trajectory_details_cluster_TKOe.m.pdf", width=30, height=5)
  temp=plot_cell_trajectory(TKOe.m, color_by = "fig.cluster", cell_size = 0.4) + facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()

  # saveRDS(TKOe.m, "monocle.TKOe.m.rds")
  # TKOe.m=readRDS("monocle.TKOe.m.rds")

  # find DEGs across branch 
  BEAM_res <- BEAM(TKOe.m, branch_point = 2, cores = 20)
  # saveRDS(BEAM_res, "monocle.BEAM_res.TKOe.m.rds")
  # BEAM_res=readRDS("monocle.BEAM_res.TKOe.m.rds")
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

  # read the subtype markers
  subtypeMarks_csv= read.csv( file = "../collectedBiomarkers_3_subtype.csv",stringsAsFactors=FALSE)
  temp_geneSymbols<- strsplit(subtypeMarks_csv$geneSymbol, ", ")
  subtype_markers=list()
  for(i in 1:nrow(subtypeMarks_csv))
  {
    temp_subtype=subtypeMarks_csv$cellName[i]
    temp_markers=temp_geneSymbols[[i]]
    if(temp_subtype %in% names(subtype_markers))
    {
      subtype_markers[[temp_subtype]]=union(temp_markers,subtype_markers[[temp_subtype]])
    }else 
    {
      subtype_markers[[temp_subtype]]=c(temp_markers)
    }
    subtype_markers[[temp_subtype]]=intersect(subtype_markers[[temp_subtype]], expressed_genes)
  }
  rm(subtypeMarks_csv)

  sig_gene_names=intersect(row.names(subset(BEAM_res, qval<1e-4)), subtype_markers$"TF_fromLieBing")
  sig_gene_names=c(sig_gene_names, "CHGA", "SYP", "ENO2", "KLK3", "AR")
  pdf("heatmap.trajectory.TKOe.m.pdf")
  plot_genes_branched_heatmap(TKOe.m[sig_gene_names,],
                                          branch_point = 2,
                                          num_clusters = 4,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = T)
  dev.off()
  # plot specified genes
  temp_genes <- c("Cbr2", "Gsdmc2", "Mmp7", "F3", "Ceacam1", "Ly6d", "Ascl1", "Neurod1", "Rmst", "Pou2f1", "Yap1")
  pdf("dotplot.trajectory.TKOe.m.pdf")
  plot_genes_branched_pseudotime(TKOe.m[temp_genes,],
                         branch_point = 2,
                         color_by = "celltype",
                         ncol = 1)
  dev.off()


  # find DEGs along Pseudotime 
  diff_test_res.TKOe.m <- differentialGeneTest(TKOe.m, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  # saveRDS(diff_test_res.TKOe.m, "diff_test_res.TKOe.m.rds")
  diff_test_res.TKOe.m=readRDS("diff_test_res.TKOe.m.rds")
  # diff_test_res.bak=diff_test_res
  diff_test_res =  diff_test_res.TKOe.m
  str(diff_test_res[,c("gene_short_name", "pval", "qval")])
  diff_test_res=subset(diff_test_res, qval < 0.1 & pval < 0.01)
  sig_gene_names=intersect(diff_test_res[,"gene_short_name"], subtype_markers$"TF_fromLieBing")
  imp.genes=c()
  imp.genes=c(imp.genes, c("BMPR2", "CCND2", "SOX4"))
  imp.genes=c(imp.genes, c("CDH2", "TWIST2", "HHIP"))
  imp.genes=c(imp.genes, c("COL5A1", "ITGA2", "MATN2", "SEMA3A"))
  imp.genes=c(imp.genes, sig_gene_names)
  imp.gen es=intersect(rownames(TKOe.m), imp.genes)
  pdf("heatmap.trajectory.TKOe.m.pdf")
  plot_pseudotime_heatmap(TKOe.m[imp.genes,],
                  num_clusters = 3,
                  cores = 1,
                  show_rownames = T)
  dev.off()

  pdf("exp.along.time.TKOe.m.pdf")
  imp.genes=c()
  imp.genes=c(imp.genes, c("BMPR2", "CCND2", "SOX4"))
  imp.genes=c(imp.genes, c("CDH2",  "HHIP"))# "TWIST2",
  imp.genes=c(imp.genes, c("COL5A1", "ITGA2", "MATN2", "SEMA3A"))
  genes <- imp.genes
  plot_genes_branched_pseudotime(TKOe.m[genes,],
                       branch_point = 1,
                       color_by = "fig.cluster",
                       ncol = 1)
  dev.off()

}
