### Seurat V3
require(ggplot2)
require(RColorBrewer)
library(BoutrosLab.plotting.general)
library(dplyr)
library(canprot)
source("../functions.seurat.R")
library(Seurat)
library(monocle3)
library(canprot)
library(GEOquery)
packageVersion("monocle3")
# [1] '1.0.0'
packageVersion("Seurat")
# [1] '4.0.1'

# working.directory:
setwd("/ ... ")

### read prim10xary PCa samples
### process from raw counts
data.directory="/ ... /cellrangered"
matrix.directory="/outs/raw_feature_bc_matrix"
dir.raw=list.dirs(path = data.directory, full.names = TRUE, recursive = F)
dir.raw=paste(dir.raw, matrix.directory, sep="")
pos.1=regexpr(data.directory, dir.raw, perl=TRUE)
pos.1=attr(pos.1, "match.length")+2
pos.2=regexpr(matrix.directory, dir.raw, perl=TRUE)-1
sample.names=substring(dir.raw, pos.1, pos.2)
inte.list=list()# the list prepared for further dataset integration
for(i in 1:length(dir.raw))
{
   temp=Read10X(data.dir = dir.raw[i])
   inte.list[i]=CreateSeuratObject(counts = temp, project = sample.names[i], min.cells = 0, min.features = 100)
}
for (i in 1:length(inte.list))
{
  inte.list[[i]][["percent.mt"]] <- PercentageFeatureSet(inte.list[[i]], pattern = "^MT-")
}
# check some qualities
# for (i in 1:length(inte.list))
# {
#   pdf(paste("temp.VlnPlot_",sample.names[i],"_quality.pdf", sep=""), width = 15, height =4)
#   temp=VlnPlot(inte.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size = 0.01, group.by=NULL)
#   temp[[1]]=temp[[1]]+geom_hline(yintercept=300)
#   print(temp)
#   dev.off()
# }
for(i in 1:length(inte.list))
{
  # BD: nFeature_RNA > 1000
  temp<- subset(inte.list[[i]], subset = nFeature_RNA > 200 & percent.mt < 30 ) # & top6k==TRUE)
  print(sample.names[i])
  print(dim(temp))
}
for(i in 1:length(inte.list))
{
  pdf(paste("temp.featureScatter_",sample.names[i],"_quality.pdf", sep=""), width = 15, height =4)
  plot1 <- FeatureScatter(inte.list[[i]], feature1 = "nFeature_RNA", feature2 = "percent.mt")
  plot1=plot1+geom_vline(xintercept=500)+geom_vline(xintercept=1000)+geom_hline(yintercept=60)+geom_hline(yintercept=40)+geom_hline(yintercept=30)
  plot2 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot3 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2 + plot3)
  dev.off()
}
### nFeature_RNA>500 and percent.mt<60 and NOT - restricted on top 6000 cells with highest nFeature_RNA
for(i in 1:length(inte.list))
{
  ### NOT - only consider cells from top 6000 cells
  inte.list[[i]]@meta.data$top6k=rep(TRUE, ncol(inte.list[[i]]))
  temp=sort(inte.list[[i]]@meta.data$nFeature_RNA, decreasing=TRUE)[6000]
  inte.list[[i]]@meta.data$top6k[inte.list[[i]]@meta.data$nFeature_RNA<=temp]=FALSE
  print(table(inte.list[[i]]@meta.data$top6k))
  # inte.list[[i]]<- subset(inte.list[[i]], subset = nFeature_RNA > 500 & percent.mt < 60 ) # & top6k==TRUE) # --->>> this is last version
  inte.list[[i]]<- subset(inte.list[[i]], subset = nFeature_RNA > 200 & percent.mt < 30 & top6k==TRUE)
  # inte.list[[i]]@meta.data$top6k=NULL
}
# filter out cells with > 2 x median nFeature_RNA
for(i in 1:length(inte.list))
{
  inte.list[[i]]@meta.data$top6k=NULL
  temp=median(inte.list[[i]]@meta.data$nFeature_RNA)*2
  inte.list[[i]]@meta.data$double=rep(FALSE, ncol(inte.list[[i]]))
  inte.list[[i]]@meta.data$double[inte.list[[i]]@meta.data$nFeature_RNA>temp]=TRUE
  print(sample.names[i])
  print(table(inte.list[[i]]@meta.data$double))
  print(temp)
  ### NO need to, the nFeature_RNA~nCount_RNA plot show continuous distribution
  # inte.list[[i]]<- subset(inte.list[[i]], subset = (top6k==FALSE) )
  inte.list[[i]]@meta.data$double=NULL
}
# NOT TO individually normalize each batch
# for(i in 1:length(inte.list))
# {
#   inte.list[[i]]=SCTransform(inte.list[[i]], vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE)
# }
prim10x=merge(x = inte.list[[1]], y = inte.list[2:length(inte.list)] )
table(prim10x@meta.data$orig.ident)
#   53   54   55   56   59   62   71   72   73   74   75   76   77
# 2369 1986 3760 5988 3937 5953 5953 5708 3396 5932 3510 5153 5173
### delete patients with too few cells, and 1 metastatic sample 
# prim10x=subset(prim10x, cells = colnames(prim10x)[!prim10x@meta.data$orig.ident%in%c("54", "71", "72")])
prim10x=subset(prim10x, cells = colnames(prim10x)[!prim10x@meta.data$orig.ident%in%c("72")])
table(prim10x@meta.data$orig.ident)
#   53   54   55   56   59   62   71   73   74   75   76   77
# 2369 1986 3760 5988 3937 5953 5953 3396 5932 3510 5153 5173
prim10x@meta.data$fig.sample=as.character(prim10x@meta.data$orig.ident)
prim10x@meta.data$fig.patient=as.character(prim10x@meta.data$orig.ident)
prim10x@meta.data$pheno=rep("prim.10x", nrow(prim10x@meta.data))
prim10x@meta.data$pheno[prim10x@meta.data$fig.patient%in%c("53", "62","71", "75", "77")]="mHSPC.10x"
prim10x@meta.data$pheno[prim10x@meta.data$fig.patient%in%c("73")]="ADT.10x"
prim10x@meta.data$pheno[prim10x@meta.data$fig.patient%in%c("72")]="Meta.10x"
pat2gleason=data.frame(patient=sort(unique(prim10x@meta.data$fig.patient)))
pat2gleason$patient
 # [1] "53" "54" "55" "56" "59" "62" "71" "73" "74" "75" "76" "77"
pat2gleason$gleason=factor(c("5", "5", "5", "5", "3", "3", "5", "5", "2", "5", "3", "5"))
pat2gleason
#    patient gleason
# 1       53       5
# 2       54       5
# 3       55       5
# 4       56       5
# 5       59       3
# 6       62       3
# 7       71       5
# 8       73       5
# 9       74       2
# 10      75       5
# 11      76       3
# 12      77       5
prim10x@meta.data$gleason=pat2gleason$gleason[match(prim10x@meta.data$fig.patient, pat2gleason$patient)]
prim10x@meta.data$gleason=as.character(prim10x@meta.data$gleason)
table(prim10x@meta.data$gleason)
 #    2     3     5
 # 5932 15043 32135
# check ingredients in meta.data
str(prim10x@meta.data)
# 'data.frame':   53110 obs. of  9 variables:
#  $ orig.ident  : chr  "53" "53" "53" "53" ...
#  $ nCount_RNA  : num  17413 371 668 17295 48754 ...
#  $ nFeature_RNA: int  3917 288 334 4462 6889 693 6357 7201 5414 2749 ...
#  $ percent.mt  : num  5.95 4.582 11.677 0.231 6.324 ...
#  $ top6k       : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ fig.sample  : chr  "53" "53" "53" "53" ...
#  $ fig.patient : chr  "53" "53" "53" "53" ...
#  $ pheno       : chr  "mHSPC.10x" "mHSPC.10x" "mHSPC.10x" "mHSPC.10x" ...
#  $ gleason     : chr  "5" "5" "5" "5" ...
prim10x=SCTransform(prim10x, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting

# check most variable genes
str(prim10x@assays$SCT@var.features)
# chr [1:3000] "APOE" "APOC1" "CST1" "TPSB2" "RGS5" "CYTL1" "NPY" "IGFBP7" ...
# top10 <- head(VariableFeatures(prim10x, assay="SCT"), 10)
# # plot variable features with and without labels
# toplot<- VariableFeaturePlot(prim10x, assay = "SCT", selection.method="sctransform")
# toplot <- LabelPoints(plot = toplot, points = top10, repel = TRUE)
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()
# str(prim10x)

DefaultAssay(prim10x)="SCT"
pdf("temp.PCA.pdf")
prim10x <- RunPCA(prim10x, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(prim10x,  ndims = 100)
dev.off()

### PCA-based clustering
{
  # prim10x <- RunUMAP(prim10x, reduction = "pca", dims = 1:25, verbose = FALSE)
  # prim10x <- FindNeighbors(prim10x, reduction = "pca", dims = 1:25, verbose = FALSE)
  # prim10x <- FindClusters(prim10x, verbose = FALSE, resolution = 0.2)
  # table(prim10x@active.ident)

  # prim10x.markers <- FindAllMarkers(prim10x, assay="SCT")
  # saveRDS(prim10x.markers, "prim10x.markers.as.fig.cluster.byPCA.rds")
  # prim10x.markers=readRDS("prim10x.markers.as.fig.cluster.byPCA.rds")
}

### Harmony-based clustering
{
  require(harmony)
  DefaultAssay(prim10x)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  prim10x <- RunHarmony(prim10x, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  dev.off()
  Sys.time()

  # prim10x <- RunUMAP(prim10x, dims = 1:25, verbose = FALSE)
  prim10x <- RunUMAP(prim10x, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prim10x <- FindNeighbors(prim10x, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prim10x <- FindClusters(prim10x, verbose = FALSE, resolution = 0.2)
  table(prim10x@active.ident)
  #    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
  # 8900 7950 2456 1694 1561 1534 1321 1301  844  709  689  499  496  473  361  248
  #   16   17
  #  226  121
  prim10x.markers <- FindAllMarkers(prim10x, assay="SCT")
  saveRDS(prim10x.markers, "prim10x.markers.as.fig.cluster.byHarmony.rds")
  # prim10x.markers =readRDS("prim10x.markers.as.fig.cluster.byHarmony.rds")
}


prim10x@meta.data$fig.cluster=prim10x@meta.data$seurat_clusters
prim10x <- SetIdent(prim10x, value = "fig.cluster")
png("temp_cluster_dimplot.png", width=5*100, height=5*100)
DimPlot(prim10x, label = TRUE, group.by="fig.cluster", pt.size=0.001)# + NoLegend()
dev.off()
png("temp_patient_dimplot.png", width=6*100, height=5*100)
DimPlot(prim10x, label = FALSE, group.by="fig.sample", pt.size=0.001)# + NoLegend()
dev.off()
png("temp_pheno_dimplot.png", width=6*100, height=5*100)
DimPlot(prim10x, label = FALSE, group.by="pheno", pt.size=0.001)# + NoLegend()
dev.off()


### check celltype markers' expression to make sure every cluster's celltype
markers=prim10x.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[markers$avg_log2FC>log2(1.5), ]
setwd("..")
source("functions.seurat.R")
prim10x@meta.data$celltype=celltypeAnno(prim10x, markers)
# [10] "Basal cell: 0.15" CST3+CLU+(Sertoli cell) CD59+C2orf40+(Ciliated cell) 
# [12] "B/Plasma cell_cell201801: -0.08" ACPP+(7FC) KLK2+(1.85FC)
# [14] "Epithelial cell: 0.14 Luminal cell: 0.14" STMN1+ HMMR+
# [15] "Epithelial cell: 0.17 Luminal cell: 0.17" KLK3+(10FC)
# [21] "Basal cell: 0.77 Fibroblast_cell201801: 0.75" COL9A3+

prim10x@meta.data$celltype[prim10x@meta.data$fig.cluster=="9"]="Epithelial cell"
prim10x@meta.data$celltype[prim10x@meta.data$fig.cluster=="11"]="Epithelial cell"
prim10x@meta.data$celltype[prim10x@meta.data$fig.cluster=="20"]="Fibroblast"

setwd("./ ... ")
# no obscure identification
# plot with cluster names
png("temp_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(prim10x, label = TRUE, group.by="celltype") + NoLegend()
dev.off()


# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=prim10x.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.prim10x.cluster","csv",sep="."), row.name=T)


# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # PCa markers
# check.genes=c("FOLH1", "SCHLAP1", "AMACR","CCL2", "CD74", "LY6D", "KLK3", "AR", "ACPP")
# check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

DefaultAssay(prim10x)="RNA"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prim10x, features = check.genes)
dev.off()
DefaultAssay(prim10x)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(prim10x, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(prim10x, "prim10x.rds")
# prim10x=readRDS("prim10x.rds")


### further analysis focused on epithelial cells
table(prim10x@meta.data$celltype)
   # B/Plasma cell       Basal cell Endothelial cell  Epithelial cell
   #           606             1766             6002            28689
   #    Fibroblast       Macrophage        Mast cell    Myofibroblast
   #          5961             2206             1157             2340
   #        T cell
   #          4383
epitypes=c("Epithelial cell", "Luminal cell", "Basal cell")
prim10xe=subset(prim10x, subset=(celltype%in%epitypes))
str(prim10xe)
prim10xe@commands=list()


dim(prim10xe)
# [1] 23913 30455
DefaultAssay(prim10xe)="SCT"
str(prim10xe@assays$SCT@var.features)
# chr chr [1:3000] "NPY" "RGS5" "IGFBP7" "APOE" "APOC1" "CD74" "PLA2G2A" "CST1" ...
# prim10xe <- FindVariableFeatures(prim10xe, assay="SCT", nfeatures = 3000)
# str(prim10xe@assays$SCT@var.features)
# top10 <- head(VariableFeatures(prim10xe, assay="SCT"), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(prim10xe, assay = "SCT")
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# toplot=plot1 + plot2
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()

prim10xe=SCTransform(prim10xe, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
prim10xe <- RunPCA(prim10xe, assay="SCT", verbose = FALSE, npcs = 100)
pdf("temp.PCA.pdf")
ElbowPlot(prim10xe,  ndims = 100)
dev.off()

### clustering basd on PCA
{
  # DefaultAssay(prim10xe)="SCT"
  # prim10xe <- RunUMAP(prim10xe, dims = 1:20, verbose = FALSE)
  # prim10xe <- FindNeighbors(prim10xe, dims = 1:20, verbose = FALSE)

  # DefaultAssay(prim10xe)="SCT" 
  # prim10xe <- FindClusters(prim10xe, verbose = FALSE, resolution = 0.2)
  # table(prim10xe@meta.data$seurat_clusters)
  # prim10xe.markers <- FindAllMarkers(prim10xe, assay="SCT")
  # saveRDS(prim10xe.markers, "prim10xe.markers.as.fig.cluster.byPCA.rds")
  # prim10xe.markers =readRDS("prim10xe.markers.as.fig.cluster.byPCA.rds")

}

### clustering basd on Harmony
{
  require(harmony)
  DefaultAssay(prim10xe)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  prim10xe <- RunHarmony(prim10xe, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
  dev.off()
  Sys.time()

  # prim10x <- RunUMAP(prim10x, dims = 1:25, verbose = FALSE)
  prim10xe <- RunUMAP(prim10xe, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prim10xe <- FindNeighbors(prim10xe, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prim10xe <- FindClusters(prim10xe, verbose = FALSE, resolution = 0.2)
  table(prim10xe@meta.data$seurat_clusters)
  #    0    1    2    3    4    5    6    7    8    9   10   11   12   13
  # 9007 7085 5444 2024 1879 1735 1079  750  533  274  210  177  153  105
  prim10xe.markers <- FindAllMarkers(prim10xe, assay="SCT")
  saveRDS(prim10xe.markers, "prim10xe.markers.as.fig.cluster.byHarmony.rds")
  # prim10xe.markers =readRDS("prim10xe.markers.as.fig.cluster.byHarmony.rds")
}

markers=prim10xe.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.prim10xe.cluster","csv",sep="."), row.name=T)
prim10xe@meta.data$pheno="Primary"
prim10xe@meta.data$pheno[prim10xe@meta.data$fig.patient%in%c("71", "75", "77")]="mHSPC"
prim10xe@meta.data$pheno[prim10xe@meta.data$fig.patient%in%c("73")]="pADT"
prim10xe@meta.data$fig.cluster=prim10xe@meta.data$seurat_clusters
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
png("dimplot.patient.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = FALSE, group.by="fig.patient")
dev.off()
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = FALSE, group.by="gleason")# + NoLegend()
dev.off()
png("dimplot.pheno.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = FALSE, group.by="pheno")# + NoLegend()
dev.off()


# prim10xe=subset(prim10xe, cells = colnames(prim10xe)[!prim10xe@meta.data$fig.cluster%in%c("6", "7", "8", "9", "10")])
table(prim10xe@meta.data$fig.cluster)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13
# 9007 7085 5444 2024 1879 1735 1079  750  533  274  210  177  153  105


# check known epithelial cell types
check.genes=c("KLK3", "KRT8", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # check genes 
# check.genes=c("CCL2", "CD74")
# check.genes=c("SIX3", "MIAT", "SG3", "MIAT", "CHD7", "KIF19", "NEB", "RGS16", "PCP4", "GPX2", "ID4", "CALB1", "SGB2A1")
# check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

prim10xe <- SetIdent(prim10xe, value = "fig.cluster")
DefaultAssay(prim10xe)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prim10xe, features = check.genes)
dev.off()
DefaultAssay(prim10xe)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(prim10xe, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# indi.plot(prim10xe, tocheck="fig.patient")
# indi.plot(prim10xe, tocheck="fig.cluster")

current.cluster.ids <- c(0:13)
new.cluster.ids <- paste(rep("l", 14), c(0:13), sep="")
new.cluster.ids[6]="bI5"
# tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
prim10xe@meta.data$epitype <- plyr::mapvalues(x = as.character(prim10xe@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
png("dimplot.epitype.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = TRUE, group.by="epitype", label.size=8)
dev.off()


### re-cluster basal & intermediate & typeC cells
{
  prim10xb=subset(prim10xe, subset=(fig.cluster%in%c("5")))
  dim(prim10xb)
  # [1] 22543  1833
  DefaultAssay(prim10xb)="SCT"
  str(prim10xb@assays$SCT@var.features)
  # chr [1:3000] "NPY" "CST1" "PLA2G2A" "MSMB" "ADIRF" "KRT17" "MMP7" "MGP" ...
  # prim10xb <- FindVariableFeatures(prim10xb, assay="SCT", nfeatures = 3000)
  # str(prim10xb@assays$SCT@var.features)
  # top10 <- head(VariableFeatures(prim10xb, assay="SCT"), 10)
  # # plot variable features with and without labels
  # plot1 <- VariableFeaturePlot(prim10xb, assay = "SCT")
  # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # toplot=plot1 + plot2
  # png("temp.png", width=10*100, height=5*100)
  # print(toplot)
  # dev.off()

  prim10xb=SCTransform(prim10xb, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  pdf("temp.PCA.pdf")
  prim10xb <- RunPCA(prim10xb, assay="SCT", verbose = FALSE, npcs = 100)
  ElbowPlot(prim10xb,  ndims = 100)
  dev.off()

  ### clustering basd on PCA
  {
    # DefaultAssay(prim10xb)="SCT"
    # prim10xb <- RunUMAP(prim10xb, dims = 1:20, verbose = FALSE)
    # prim10xb <- FindNeighbors(prim10xb, dims = 1:20, verbose = FALSE)

    # DefaultAssay(prim10xb)="SCT" 
    # prim10xb <- FindClusters(prim10xb, verbose = FALSE, resolution = 0.2)
    # table(prim10xb@meta.data$seurat_clusters)
    # prim10xb.markers <- FindAllMarkers(prim10xb, assay="SCT")
    # saveRDS(prim10xb.markers, "prim10xb.markers.as.fig.cluster.byPCA.rds")
    # prim10xb.markers =readRDS("prim10xb.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(prim10xb)
    # [1] "SCT"
    Sys.time()
    pdf("temp.runHarmony.pdf")
    prim10xb <- RunHarmony(prim10xb, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    dev.off()
    Sys.time()

    # prim10x <- RunUMAP(prim10x, dims = 1:25, verbose = FALSE)
    prim10xb <- RunUMAP(prim10xb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    prim10xb <- FindNeighbors(prim10xb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    prim10xb <- FindClusters(prim10xb, verbose = FALSE, resolution = 0.2)
    table(prim10xb@meta.data$seurat_clusters)
    #   0   1   2   3   4   5   6
    # 529 422 417 327  73  39  26
    # prim10xb.markers <- FindAllMarkers(prim10xb, assay="SCT")
    # saveRDS(prim10xb.markers, "prim10xb.markers.as.fig.cluster.byHarmony.rds")
    # prim10xb.markers =readRDS("prim10xb.markers.as.fig.cluster.byHarmony.rds")
  }

  prim10xb@meta.data$fig.cluster=prim10xb@meta.data$seurat_clusters
  png("dimplot.cluster.png", height=5*100, width=5*100)
  DimPlot(prim10xb, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  dev.off()
  png("dimplot.patient.png", height=5*100, width=5*100)
  DimPlot(prim10xb, label = TRUE, group.by="fig.patient")
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  DimPlot(prim10xb, label = TRUE, group.by="gleason")
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  DimPlot(prim10xb, label = TRUE, group.by="pheno")
  dev.off()

  # check known epithelial cell types
  check.genes=c("KLK3", "KRT8", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
  # # prostate/urinary-cancer-specific marker
  # check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
  # # tumor-initiating markers
  # check.genes=c("NKX3-1", "LY6D")
  # # zhangbo's markers
  # check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
  # # check genes 
  # check.genes=c("CCL2", "CD74")
  # check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

  prim10xb <- SetIdent(prim10xb, value = "fig.cluster")
  DefaultAssay(prim10xb)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(prim10xb, features = check.genes)
  dev.off()
  DefaultAssay(prim10xb)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=8*100)
  temp=VlnPlot(prim10xb, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  ### name clusters
  current.cluster.ids <- c(0:5)
  new.cluster.ids <- c("CB50", "B51", "I52", "B53", "I54", "B55")
  # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
  prim10xb@meta.data$epitype <- plyr::mapvalues(x = as.character(prim10xb@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  png("dimplot.epitype.png", height=5*100, width=5*100)
  DimPlot(prim10xb, label = TRUE, group.by="epitype", label.size=8)
  dev.off()

  btype=prim10xb@meta.data$epitype
  names(btype)=rownames(prim10xb@meta.data)
  table(btype)
  # btype
  #  b51  b53  b55 cb50  i52  i54
  #  353  222  102  525  343  190

  # saveRDS(prim10xb, "prim10xb.rds")
  # prim10xb=readRDS("prim10xb.rds")
}


### re-cluster basal & typeC cells
{
  prim10xc=subset(prim10xb, subset=(fig.cluster%in%c("0")))
  dim(prim10xb)
  # [1] 22543  1833
  DefaultAssay(prim10xc)="SCT"
  str(prim10xc@assays$SCT@var.features)
  # chr [1:3000] "NPY" "CST1" "PLA2G2A" "MSMB" "ADIRF" "KRT17" "MMP7" "MGP" ...
  # prim10xc <- FindVariableFeatures(prim10xc, assay="SCT", nfeatures = 3000)
  # str(prim10xc@assays$SCT@var.features)
  # top10 <- head(VariableFeatures(prim10xc, assay="SCT"), 10)
  # # plot variable features with and without labels
  # plot1 <- VariableFeaturePlot(prim10xc, assay = "SCT")
  # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # toplot=plot1 + plot2
  # png("temp.png", width=10*100, height=5*100)
  # print(toplot)
  # dev.off()

  prim10xc=SCTransform(prim10xc, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  pdf("temp.PCA.pdf")
  prim10xc <- RunPCA(prim10xc, assay="SCT", verbose = FALSE, npcs = 100)
  ElbowPlot(prim10xc,  ndims = 100)
  dev.off()

  ### clustering basd on PCA
  {
    # DefaultAssay(prim10xc)="SCT"
    # prim10xc <- RunUMAP(prim10xc, dims = 1:20, verbose = FALSE)
    # prim10xc <- FindNeighbors(prim10xc, dims = 1:20, verbose = FALSE)

    # DefaultAssay(prim10xc)="SCT" 
    # prim10xc <- FindClusters(prim10xc, verbose = FALSE, resolution = 0.2)
    # table(prim10xc@meta.data$seurat_clusters)
    # prim10xc.markers <- FindAllMarkers(prim10xc, assay="SCT")
    # saveRDS(prim10xc.markers, "prim10xc.markers.as.fig.cluster.byPCA.rds")
    # prim10xc.markers =readRDS("prim10xc.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(prim10xc)
    # [1] "SCT"
    Sys.time()
    pdf("temp.runHarmony.pdf")
    prim10xc <- RunHarmony(prim10xc, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    dev.off()
    Sys.time()

    # prim10x <- RunUMAP(prim10x, dims = 1:25, verbose = FALSE)
    # prim10xc <- RunUMAP(prim10xc, reduction = "harmony", dims = 1:25, verbose = FALSE)
    prim10xc <- FindNeighbors(prim10xc, reduction = "harmony", dims = 1:25, verbose = FALSE)
    prim10xc <- FindClusters(prim10xc, verbose = FALSE, resolution = 0.2)
    table(prim10xc@meta.data$seurat_clusters)
    #   0   1
    # 279 246
    # prim10xc.markers <- FindAllMarkers(prim10xc, assay="SCT")
    # saveRDS(prim10xc.markers, "prim10xc.markers.as.fig.cluster.byHarmony.rds")
    # prim10xc.markers =readRDS("prim10xc.markers.as.fig.cluster.byHarmony.rds")
  }

  prim10xc@meta.data$fig.cluster=prim10xc@meta.data$seurat_clusters
  png("dimplot.cluster.png", height=5*100, width=5*100)
  DimPlot(prim10xc, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  dev.off()
  png("dimplot.patient.png", height=5*100, width=5*100)
  DimPlot(prim10xc, label = TRUE, group.by="fig.patient")
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  DimPlot(prim10xc, label = TRUE, group.by="gleason")
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  DimPlot(prim10xc, label = TRUE, group.by="pheno")
  dev.off()

  # check known epithelial cell types
  check.genes=c("KLK3", "KRT8", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
  # # prostate/urinary-cancer-specific marker
  # check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
  # # tumor-initiating markers
  # check.genes=c("NKX3-1", "LY6D")
  # # zhangbo's markers
  # check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
  # # check genes 
  # check.genes=c("CCL2", "CD74")
  # check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

  prim10xc <- SetIdent(prim10xc, value = "fig.cluster")
  DefaultAssay(prim10xc)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(prim10xc, features = check.genes)
  dev.off()
  DefaultAssay(prim10xc)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=8*100)
  temp=VlnPlot(prim10xc, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  ### name clusters
  current.cluster.ids <- c(0:1)
  new.cluster.ids <- c("C500", "B501")
  # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
  prim10xc@meta.data$epitype <- plyr::mapvalues(x = as.character(prim10xc@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  png("dimplot.epitype.png", height=5*100, width=5*100)
  DimPlot(prim10xc, label = TRUE, group.by="epitype", label.size=8)
  dev.off()

  ctype=prim10xc@meta.data$epitype
  names(ctype)=rownames(prim10xc@meta.data)
  table(ctype)
  # ctype
  # B501 C500
  #  246  279

  # saveRDS(prim10xc, "prim10xc.rds")
  # prim10xc=readRDS("prim10xc.rds")
}


### epitype annotation
current.cluster.ids <- c(0:13)
new.cluster.ids <- paste("L", c(0:13), sep="")

# tL: tumor luminal, tB: tumor basal, tCy: tumor CellCycle
prim10xe@meta.data$epitype <- plyr::mapvalues(x = as.character(prim10xe@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot
prim10xe@meta.data[names(btype), "epitype"]=btype
prim10xe@meta.data[names(ctype), "epitype"]=ctype
pdf("dimplot.prim10xe.epitype.pdf", height=5, width=5)
DimPlot(prim10xe, label = FALSE, group.by="epitype", label.size=5, raster=TRUE)
dev.off()



### celltype annotation # add refined basal/typeC/interm types into epitype slot
current.cluster.ids <- unique(prim10xe@meta.data$epitype)
current.cluster.ids=sort(current.cluster.ids)
current.cluster.ids
#  [1] "b501" "b51"  "b53"  "b55"  "c500" "i52"  "i54"  "L0"   "L1"   "L10"
# [11] "L11"  "L12"  "L13"  "L2"   "L3"   "L4"   "L6"   "L7"   "L8"   "L9"
new.cluster.ids <- rep("Luminal", length(current.cluster.ids))
new.cluster.ids[grepl("B", current.cluster.ids)]="Basal"
new.cluster.ids[grepl("C", current.cluster.ids)]="TypeC"
new.cluster.ids[grepl("I", current.cluster.ids)]="Intermediate"
prim10xe@meta.data$celltype <- plyr::mapvalues(x = as.character(prim10xe@meta.data$epitype), from = current.cluster.ids, to = new.cluster.ids)
pdf("dimplot.prim10xe.celltype.pdf", height=5, width=5)
DimPlot(prim10xe, label = FALSE, group.by="celltype", label.size=5, raster=TRUE)
dev.off()


# check known epithelial cell types
check.genes=c("KLK3", "KRT8", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # check genes 
# check.genes=c("CCL2", "CD74")
# check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

prim10xe <- SetIdent(prim10xe, value = "celltype")
DefaultAssay(prim10xe)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prim10xe, features = check.genes)
dev.off()
DefaultAssay(prim10xe)="SCT"
png("temp.VlnPlot.png", width=10*100, height=12*100)
temp=VlnPlot(prim10xe, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


### barplot the percentage of each samples in cell types
{
  # give sample ID as phenotype + number
  table(prim10xe@meta.data$pheno)
  # mHSPC    pADT Primary
  #  8982    1322   20151
  prim10xe@meta.data$fig.pheno=factor(prim10xe@meta.data$pheno, levels=c("Primary", "mHSPC", "pADT"))
  pheno.sample=prim10xe@meta.data[, c("fig.pheno", "fig.sample")]
  pheno.sample$fig.sample=as.character(pheno.sample$fig.sample)
  # 'data.frame':   39463 obs. of  2 variables:
  #  $ fig.pheno : Factor w/ 8 levels "Normal prostate",..: 1 1 1 1 1 1 1 1 1 1 ...
  #  $ fig.sample: chr  "waizhouyouCBH" "waizhouyouCBH" "waizhouyouCBH" "waizhouyouCBH" ...
  pheno.sample=pheno.sample[order(pheno.sample$fig.pheno, pheno.sample$fig.sample, decreasing=c(FALSE, FALSE)), ]
  unique.sample=unique(pheno.sample$fig.sample)
  pheno.sample$fig.sample=factor(pheno.sample$fig.sample, levels=unique.sample)
  pheno.sample$fig.sample.id=as.numeric(pheno.sample$fig.sample)
  pheno.sample$pheno.sample=paste(as.character(pheno.sample$fig.pheno), pheno.sample$fig.sample.id)
  unique.pheno.sample=unique(pheno.sample$pheno.sample)
  pheno.sample$pheno.sample=factor(pheno.sample$pheno.sample, levels=unique.pheno.sample)
  prim10xe@meta.data$pheno.sample=pheno.sample[rownames(prim10xe@meta.data), ]$pheno.sample

  celltype.sta=table(prim10xe@meta.data$pheno.sample, prim10xe@meta.data$celltype)
  celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
  # celltype.sta=celltype.sta[ c("Clike", "TypeC", "Basal", "NE"), !grepl("incidental|Normal", colnames(celltype.sta), perl=TRUE) ]
  celltype.sta=celltype.sta[ c("TypeC", "Basal", "Intermediate"),]
  barN=ncol(celltype.sta)
  n_top = nrow(celltype.sta) 
  pdf(paste("prim10xe.barplot.celltypeInSample.pdf",sep=""))
  # the width between paper edge and plot ( down,left,up,right : four directions)
  par(mar = c(7, 4, 1, 8), xpd=T)
  # color=DiscretePalette(n_top, palette = "polychrome")
  require("scales")
  color=hue_pal()(n_top)
  barplot( celltype.sta, xlab="", ylab="Percentage", col=color, names.arg=colnames(celltype.sta), las=2, cex.names=0.5) # , xaxt="n"
  # facts=barN+1
  end_point = barN #0.5 + n_top *facts-1
  # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
  legend( barN+2,60, rownames(celltype.sta), fill=color );
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
  dev.off()
}


# saveRDS(prim10xe, "prim10xe.rds")
# prim10xe=readRDS("prim10xe.rds")





### Not employed --->>>
  ### CCA to identify epithelial cell types WithiN luminal tumor cells
  # norme=readRDS("norme.rds")
  prim10xl=subset(prim10xe, subset=celltype=="Luminal")
  anchors <- FindTransferAnchors(reference = norme, query = prim10xl, 
      dims = 1:30, normalization.method="LogNormalize", reference.assay="SCT")
  predictions <- TransferData(anchorset = anchors, refdata = norme$celltype, 
      dims = 1:30)
  prim10xl@meta.data$predtype=predictions$predicted.id 
  png("dimplot.predtype.png", height=5*100, width=5*100)
  DimPlot(prim10xl, label = TRUE, group.by="predtype")
  dev.off()
### Not employed <<<---


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(prim10xe@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(prim10xe@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(prim10xe@meta.data$fig.patient, prim10xe@meta.data$celltype)
temp=apply(pct.matrix, 1, sum)
pct.matrix=round(apply(pct.matrix, 2, function(x) x/temp), 2)
temp.name=rownames(pct.matrix)
temp.name
#  [1] "CZK-T"      "DHB-T"      "HYQ"        "LPX"        "PXL"
#  [6] "QLX-M"      "R-M"        "RSC1023"    "S6-TUMOR-M" "SHSCR"
# [11] "SJ-M"       "T"          "ZQF-M"
### 9.5: N1.Gleason5
GS=c( NA, NA, 9.5, 5, 9.5, 9.5, 2, 9.4, 5, 5, 2, NA, 3)
pct.matrix=pct.matrix[order(GS, decreasing=TRUE), ]
write.csv(pct.matrix, file="part4__prim10xe__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(prim10xe@meta.data$epitype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(prim10xe@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(prim10xe@meta.data$fig.patient, prim10xe@meta.data$epitype)
temp=apply(pct.matrix, 1, sum)
pct.matrix=round(apply(pct.matrix, 2, function(x) x/temp), 2)
temp.name=rownames(pct.matrix)
temp.name
#  [1] "CZK-T"      "DHB-T"      "HYQ"        "LPX"        "PXL"
#  [6] "QLX-M"      "R-M"        "RSC1023"    "S6-TUMOR-M" "SHSCR"
# [11] "SJ-M"       "T"          "ZQF-M"
### 9.5: N1.Gleason5
GS=c( NA, NA, 9.5, 5, 9.5, 9.5, 2, 9.4, 5, 5, 2, NA, 3)
pct.matrix=pct.matrix[order(GS, decreasing=TRUE), ]
write.csv(pct.matrix, file="part4__prim10xe__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


# saveRDS(prim10xe, "prim10xe.rds")
# prim10xe=readRDS("prim10xe.rds")


### celltype-specific DEG identification
table(prim10xe@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   12893     925
DefaultAssay(prim10xe)="SCT"
### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(prim10xe@meta.data$celltype)
celltype=sort(celltype)
prim10xe = SetIdent(prim10xe, value = "celltype")
for(i in celltype)
{
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        markerlist[[temp.name]]=FindMarkers(prim10xe, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
      }
  }
}
# saveRDS(markerlist, "prim10xe.celltype.pairwise.rds")
# markerlist=readRDS("prim10xe.celltype.pairwise.rds")

sp.markerlist=list()
celltype=unique(prim10xe@meta.data$celltype)
celltype=sort(celltype)
for(i in celltype)
{
  temp.genes=rownames(prim10xe@assays$SCT@data)
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        cur.list=markerlist[[temp.name]]
        cur.list=cur.list[order(cur.list$avg_log2FC, decreasing=TRUE), ]
        cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(2) & cur.list$pct.1>0.3 & cur.list$pct.2<0.1,])
        temp.genes=intersect(temp.genes, cur.genes)
      }
  }
  sp.markerlist[[i]]=temp.genes
}
str(sp.markerlist)
sp.markerlist[["Luminal"]]=c(sp.markerlist[["Luminal"]], "KLK3", "AR")
sp.markerlist[["TypeC"]]=c(sp.markerlist[["TypeC"]], "KRT4", "TACSTD2", "PSCA")

temp=subset(prim10xe, downsample=300)
sp.markerlist=sp.markerlist[c("Basal", "Intermediate", "Luminal", "TypeC")]
tocheck=Reduce(union, sp.markerlist)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("Basal", "Intermediate", "Luminal", "TypeC")
pdf(paste("heatmap.prim10xe.typeC.pdf", sep=""),  width=10, height=length(tocheck)/3)
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

table(normale@assays$RNA@counts["FABP4", ]==0)
# FALSE  TRUE
#     4 11941

best=c("CD55", "IL1RN", "IL1A", "UPK1A", "UPK1B", "PIM1", "PMAIP1")
tocheck=best

# created by part11__mergeEpi.V3.excludeIncidente.R
tocheck=epi.tumor.clike.markers


### check specific marker across all cell types
# prim10x=readRDS("prim10x.rds")
check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
clike.cells=rownames(prim10xe@meta.data)[prim10xe@meta.data$celltype=="Clike"]
prim10x@meta.data$celltype[rownames(prim10x@meta.data)%in%clike.cells]="Clike"
prim10x <- SetIdent(prim10x, value = "celltype")
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prim10x, features = check.genes)
dev.off()
DefaultAssay(prim10x)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(prim10x, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

prim10x <- SetIdent(prim10x, value = "celltype")
clike.allcelltype=FindMarkers(prim10x, ident.1 = "Clike", assay="SCT")
# saveRDS(clike.allcelltype, "markers.clikeVSallcelltypes.rds")
# clike.allcelltype=readRDS("markers.clikeVSallcelltypes.rds")


### DEG identification, across all clusters
prim10xe <- SetIdent(prim10xe, value = "celltype")
prim10xe.markers <- FindAllMarkers(prim10xe, slot="data", assay="SCT")
# saveRDS(prim10xe.markers, "prim10xe.markers.celltype.rds")
# prim10xe.markers=readRDS("prim10xe.markers.celltype.rds")

### plot heatmap
# prim10xe.markers=readRDS("prim10xe.indivSCT.markers.as.fig.cluster.rds")
# prim10xe <- NormalizeData(prim10xe, normalization.method = "LogNormalize", scale.factor = 10000, assay="SCT" )
# prim10xe <- ScaleData(prim10xe, features = rownames(prim10xe), assay="RNA")

# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=prim10xe.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.prim10xe.cluster","csv",sep="."), row.name=T)


# heatmap top 10 DEGs for every cluster
markers=prim10xe.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(2) & markers$pct.1>0.3 & markers$pct.2<0.1, ]
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
top10$cluster=factor(top10$cluster, levels=c("Clike", "TypeC", "Luminal", "Basal", "CellCycle"))
top10=top10[order(top10$cluster), ]
pdf(paste("heatmap.prim10xe.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=nrow(top10)/8)
temp=subset(prim10xe, downsample=300)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
tocheck=c(top10$gene, "KLK3", "AR")
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


# check cell-surface markers of each celltype
markers=prim10xe.markers
markers$cluster=factor(markers$cluster, levels=c("Clike", "TypeC", "Luminal", "Basal", "CellCycle"))
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, abs(markers$avg_log2FC),decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[ (abs(markers$avg_log2FC)>log2(2))  , ]
require(msigdbr)
cc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
cc = cc %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
cdgenes=sort(unique(cc$gene_symbol[cc$gs_name=="GO_CELL_SURFACE"]))
tempcheck=cdgenes
temp.markers=markers[markers$gene%in%tempcheck, ]
top10 <- temp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = (avg_log2FC) )
if(nrow(top10)>0)
{
  DefaultAssay(prim10xe)="SCT"
  pdf(paste("temp.heatmap.CD_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
  temp=subset(prim10xe, downsample=1000)
  levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
  temp.heatmap<-DoHeatmap(temp , features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
  DefaultAssay(prim10xe)="SCT"
}


### pairwise DEGs to find three surface marker for Clike cells
markerlist=readRDS("prim10xe.celltype.pairwise.rds")
require(msigdbr)
cc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
cc = cc %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
cdgenes=sort(unique(cc$gene_symbol[cc$gs_name=="GO_CELL_SURFACE"]))
tempcheck=cdgenes
tocheck=c("Clike__TypeC", "Clike__Luminal", "Clike__Basal")
temp.marker=c()
for(i in tocheck)
{
  markers=markerlist[[i]]
  markers=markers[rownames(markers)%in%cdgenes, ]
  # markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[order(markers$avg_log2FC,decreasing = TRUE), ]
  markers$foldchange=2^(markers$avg_log2FC)
  markers=markers[  markers$pct.1>0.5 & markers$pct.2<0.2, ]
  print(i)
  print(markers[1:5,])
  temp.marker=c(temp.marker, rownames(markers[1:5,]))
}
#                p_val avg_log2FC pct.1 pct.2     p_val_adj foldchange
# [1] "Clike__TypeC"
# FGFR3  4.312549e-61  0.8725208 0.585 0.086 8.653562e-57   1.830859
# TGFBR3 3.913872e-37  0.7205307 0.591 0.159 7.853575e-33   1.647788
# [1] "Clike__Luminal"
# UPK1A  0.000000e+00   4.588778 0.895 0.025  0.000000e+00  24.063562
# CD55  2.705167e-158   3.875505 0.877 0.194 5.428189e-154  14.677200
# BMP2   0.000000e+00   2.967132 0.784 0.020  0.000000e+00   7.819803
# [1] "Clike__Basal"
# UPK1A   4.265789e-99  3.8899713 0.895 0.149 8.559732e-95  14.825114
# WNT5A   1.843048e-44  1.0325982 0.626 0.129 3.698259e-40   2.045705
# SRPX2   9.797940e-41  0.6974380 0.561 0.107 1.966055e-36   1.621622
temp.marker=c("FGFR3", "TGFBR3", "UPK1A", "CD55", "SRPX2", "WNT5A")
temp.marker=c("TGFBR3", "UPK1A", "SRPX2")
temp.marker=c("UPK1A", "CD55", "TGFBR3")
temp=subset(prim10xe, downsample=1000)
# temp=NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA" )
# temp=ScaleData(temp, features = rownames(temp), assay="RNA" )
levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
pdf(paste("temp.heatmap.pdf", sep=""), width=(length(unique(prim10xe@active.ident))/4+5), height=length(temp.marker)/10+3)
temp.heatmap<-DoHeatmap(temp , features = temp.marker, slot="scale.data", assay="SCT", angle = 0, size=3) 
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### plot heatmap for special categories of DEGs
setwd("..")
subtype_markers=subtypeMarker(prim10xe)
setwd("./ ... ")
prim10xe.markers=readRDS("prim10xe.markers.celltype.rds")
markers=prim10xe.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
# for TFs, change log2(1.5) to log(1.2)
# prim10xe <- ScaleData(prim10xe, features = rownames(prim10xe), assay="RNA" )
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
   top10$cluster=factor(top10$cluster, levels=c("Clike", "TypeC", "Luminal", "Basal", "CellCycle"))
  top10=top10[order(top10$cluster), ]
  if(nrow(top10)>0)
  {
    DefaultAssay(prim10xe)="SCT"
    pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
    temp=subset(prim10xe, downsample=300)
    temp <- SetIdent(temp, value = "celltype")
    levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
    temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()
    DefaultAssay(prim10xe)="SCT"
  }
}

# for too few cluster DExpressing genes, do vlnplot
tocheck=c("CCL2", "CXCL11", "CXCL1", "CXCL8", "CXCL17", "AR")
pdf("VlnPlot.chemokines_celltype.pdf", width = 10, height =10)
temp=VlnPlot(prim10xe, features = tocheck, pt.size = 0, group.by="celltype", ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# enrichment of Clike compared with all other cells
markers=prim10xe.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)

# enrichment of Clike compared with TypeC
markerlist=readRDS("prim10xe.celltype.pairwise.rds")
markers=markerlist[["Clike__TypeC"]]
markers=markers[order(markers$avg_log2FC,decreasing = TRUE), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(rownames(markers), topNterm=20)


# Qusage to enrichment
require(qusage)
require(msigdbr)
m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
m_t2g.reactome = m_df.reactome %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.kegg = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
m_t2g.kegg = m_df.kegg %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.hall = msigdbr(species = "Homo sapiens", category = "H")
m_t2g.hall = m_df.hall %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.regulate = msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
m_t2g.regulate = m_df.regulate %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
# convert annotation from data.frame to list format
kegg.list=list()
kegg.name=unique(m_t2g.kegg$gs_name)
for(i in kegg.name)
{
  kegg.list[[i]]=m_t2g.kegg$gene_symbol[m_t2g.kegg$gs_name==i]
}
kegg=read.gmt("../KEGG_human_20190613_fromZhangbo.gmt")
hall.list=list()
hall.name=unique(m_t2g.hall$gs_name)
for(i in hall.name)
{
  hall.list[[i]]=m_t2g.hall$gene_symbol[m_t2g.hall$gs_name==i]
}
packageVersion("qusage")
# [1] '2.24.0'
require(qusage)
  #seurat.all = a seurat object from scran
  #nm = name of the run
  #gs = gene set to test, a list
# the "run_qusage_heatmap" function needs that @meta.data$for_temp_color is a numeric format vector with MIN=1 rather than MIN=0, whic represent the cluster.ID
prim10xe@meta.data$for_temp_color=as.numeric(factor(prim10xe@meta.data$celltype, levels=c("Clike", "TypeC", "Luminal", "Basal")))
run_qusage_heatmap.seurat3(prim10xe, nm = 'prim10xe.celltype', kegg, my.seed=100)
run_qusage_heatmap.seurat3(prim10xe, nm = 'prim10xe.celltype', hall.list, my.seed=100)


### NOT RUN yet --->>>
### individually analysis of basal cells
bas13=subset(prim10xe, subset=(celltype=="tBI"))
bas13=ScaleData(bas13, features = rownames(bas13), assay="RNA" )

DefaultAssay(bas13)="SCT"
bas13 <- FindVariableFeatures(bas13, assay="SCT", nfeatures = 3000)
str(bas13@assays$SCT@var.features)
# chr [1:3000] "LTF" "MMP10" "KLK3" "PLA2G2A" "CCL2" "OGN" "SCGB1A1" "KRT13" ...
top10 <- head(VariableFeatures(bas13, assay="SCT"), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bas13, assay = "SCT")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
toplot=plot1 + plot2
png("temp.png", width=10*100, height=5*100)
print(toplot)
dev.off()

bas13 <- RunPCA(bas13, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(bas13,  ndims = 100)
dev.off()

DefaultAssay(bas13)="SCT"
bas13 <- RunUMAP(bas13, dims = 1:9, verbose = FALSE)
bas13 <- FindNeighbors(bas13, dims = 1:9, verbose = FALSE)

DefaultAssay(bas13)="SCT" 
bas13 <- FindClusters(bas13, verbose = FALSE, resolution = 0.5)
table(bas13@meta.data$seurat_clusters)
#   0   1   2   3
# 373 217 166 152
bas13@meta.data$fig.cluster=bas13@meta.data$seurat_clusters
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(bas13, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
png("dimplot.patient.png", height=5*100, width=5*100)
DimPlot(bas13, label = TRUE, group.by="fig.patient")
dev.off()
png("dimplot.zone.png", height=5*100, width=5*100)
DimPlot(bas13, label = TRUE, group.by="fig.zone")
dev.off()

# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # check genes 
# check.genes=c("CCL2", "CD74")

DefaultAssay(bas13)="RNA"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(bas13, features = check.genes)
dev.off()
DefaultAssay(bas13)="SCT"
png("temp.VlnPlot.png", width=10*100, height=8*100)
temp=VlnPlot(bas13, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2, assay = "RNA")
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95, assay="RNA")
}
print(temp)
dev.off()

# check every cluster's and patient's cells
indi.plot(bas13, tocheck="fig.patient")
indi.plot(bas13, tocheck="fig.cluster")

### DEG identification
bas13 <- SetIdent(bas13, value = "fig.cluster")
bas13.markers <- FindAllMarkers(bas13, assay="SCT")
# saveRDS(bas13.markers, "bas13.indivSCT.markers.as.fig.cluster.rds")
# bas13.markers=readRDS("bas13.indivSCT.markers.as.fig.cluster.rds")

### plot heatmap
# bas13.markers=readRDS("bas13.indivSCT.markers.as.fig.cluster.rds")
# bas13 <- ScaleData(bas13, features = rownames(bas13), assay="RNA" )
# scale.data.bak=bas13@assays$SCT@scale.data
# bas13 <- ScaleData(bas13, features = rownames(bas13), assay="SCT" )
# bas13@assays$SCT@scale.data=scale.data.bak

# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=bas13.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.bas13.cluster","csv",sep="."), row.name=T)

# heatmap top 10 DEGs for every cluster
markers=bas13.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
pdf(paste("heatmap.bas13.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=nrow(top10)/8)
temp.heatmap<-DoHeatmap( subset(bas13, downsample=300), features = top10$gene, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

# heatmap KRT genes for every cluster
markers=bas13.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log(1.1) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
krt.iden=grepl("^KRT", markers$gene, perl=TRUE)
features=markers$gene[krt.iden]
features=c("KLK3","AR", "KRT13", "CCL2", features)
pdf(paste("heatmap.bas13.cluster.KRT.pdf", sep=""), width=(length(unique(markers$cluster))/4+5), height=length(features)/8)
temp.heatmap<-DoHeatmap( subset(bas13, downsample=300), features = features, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

### name clusters
current.cluster.ids <- c(0:5)
new.cluster.ids <- c("tBI0", "tB1", "tB2", "tBC3", "tB4", "tBI5")
# tBI: tumor intermediate cells in basal cluster, tBC: tumor type-C cells in basal cluster, tB: tumor basal cells
bas13@meta.data$celltype <- plyr::mapvalues(x = as.character(bas13@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)

# save obj
bas13@assays$RNA@scale.data=matrix(0, 0, 0)
bas13@commands=list()
bas13@meta.data$basal.sig=NULL
bas13@meta.data$cycle.sig=NULL
bas13@meta.data$SCT_snn_res.0.5=NULL
bas13@meta.data$SCT_snn_res.0.4=NULL
bas13@meta.data$SCT_snn_res.0.3=NULL
bas13@meta.data$SCT_snn_res.0.2=NULL
bas13@meta.data$SCT_snn_res.0.1=NULL
table(bas13@meta.data$celltype)
 # tB1  tB2  tB4 tBC3 tBI0 tBI5
 # 213  165  109  115  219   87
# saveRDS(bas13, "bas13.seurat3.indivSCTed.celltyped.rds")
# bas13=readRDS("bas13.seurat3.indivSCTed.celltyped.rds")


### individually analysis of basal cells
lum13=subset(prim10xe, subset=(celltype%in%c("tLum", "tInt")))
lum13=ScaleData(lum13, features = rownames(lum13), assay="RNA" )

DefaultAssay(lum13)="SCT"
lum13 <- FindVariableFeatures(lum13, assay="SCT", nfeatures = 3000)
str(lum13@assays$SCT@var.features)
# chr [1:3000] "CST1" "ACY3" "CHGB" "LYZ" "AL031708.1" "NPY" "AREG" "FN1" ...
top10 <- head(VariableFeatures(lum13, assay="SCT"), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(lum13, assay = "SCT")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
toplot=plot1 + plot2
png("temp.png", width=10*100, height=5*100)
print(toplot)
dev.off()

lum13 <- RunPCA(lum13, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(lum13,  ndims = 100)
dev.off()

DefaultAssay(lum13)="SCT"
lum13 <- RunUMAP(lum13, dims = 1:21, verbose = FALSE)
lum13 <- FindNeighbors(lum13, dims = 1:21, verbose = FALSE)

DefaultAssay(lum13)="SCT" 
lum13 <- FindClusters(lum13, verbose = FALSE, resolution = 0.1)
table(lum13@meta.data$seurat_clusters)
#   0   1   2   3
# 373 217 166 152
lum13@meta.data$fig.cluster=lum13@meta.data$seurat_clusters
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(lum13, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
png("dimplot.patient.png", height=5*100, width=5*100)
DimPlot(lum13, label = TRUE, group.by="fig.patient")
dev.off()
png("dimplot.zone.png", height=5*100, width=5*100)
DimPlot(lum13, label = TRUE, group.by="fig.zone")
dev.off()
png("dimplot.celltype.png", height=5*100, width=5*100)
DimPlot(lum13, label = TRUE, group.by="celltype")
dev.off()

# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # check genes 
# check.genes=c("CCL2", "CD74")

DefaultAssay(lum13)="RNA"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(lum13, features = check.genes)
dev.off()
DefaultAssay(lum13)="SCT"
png("temp.VlnPlot.png", width=10*100, height=8*100)
temp=VlnPlot(lum13, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2, assay = "RNA")
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95, assay="RNA")
}
print(temp)
dev.off()

# check every cluster's and patient's cells
indi.plot(lum13, tocheck="fig.patient")
indi.plot(lum13, tocheck="fig.cluster")

### DEG identification
lum13 <- SetIdent(lum13, value = "fig.cluster")
lum13.markers <- FindAllMarkers(lum13, assay="SCT")
# saveRDS(lum13.markers, "lum13.indivSCT.markers.as.fig.cluster.rds")
# lum13.markers=readRDS("lum13.indivSCT.markers.as.fig.cluster.rds")

### plot heatmap
# lum13.markers=readRDS("lum13.indivSCT.markers.as.fig.cluster.rds")
# lum13 <- ScaleData(lum13, features = rownames(lum13), assay="RNA" )
# scale.data.bak=lum13@assays$SCT@scale.data
# lum13 <- ScaleData(lum13, features = rownames(lum13), assay="SCT" )
# lum13@assays$SCT@scale.data=scale.data.bak

# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=lum13.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.lum13.cluster","csv",sep="."), row.name=T)

# heatmap top 10 DEGs for every cluster
markers=lum13.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = (avg_log2FC) )
pdf(paste("heatmap.lum13.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=nrow(top10)/8)
temp.heatmap<-DoHeatmap( subset(lum13, downsample=300), features = top10$gene, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

# heatmap KRT genes for every cluster
markers=lum13.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log(1.1) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
krt.iden=grepl("^KRT", markers$gene, perl=TRUE)
features=markers$gene[krt.iden]
features=c("KLK3","AR", "KRT13", "CCL2", features)
pdf(paste("heatmap.lum13.cluster.KRT.pdf", sep=""), width=(length(unique(markers$cluster))/4+5), height=length(features)/8)
temp.heatmap<-DoHeatmap( subset(lum13, downsample=300), features = features, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(lum13@meta.data$fig.cluster)
#     0     1     2     3
# 11931  6763  4929   446
pct.matrix=table(lum13@meta.data$fig.patient, lum13@meta.data$fig.cluster)
temp=apply(pct.matrix, 1, sum)
pct.matrix=round(apply(pct.matrix, 2, function(x) x/temp), 2)
temp.name=rownames(pct.matrix)
# temp.name
# [1] "53" "54" "55" "56" "59" "62" "71" "73" "74" "75" "76" "77"
GS=c(9, 10, 9, 9, 7.2, 7.2, 9, 9, 7.1, 9, 7.2, 9)
pct.matrix=pct.matrix[order(GS, decreasing=TRUE), ]
write.csv(pct.matrix, file="cell.ori.seurat3.merge__PCa13.samples.diffGS_part3________201217__plotCellPct__luminal.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:8, "1"]), as.numeric(pct.matrix[9:12, "1"]))
str(res)


# Qusage to enrichment
require(qusage)
require(msigdbr)
m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
m_t2g.reactome = m_df.reactome %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.kegg = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
m_t2g.kegg = m_df.kegg %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.hall = msigdbr(species = "Homo sapiens", category = "H")
m_t2g.hall = m_df.hall %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_df.regulate = msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
m_t2g.regulate = m_df.regulate %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
# convert annotation from data.frame to list format
kegg.list=list()
kegg.name=unique(m_t2g.kegg$gs_name)
for(i in kegg.name)
{
  kegg.list[[i]]=m_t2g.kegg$gene_symbol[m_t2g.kegg$gs_name==i]
}
kegg=read.gmt("KEGG_human_20190613_fromZhangbo.gmt")

packageVersion("qusage")
# [1] '2.24.0'
require(qusage)
  #seurat.all = a seurat object from scran
  #nm = name of the run
  #gs = gene set to test, a list
# the "run_qusage_heatmap" function needs that @meta.data$for_temp_color is a numeric format vector with MIN=1 rather than MIN=0, whic represent the cluster.ID
lum13@meta.data$for_temp_color=as.numeric(lum13@meta.data$fig.cluster)
run_qusage_heatmap.seurat3(lum13, nm = 'lum13.celltype', kegg, my.seed=100)


### name clusters
current.cluster.ids <- c(0:3)
new.cluster.ids <- c("tL0", "tL1", "tL2", "tLI3")
# tL: tumor luminal, tLI: tumor intemediate in luminal clusters
lum13@meta.data$celltype <- plyr::mapvalues(x = as.character(lum13@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)

# save obj
lum13@assays$RNA@scale.data=matrix(0, 0, 0)
lum13@commands=list()
lum13@meta.data$basal.sig=NULL
lum13@meta.data$cycle.sig=NULL
lum13@meta.data$SCT_snn_res.0.5=NULL
lum13@meta.data$SCT_snn_res.0.4=NULL
lum13@meta.data$SCT_snn_res.0.3=NULL
lum13@meta.data$SCT_snn_res.0.2=NULL
lum13@meta.data$SCT_snn_res.0.1=NULL
table(lum13@meta.data$celltype)
#   tL0   tL1   tL2  tLI3
# 11931  6763  4929   446
# saveRDS(lum13, "lum13.seurat3.indivSCTed.celltyped.rds")
# lum13=readRDS("lum13.seurat3.indivSCTed.celltyped.rds")


### plot heatmap for special categories of DEGs
subtype_markers=subtypeMarker(lum13)

markers=lum13.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) & (markers$pct.1>0.2 | markers$pct.2>0.5), ]
# for TFs, change log2(1.5) to log(1.2)
# lum13 <- ScaleData(lum13, features = rownames(lum13), assay="RNA" )
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
  if(nrow(top10)>0)
  {
    DefaultAssay(lum13)="RNA"
    pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
    temp.lum13=subset(lum13, downsample=300)
    temp.heatmap<-DoHeatmap(temp.lum13 , features = top10$gene, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0, size=3) 
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()
    DefaultAssay(lum13)="SCT"
  }
}

# for too few cluster DExpressing genes, do vlnplot
tocheck=c("CCL2", "CXCL11", "CXCL1", "CXCL8", "CXCL17", "AR")
pdf("VlnPlot.chemokines_celltype.pdf", width = 10, height =10)
temp=VlnPlot(lum13, features = tocheck, pt.size = 0, group.by="celltype", ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# check cell-surface markers of each cluster
markers=lum13.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ (markers$avg_log2FC>log2(1.5)) & (markers$pct.1>0.50) , ]
require(msigdbr)
cc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
cc = cc %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
cdgenes=sort(unique(cc$gene_symbol[cc$gs_name=="GO_CELL_SURFACE"]))
tempcheck=cdgenes
temp.markers=markers[markers$gene%in%tempcheck, ]
top10 <- temp.markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
if(nrow(top10)>0)
{
  DefaultAssay(lum13)="RNA"
  pdf(paste("temp.heatmap.CD_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
  temp.lum13=subset(lum13, downsample=300)
  temp.heatmap<-DoHeatmap(temp.lum13 , features = top10$gene, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0, size=3) 
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
  DefaultAssay(lum13)="SCT"
}


### combine all finer clusters in original clusters with original clusters
basal.anno=bas13@meta.data$celltype
names(basal.anno)=rownames(bas13@meta.data)
table(basal.anno)
   # tB1  tB2  tB4 tBC3 tBI0 tBI5
   # 213  165  109  115  219   87
luminal.anno=lum13@meta.data$celltype
names(luminal.anno)=rownames(lum13@meta.data)
table(luminal.anno)
  #   tL0   tL1   tL2  tLI3
  # 11931  6763  4929   446
table(prim10xe@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
# to detect pairwise DEGs
prim10xe@meta.data$celltype2compare=prim10xe@meta.data$celltype
prim10xe@meta.data[names(basal.anno), ]$celltype2compare=basal.anno
prim10xe@meta.data[names(luminal.anno), ]$celltype2compare=luminal.anno
table(prim10xe@meta.data$celltype2compare)
  # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
  # 213   165   109   115   219    87   632 11931  6763  4929   446


### get gleason annotation
table(prim10xe@meta.data$fig.patient)
 #  53   54   55   56   59   62   71   73   74   75   76   77
 # 990  267 1335 5820 1677 4246  994  340 2617 1000  932 5391
current.cluster.ids <- sort(unique(prim10xe@meta.data$fig.patient))
new.cluster.ids <- c("9", "10", "9", "9", "7", "7", "9", "9", "3+4", "9", "7", "9")
prim10xe@meta.data$gleason <- plyr::mapvalues(x = as.character(prim10xe@meta.data$fig.patient), from = current.cluster.ids, to = new.cluster.ids)
table(prim10xe@meta.data$gleason)
  #  10   3+4     7     9
  # 267  2617  6855 15870
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = FALSE, group.by="gleason", cols=c("green", "blue", "purple", "red"), order=c( "10", "9", "7", "3+4"), pt.size=0.1)
dev.off()

# saveRDS(prim10xe, "prim10xe.seurat3.indivSCTed.celltyped.rds")
# prim10xe=readRDS("prim10xe.seurat3.indivSCTed.celltyped.rds")


### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(prim10xe@meta.data$celltype2compare)
prim10xe = SetIdent(prim10xe, value = "celltype2compare")
for(i in celltype)
{
  for(j in celltype)
  {
    if(i!=j)
    {
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(prim10xe, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "prim10xe.markerlist.pairwise.rds")
# markerlist=readRDS("prim10xe.markerlist.pairwise.rds")


### pairwise fig.cluster's DEGs
markerlist=list()
celltype=unique(prim10xe@meta.data$fig.cluster)
prim10xe = SetIdent(prim10xe, value = "fig.cluster")
for(i in celltype)
{
  for(j in celltype)
  {
    if(i!=j)
    {
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(prim10xe, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "prim10xe.markerlist.pairwise.fig.cluster.rds")
# markerlist=readRDS("prim10xe.markerlist.pairwise.fig.cluster.rds")


### go to "_parts1and2_201216" to combine normal and PCa samples


### check basal/intermediate and luminal/intermediate signature in different-grade bulk samples
# identify signatures
table(prim10xe@meta.data$celltype2compare)
  # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
  # 213   165   109   115   219    87   632 11931  6763  4929   446
png("temp_cluster_dimplot.png", width=5*100, height=5*100)
DimPlot(prim10xe, label = TRUE, group.by="celltype2compare", pt.size=0.01)# + NoLegend()
dev.off()


### get gleason annotation
table(prim10xe@meta.data$fig.patient)
 #  53   54   55   56   59   62   71   73   74   75   76   77
 # 990  267 1335 5820 1677 4246  994  340 2617 1000  932 5391
current.cluster.ids <- sort(unique(prim10xe@meta.data$fig.patient))
new.cluster.ids <- c("9", "10", "9", "9", "7", "7", "9", "9", "3+4", "9", "7", "9")
prim10xe@meta.data$gleason <- plyr::mapvalues(x = as.character(prim10xe@meta.data$fig.patient), from = current.cluster.ids, to = new.cluster.ids)
table(prim10xe@meta.data$gleason)
  #  10   3+4     7     9
  # 267  2617  6855 15870
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = FALSE, group.by="gleason", cols=c("green", "blue", "purple", "red"), order=c( "10", "9", "7", "3+4"), pt.size=0.1)
dev.off()

###--->>>originating cell type based gene signature
# get originating-cell-type annotations
table(prim10xe@meta.data$celltype2compare)
  # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
  # 213   165   109   115   219    87   632 11931  6763  4929   446
current.cluster.ids <- sort(unique(prim10xe@meta.data$celltype2compare))
new.cluster.ids <- c("TypeC", "TypeC", "Basal", "TypeC", "Basal", "TypeC", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal")
prim10xe@meta.data$ori.type <- plyr::mapvalues(x = as.character(prim10xe@meta.data$celltype2compare), from = current.cluster.ids, to = new.cluster.ids)
table(prim10xe@meta.data$ori.type)
  # Basal Luminal   TypeC
  #   328   24701     580
png("dimplot.type.originating.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = FALSE, group.by="ori.type", pt.size=0.1)
dev.off()

### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(prim10xe@meta.data$ori.type)
  # Basal  Luminal   TypeC
  #   328    24701     580
pct.matrix=table(prim10xe@meta.data$fig.patient, prim10xe@meta.data$ori.type)
temp=apply(pct.matrix, 1, sum)
pct.matrix=round(apply(pct.matrix, 2, function(x) x/temp), 2)
temp.name=rownames(pct.matrix)
# temp.name
# [1] "53" "54" "55" "56" "59" "62" "71" "73" "74" "75" "76" "77"
GS=c(9, 10, 9, 9, 7.2, 7.2, 9, 9, 7.1, 9, 7.2, 9)
pct.matrix=pct.matrix[order(GS, decreasing=TRUE), ]
write.csv(pct.matrix, file="cell.ori.seurat3.merge__PCa13.samples.diffGS_part3________201217__plotCellPct.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:8, "TypeC"]), as.numeric(pct.matrix[9:12, "TypeC"]))
str(res)
# go to Excel to plot barplot and boxplot

# identify speficific pairwise DEGs for each originating-cell-type
table(prim10xe@meta.data$ori.type)
  # Basal Luminal   TypeC
  #   328   24701     580
ori.type=prim10xe@meta.data$ori.type
names(ori.type)=rownames(prim10xe@meta.data)
table(prim10x13@meta.data$ori.type)
  #      Basal Endothelial cell       Fibroblast          Luminal
  #        328             3212              408            24701
  # Macrophage        Mast cell    Myofibroblast           T cell
  #       2142             1443             1798             3135
  #      TypeC
  #        580
prim10x13@meta.data$ori.type=prim10x13@meta.data$celltype
prim10x13@meta.data[names(ori.type), ]$ori.type=ori.type
table(prim10x13@meta.data$ori.type)
cellnum=table(prim10x13@meta.data$ori.type)

markerlist=list()
celltype=unique(prim10x13@meta.data$ori.type)
prim10x13 = SetIdent(prim10x13, value = "ori.type")
for(i in celltype)
{
  for(j in celltype)
  {
    if(i!=j)
    {
      print("i:")
      print(i)
      print("j:")
      print(j)
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(prim10x13, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "prim10x13.markerlist.pairwise.ori.type.rds")
# markerlist=readRDS("prim10x13.markerlist.pairwise.ori.type.rds")

### get originating-cell-type signature
# markerlist=readRDS("prim10x13.markerlist.pairwise.ori.type.rds")
celltype=unique(prim10x13@meta.data$ori.type)
# 
signature=list()
for(j in celltype)
{
  geneset=c(rownames(prim10x13@assays$RNA))
  for(i in names(markerlist))
  {
    if(grepl(paste(j, "__", sep=""), i))
    {
      # print("j:")
      # print(j)
      # print("i:")
      # print(i)
      temp=markerlist[[i]]
      geneset=intersect(geneset, rownames(temp)[ (temp$p_val_adj<0.05) & (temp$avg_log2FC>0) & (temp$pct.2<0.05) ])# & (temp$pct.1>0.2) 
    }
  }
  print(str(geneset))
  signature[[j]]=geneset
}

# prim10x13 <- ScaleData(prim10x13, features = rownames(prim10x13), assay="RNA" )
# for( i in names(signature) )
prim10x13 = SetIdent(prim10x13, value = "ori.type")
for( i in c("Luminal", "TypeC") )
{
  temp = colMeans( prim10x13@assays$RNA@scale.data[signature[[i]] ,] )
  prim10x13@meta.data$for_temp_color=temp
  check.genes=c("for_temp_color")
  DefaultAssay(prim10x13)="RNA"
  png(paste("temp.featureplot", i, "png", sep="."), height=10*100, width=14*100)
  temp=FeaturePlot(prim10x13, features = check.genes)
  print(temp)
  dev.off()
  DefaultAssay(prim10x13)="SCT"
  png(paste("temp.VlnPlot", i, "png", sep="."), width=10*100, height=8*100)
  temp=VlnPlot(prim10x13, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2, assay = "RNA")
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95, assay="RNA")
  }
  print(temp)
  dev.off()
}

# write signatures to .txt files
write.table(signature[["TypeC"]], file="marker.TypeC.originating.txt", col.name=F, row.name=F, quote=F, sep='\t')
write.table(signature[["Luminal"]], file="marker.luminal.originating.txt", col.name=F, row.name=F, quote=F, sep='\t')
###<<<---originating cell type based gene signature

###--->>>easy cell type based gene signature
# get originating-cell-type annotations
table(prim10xe@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
png("dimplot.celltype.png", height=5*100, width=5*100)
DimPlot(prim10xe, label = FALSE, group.by="celltype", pt.size=0.1)
dev.off()

### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(prim10xe@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
pct.matrix=table(prim10xe@meta.data$fig.patient, prim10xe@meta.data$celltype)
temp=apply(pct.matrix, 1, sum)
pct.matrix=round(apply(pct.matrix, 2, function(x) x/temp), 2)
temp.name=rownames(pct.matrix)
# temp.name
# [1] "53" "54" "55" "56" "59" "62" "71" "73" "74" "75" "76" "77"
GS=c(9, 10, 9, 9, 7.2, 7.2, 9, 9, 7.1, 9, 7.2, 9)
pct.matrix=pct.matrix[order(GS, decreasing=TRUE), ]
write.csv(pct.matrix, file="cell.ori.seurat3.merge__PCa13.samples.diffGS_part3________201217__plotCellPct.csv", row.name=T)
res1 <- t.test(as.numeric(pct.matrix[1:8, "tBI"]), as.numeric(pct.matrix[9:12, "tBI"]))
# p.value=0.036
res2 <- t.test(as.numeric(pct.matrix[1:8, "tInt"]), as.numeric(pct.matrix[9:12, "tInt"]))
# p.value=0.571
str(res1)
str(res2)
# go to Excel to plot barplot and boxplot

# identify speficific pairwise DEGs for each originating-cell-type
table(prim10xe@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
celltype=prim10xe@meta.data$celltype
names(celltype)=rownames(prim10xe@meta.data)
table(prim10x13@meta.data$celltype)
    #   Basal cell Endothelial cell  Epithelial cell       Fibroblast
    #          964             3212              651              408
    # Luminal cell       Macrophage        Mast cell    Myofibroblast
    #        23994             2142             1443             1798
    #       T cell
    #         3135
prim10x13@meta.data$finer.type=prim10x13@meta.data$celltype
prim10x13@meta.data[names(celltype), ]$finer.type=celltype
table(prim10x13@meta.data$finer.type)


markerlist=list()
celltype=unique(prim10x13@meta.data$finer.type)
prim10x13 = SetIdent(prim10x13, value = "finer.type")
for(i in celltype)
{
  for(j in celltype)
  {
    if(i!=j)
    {
      print("i:")
      print(i)
      print("j:")
      print(j)
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(prim10x13, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "prim10x13.markerlist.pairwise.finer.type.rds")
# markerlist=readRDS("prim10x13.markerlist.pairwise.finer.type.rds")

### get originating-cell-type signature
# markerlist=readRDS("prim10x13.markerlist.pairwise.finer.type.rds")
celltype=unique(prim10x13@meta.data$finer.type)
# 
signature=list()
for(j in celltype)
{
  geneset=c(rownames(prim10x13@assays$RNA))
  for(i in names(markerlist))
  {
    if(grepl(paste(j, "__", sep=""), i))
    {
      # print("j:")
      # print(j)
      # print("i:")
      # print(i)
      temp=markerlist[[i]]
      geneset=intersect(geneset, rownames(temp)[ (temp$p_val_adj<0.05) & (temp$avg_log2FC>0) & (temp$pct.2<0.05) ])# & (temp$pct.1>0.2) 
    }
  }
  print(str(geneset))
  signature[[j]]=geneset
}

# prim10x13 <- ScaleData(prim10x13, features = rownames(prim10x13), assay="RNA" )
# for( i in names(signature) )
prim10x13 = SetIdent(prim10x13, value = "ori.type")
for( i in c("tBI", "tInt", "tLum") )
{
  temp = colMeans( prim10x13@assays$RNA@scale.data[signature[[i]] ,] )
  prim10x13@meta.data$for_temp_color=temp
  check.genes=c("for_temp_color")
  DefaultAssay(prim10x13)="RNA"
  png(paste("temp.featureplot", i, "png", sep="."), height=10*100, width=14*100)
  temp=FeaturePlot(prim10x13, features = check.genes)
  print(temp)
  dev.off()
  DefaultAssay(prim10x13)="SCT"
  png(paste("temp.VlnPlot", i, "png", sep="."), width=10*100, height=8*100)
  temp=VlnPlot(prim10x13, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2, assay = "RNA")
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95, assay="RNA")
  }
  print(temp)
  dev.off()
}

("tBI", "tInt", "tLum") )
# write signatures to .txt files
write.table(signature[["tInt"]], file="marker.tInt.txt", col.name=F, row.name=F, quote=F, sep='\t')
write.table(signature[["tLum"]], file="marker.tLum.txt", col.name=F, row.name=F, quote=F, sep='\t')
write.table(signature[["tBI"]], file="marker.tBI.txt", col.name=F, row.name=F, quote=F, sep='\t')
###<<<---easy cell type based gene signature


### check bulk-high-gleason signature in scRNA data
gleason.genes.up=readRDS("gleason.genes.up.rds")
gleason.genes.up=gleason.genes.up[gleason.genes.up%in%rownames(prim10xe@assays$SCT@scale.data)]
# gleason.genes.down=readRDS("gleason.genes.down.rds")
# gleason.genes.down=gleason.genes.down[gleason.genes.down%in%rownames(prim10xe@assays$SCT@scale.data)]
# prim10xe=ScaleData(prim10xe)
prim10xe@meta.data$bulk.gleason=colMeans(prim10xe@assays$SCT@scale.data[gleason.genes.up, ])#-colMeans(prim10xe@assays$SCT@scale.data[gleason.genes.down, ])
DefaultAssay(prim10xe)="RNA"
png("prim10xe.bulk.gleason.featureplot.png", width=4*100, height=5*100)
FeaturePlot(prim10xe, features = "bulk.gleason", pt.size=0.01)
dev.off()
DefaultAssay(prim10xe)="SCT"
# as celltype2compare
pdf("temp.VlnPlot.fig.cluster.pdf", width = 11, height =4)
tocheck=c("bulk.gleason")
temp=VlnPlot(prim10xe, features = tocheck, ncol = 2,  pt.size = 0, group.by="fig.cluster", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()
# as fig.cluster
pdf("temp.VlnPlot.celltype2compare.pdf", width = 11, height =4)
tocheck=c("bulk.gleason")
temp=VlnPlot(prim10xe, features = tocheck, ncol = 2,  pt.size = 0, group.by="celltype2compare", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()
# as celltype
pdf("temp.VlnPlot.celltype.pdf", width = 11, height =4)
tocheck=c("bulk.gleason")
temp=VlnPlot(prim10xe, features = tocheck, ncol = 2,  pt.size = 0, group.by="celltype", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


### survival and bulk-different-Gleason comparision

prim10xe@assays$RNA@scale.data=matrix(0, 0, 0)
prim10xe@commands=list()
prim10xe@meta.data$basal.sig=NULL
prim10xe@meta.data$cycle.sig=NULL
prim10xe@meta.data$SCT_snn_res.0.5=NULL
prim10xe@meta.data$SCT_snn_res.0.3=NULL
saveRDS(prim10xe, "prim10xe.seurat3.indivSCTed.celltyped.rds")
# prim10xe=readRDS("prim10xe.seurat3.indivSCTed.celltyped.rds")


### check indolent signature from Hansen

### NOT RUN yet <<<---

