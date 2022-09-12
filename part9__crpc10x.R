### Seurat V3
source("../functions.seurat.R")
require(ggplot2)
require(RColorBrewer)
require(BoutrosLab.plotting.general)
require(dplyr)
require(canprot)
require(Seurat)
require(monocle)
require(canprot)
require(GEOquery)
packageVersion("Seurat")
# [1] '4.0.1'

# working.directory:
setwd("/ ... ")

### read crpc10xary PCa samples
### process from raw counts
data.directory="/ ... /counts/crpc10x.10x"
dir.raw=list.files(path = data.directory, full.names = TRUE, recursive = F)
pos.1=regexpr(data.directory, dir.raw, perl=TRUE)
pos.1=attr(pos.1, "match.length")+2
pos.2=nchar(dir.raw)
sample.names=substring(dir.raw, pos.1, pos.2)
m <- regexpr("P[0-9]", sample.names, perl=TRUE)
sample.names=regmatches(sample.names, m)
inte.list=list()# the list prepared for further dataset integration
for(i in 1:length(dir.raw))
{
  temp=read.table(dir.raw[i], header = TRUE, row.names=1, stringsAsFactors=FALSE)
  temp=aggregate(temp, by = list(temp$Symbol), FUN = max, na.rm = TRUE)
  rownames(temp)=temp$Symbol
  temp=temp[, 3:ncol(temp)]
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
  temp<- subset(inte.list[[i]], subset = nFeature_RNA > 1000 & percent.mt < 30 ) # & top6k==TRUE)
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
  # inte.list[[i]]@meta.data$top6k=rep(TRUE, ncol(inte.list[[i]]))
  # temp=sort(inte.list[[i]]@meta.data$nFeature_RNA, decreasing=TRUE)[6000]
  # inte.list[[i]]@meta.data$top6k[inte.list[[i]]@meta.data$nFeature_RNA<=temp]=FALSE
  # print(table(inte.list[[i]]@meta.data$top6k))
  # inte.list[[i]]<- subset(inte.list[[i]], subset = nFeature_RNA > 500 & percent.mt < 60 ) # & top6k==TRUE) # --->>> this is last version
  inte.list[[i]]<- subset(inte.list[[i]], subset = nFeature_RNA > 1000 & percent.mt < 30 ) # & top6k==TRUE)
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
crpc10x=merge(x = inte.list[[1]], y = inte.list[2:length(inte.list)] )
crpc10x@meta.data$fig.sample=as.character(crpc10x@meta.data$orig.ident)
crpc10x@meta.data$fig.patient=as.character(crpc10x@meta.data$orig.ident)
crpc10x@meta.data$pheno=rep("crpc10x", nrow(crpc10x@meta.data))
# crpc10x@meta.data$pheno[crpc10x@meta.data$fig.patient%in%c("QJZ")]="NEPC"
# pat2gleason=data.frame(patient=sort(unique(crpc10x@meta.data$fig.patient)))
# pat2gleason$patient
# pat2gleason$gleason=factor(c(NA, NA, "5", "5", "5", "5", "2", "4", "5", "5", "2", NA, "3"))
# pat2gleason
# crpc10x@meta.data$gleason=pat2gleason$gleason[match(crpc10x@meta.data$fig.patient, pat2gleason$patient)]
# crpc10x@meta.data$gleason=as.character(crpc10x@meta.data$gleason)
# table(crpc10x@meta.data$gleason)
str(crpc10x@meta.data)
# 'data.frame':   22301 obs. of  13 variables:
#  $ orig.ident     : chr  "P1" "P1" "P1" "P1" ...
#  $ nCount_RNA     : num  8656 11793 32037 13953 16081 ...
#  $ nFeature_RNA   : int  2334 2757 4921 2694 3369 3008 4419 3154 3799 4757 ...
#  $ percent.mt     : num  3.06 2.48 5.04 6.31 2.67 ...
#  $ fig.sample     : chr  "P1" "P1" "P1" "P1" ...
#  $ fig.patient    : chr  "P1" "P1" "P1" "P1" ...
#  $ pheno          : chr  "crpc10x" "crpc10x" "crpc10x" "crpc10x" ...
#  $ nCount_SCT     : num  10970 11630 11265 11930 12145 ...
#  $ nFeature_SCT   : int  2334 2756 2844 2692 3367 3008 3783 3035 3284 2274 ...
#  $ SCT_snn_res.0.2: Factor w/ 14 levels "0","1","2","3",..: 4 2 5 1 2 2 2 1 ...
#  $ seurat_clusters: Factor w/ 14 levels "0","1","2","3",..: 4 2 5 1 2 2 2 1 ...
#  $ fig.cluster    : Factor w/ 14 levels "0","1","2","3",..: 4 2 5 1 2 2 2 1 ...
crpc10x=SCTransform(crpc10x, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting

# check most variable genes
str(crpc10x@assays$SCT@var.features)
# chr [1:3000] "TPSAB1" "TPSB2" "C6orf15" "TAGLN" "MGP" "APOD" "CHGA" ...
# top10 <- head(VariableFeatures(crpc10x, assay="SCT"), 10)
# # plot variable features with and without labels
# toplot<- VariableFeaturePlot(crpc10x, assay = "SCT", selection.method="sctransform")
# toplot <- LabelPoints(plot = toplot, points = top10, repel = TRUE)
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()
# str(crpc10x)

DefaultAssay(crpc10x)="SCT"
pdf("temp.PCA.pdf")
crpc10x <- RunPCA(crpc10x, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(crpc10x,  ndims = 100)
dev.off()

### PCA-based clustering
{
  # crpc10x <- RunUMAP(crpc10x, reduction = "pca", dims = 1:25, verbose = FALSE)
  # crpc10x <- FindNeighbors(crpc10x, reduction = "pca", dims = 1:25, verbose = FALSE)
  # crpc10x <- FindClusters(crpc10x, verbose = FALSE, resolution = 0.2)
  # table(crpc10x@active.ident)

  # crpc10x.markers <- FindAllMarkers(crpc10x, assay="SCT")
  # saveRDS(crpc10x.markers, "crpc10x.markers.as.fig.cluster.byPCA.rds")
  # crpc10x.markers=readRDS("crpc10x.markers.as.fig.cluster.byPCA.rds")
}

### Harmony-based clustering
{
  require(harmony)
  DefaultAssay(crpc10x)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  crpc10x <- RunHarmony(crpc10x, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  dev.off()
  Sys.time()

  # crpc10x <- RunUMAP(crpc10x, dims = 1:25, verbose = FALSE)
  crpc10x <- RunUMAP(crpc10x, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpc10x <- FindNeighbors(crpc10x, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpc10x <- FindClusters(crpc10x, verbose = FALSE, resolution = 0.2)
  table(crpc10x@active.ident)
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 11897 10508  3867  2841  1996  1815  1803  1503  1345   933   861   555   511
#    13    14    15    16
#   322   316   206   198
  crpc10x.markers <- FindAllMarkers(crpc10x, assay="SCT")
  saveRDS(crpc10x.markers, "crpc10x.markers.as.fig.cluster.byHarmony.rds")
  # crpc10x.markers =readRDS("crpc10x.markers.as.fig.cluster.byHarmony.rds")
}


crpc10x@meta.data$fig.cluster=crpc10x@meta.data$seurat_clusters
crpc10x <- SetIdent(crpc10x, value = "fig.cluster")
png("temp_cluster_dimplot.png", width=5*100, height=5*100)
DimPlot(crpc10x, label = TRUE, group.by="fig.cluster", pt.size=0.001)# + NoLegend()
dev.off()
png("temp_patient_dimplot.png", width=6*100, height=5*100)
DimPlot(crpc10x, label = FALSE, group.by="fig.sample", pt.size=0.001)# + NoLegend()
dev.off()
png("temp_pheno_dimplot.png", width=6*100, height=5*100)
DimPlot(crpc10x, label = FALSE, group.by="pheno", pt.size=0.001)# + NoLegend()
dev.off()


### check celltype markers' expression to make sure every cluster's celltype
markers=crpc10x.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[markers$avg_log2FC>log2(1.5), ]
setwd("..")
source("functions.seurat.R")
crpc10x@meta.data$celltype=celltypeAnno(crpc10x, markers)
setwd("./ ... ")
  # [3] "MDSCfunction_profRen: 0.57" EPCAM+KRT8+ -->epithelial
  # [5] "MDSCfunction_profRen: 0.3" TOP2A+STMN1+SOX2+ -->epithelial(CellCycle)
  # [13] "MDSCfunction_profRen: 1.68" EPCAM+KRT8+KRT19+ -->epithelial
  # [14] "Basal cell: 0.39 Fibroblast_cell201801: 0.38" VIM+SOX2+ -->stem like
# 
crpc10x@meta.data$celltype[crpc10x@meta.data$fig.cluster%in%c("2", "4", "12")]="Epithelial cell"
# plot with cluster names
png("temp_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(crpc10x, label = TRUE, group.by="celltype") + NoLegend()
dev.off()


# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=crpc10x.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.crpc10x.cluster","csv",sep="."), row.name=T)


# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1", "UPK2", "AR")
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # PCa markers
# check.genes=c("FOLH1", "SCHLAP1", "AMACR","CCL2", "CD74", "LY6D", "KLK3", "AR", "ACPP")
# check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

# crpc10x=SetIdent(crpc10x, value = "fig.cluster")
crpc10x=SetIdent(crpc10x, value = "celltype")
DefaultAssay(crpc10x)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpc10x, features = check.genes)
dev.off()
DefaultAssay(crpc10x)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(crpc10x, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(crpc10x, "crpc10x.rds")
# crpc10x=readRDS("crpc10x.rds")


### further analysis focused on epithelial cells
table(crpc10x@meta.data$celltype)
   # B/Plasma cell       Basal cell Endothelial cell  Epithelial cell
   #           429             2273             1253             3721
   #    Fibroblast     Luminal cell       Macrophage        Mast cell
   #          1825             8258             1267              905
   # Myofibroblast           T cell
   #          1289             1081
epitypes=c("Epithelial cell", "Luminal cell", "Basal cell")
crpc10xe=subset(crpc10x, subset=(celltype%in%epitypes))
# str(crpc10xe)
crpc10xe@commands=list()
crpc10xe@meta.data$pheno=rep("CRPC", nrow(crpc10xe@meta.data))
crpc10xe@meta.data$pheno[crpc10xe@meta.data$fig.patient%in%c("P2", "P5", "P6")]="NEPC"

dim(crpc10xe)
# [1] 29262 14252
DefaultAssay(crpc10xe)="SCT"

crpc10xe=SCTransform(crpc10xe, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
pdf("temp.PCA.pdf")
crpc10xe <- RunPCA(crpc10xe, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(crpc10xe,  ndims = 100)
dev.off()

### clustering basd on PCA
{
  # DefaultAssay(crpc10xe)="SCT"
  # crpc10xe <- RunUMAP(crpc10xe, dims = 1:20, verbose = FALSE)
  # crpc10xe <- FindNeighbors(crpc10xe, dims = 1:20, verbose = FALSE)

  # DefaultAssay(crpc10xe)="SCT" 
  # crpc10xe <- FindClusters(crpc10xe, verbose = FALSE, resolution = 0.2)
  # table(crpc10xe@meta.data$seurat_clusters)
  # crpc10xe.markers <- FindAllMarkers(crpc10xe, assay="SCT")
  # saveRDS(crpc10xe.markers, "crpc10xe.markers.as.fig.cluster.byPCA.rds")
  # crpc10xe.markers =readRDS("crpc10xe.markers.as.fig.cluster.byPCA.rds")

}

### clustering basd on Harmony
{
  require(harmony)
  DefaultAssay(crpc10xe)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  crpc10xe <- RunHarmony(crpc10xe, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
  dev.off()
  Sys.time()

  # crpc10x <- RunUMAP(crpc10x, dims = 1:25, verbose = FALSE)
  crpc10xe <- RunUMAP(crpc10xe, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpc10xe <- FindNeighbors(crpc10xe, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpc10xe <- FindClusters(crpc10xe, verbose = FALSE, resolution = 0.2)
  table(crpc10xe@meta.data$seurat_clusters)
  #    0    1    2    3    4    5    6    7    8
  # 3813 2438 2071 1920 1519 1123  568  558  242
  crpc10xe.markers <- FindAllMarkers(crpc10xe, assay="SCT")
  saveRDS(crpc10xe.markers, "crpc10xe.markers.as.fig.cluster.byHarmony.rds")
  # crpc10xe.markers =readRDS("crpc10xe.markers.as.fig.cluster.byHarmony.rds")
}

markers=crpc10xe.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.crpc10xe.cluster","csv",sep="."), row.name=T)

crpc10xe@meta.data$fig.cluster=crpc10xe@meta.data$seurat_clusters
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(crpc10xe, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
png("dimplot.patient.png", height=5*100, width=5*100)
DimPlot(crpc10xe, label = FALSE, group.by="fig.patient")
dev.off()
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(crpc10xe, label = FALSE, group.by="gleason")# + NoLegend()
dev.off()
png("dimplot.pheno.png", height=5*100, width=5*100)
DimPlot(crpc10xe, label = FALSE, group.by="pheno")# + NoLegend()
dev.off()


### analysis of crpc10xtasis STOP here, because the small sample number (206 cells)
# saveRDS(crpc10xe, "crpc10xe.rds")
# crpc10xe=readRDS("crpc10xe.rds")


# check known epithelial cell types
check.genes=c("KLK3", "KRT8", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1", "UPK2", "AR")
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # check genes 
# check.genes=c("CCL2", "CD74")
# check.genes=c("FGFR3", "TGFBR3", "UPK1A", "CD55", "SRPX2", "WNT5A")
# check.genes=c("TGFBR3", "UPK1A", "SRPX2")
# check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3")
# check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")

crpc10xe <- SetIdent(crpc10xe, value = "fig.cluster")
DefaultAssay(crpc10xe)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpc10xe, features = check.genes)
dev.off()
DefaultAssay(crpc10xe)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(crpc10xe, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# indi.plot(crpc10xe, tocheck="fig.patient")
# indi.plot(crpc10xe, tocheck="fig.cluster")

### re-cluster basal & typeC cells
{
  crpc10xb=subset(crpc10xe, subset=(fig.cluster%in%c("2")))
  dim(crpc10xb)
  # [1] 27945  2071
  DefaultAssay(crpc10xb)="SCT"

  crpc10xb=SCTransform(crpc10xb, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  pdf("temp.PCA.pdf")
  crpc10xb <- RunPCA(crpc10xb, assay="SCT", verbose = FALSE, npcs = 100)
  ElbowPlot(crpc10xb,  ndims = 100)
  dev.off()

  ### clustering basd on PCA
  {
    # DefaultAssay(crpc10xb)="SCT"
    # crpc10xb <- RunUMAP(crpc10xb, dims = 1:20, verbose = FALSE)
    # crpc10xb <- FindNeighbors(crpc10xb, dims = 1:20, verbose = FALSE)

    # DefaultAssay(crpc10xb)="SCT" 
    # crpc10xb <- FindClusters(crpc10xb, verbose = FALSE, resolution = 0.2)
    # table(crpc10xb@meta.data$seurat_clusters)
    # crpc10xb.markers <- FindAllMarkers(crpc10xb, assay="SCT")
    # saveRDS(crpc10xb.markers, "crpc10xb.markers.as.fig.cluster.byPCA.rds")
    # crpc10xb.markers =readRDS("crpc10xb.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(crpc10xb)
    # [1] "SCT"
    Sys.time()
    pdf("temp.runHarmony.pdf")
    crpc10xb <- RunHarmony(crpc10xb, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    dev.off()
    Sys.time()

    # crpc10x <- RunUMAP(crpc10x, dims = 1:25, verbose = FALSE)
    crpc10xb <- RunUMAP(crpc10xb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    crpc10xb <- FindNeighbors(crpc10xb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    crpc10xb <- FindClusters(crpc10xb, verbose = FALSE, resolution = 0.2)
    table(crpc10xb@meta.data$seurat_clusters)
    #   0   1   2   3   4   5
    # 811 604 363 201  63  29
    # crpc10xb.markers <- FindAllMarkers(crpc10xb, assay="SCT")
    # saveRDS(crpc10xb.markers, "crpc10xb.markers.as.fig.cluster.byHarmony.rds")
    # crpc10xb.markers =readRDS("crpc10xb.markers.as.fig.cluster.byHarmony.rds")
  }

  crpc10xb@meta.data$fig.cluster=crpc10xb@meta.data$seurat_clusters
  png("dimplot.cluster.png", height=5*100, width=5*100)
  DimPlot(crpc10xb, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  dev.off()
  png("dimplot.patient.png", height=5*100, width=5*100)
  DimPlot(crpc10xb, label = TRUE, group.by="fig.patient")
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  DimPlot(crpc10xb, label = TRUE, group.by="gleason")
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  DimPlot(crpc10xb, label = TRUE, group.by="fig.pheno")
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
  # check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
  # check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3")
  # check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "UPK1A", "CD55", "TGFBR3", "FGFR3")

  crpc10xb <- SetIdent(crpc10xb, value = "fig.cluster")
  DefaultAssay(crpc10xb)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(crpc10xb, features = check.genes)
  dev.off()
  DefaultAssay(crpc10xb)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=8*100)
  temp=VlnPlot(crpc10xb, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  ### name clusters
  current.cluster.ids <- c(0:5)
  new.cluster.ids <- c("B20", "I21", "B22", "L23", "I24", "I25")
  # b: Basal, i: Intermediate, l: Luminal
  crpc10xb@meta.data$epitype <- plyr::mapvalues(x = as.character(crpc10xb@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  png("dimplot.epitype.png", height=5*100, width=5*100)
  DimPlot(crpc10xb, label = TRUE, group.by="epitype", label.size=8)
  dev.off()

  btype=crpc10xb@meta.data$epitype
  names(btype)=rownames(crpc10xb@meta.data)
  table(btype)
  # btype
  # B20 B22 I21 I24 I25 L23
  # 811 363 604  63  29 201

  # indi.plot(crpc10xb, tocheck="fig.patient")
  # indi.plot(crpc10xb, tocheck="fig.cluster")

  # saveRDS(crpc10xb, "crpc10xb.rds")
  # crpc10xb=readRDS("crpc10xb.rds")
}


### epitype annotation
current.cluster.ids <- c(0:8)
new.cluster.ids <- c("L0", "L1", "B2", "L3","NE4", "NE5", "I6", "NE7", "L8")
# i: Intermediate, l: Luminal, ne: Neuroendocrine
crpc10xe@meta.data$epitype <- plyr::mapvalues(x = as.character(crpc10xe@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot
table(names(btype)%in%rownames(crpc10xe@meta.data))
crpc10xe@meta.data[names(btype), "epitype"]=btype
pdf("dimplot.crpc10xe.epitype.pdf", height=5, width=5)
DimPlot(crpc10xe, label = FALSE, group.by="epitype", label.size=8, raster=TRUE)
dev.off()


### celltype annotation # add refined basal/typeC/interm types into epitype slot
current.cluster.ids <- unique(crpc10xe@meta.data$epitype)
current.cluster.ids=sort(current.cluster.ids)
current.cluster.ids
#  [1] "B20" "B22" "I21" "I24" "I25" "I6"  "L0"  "L1"  "L23" "L3"  "L8"  "NE4"
# [13] "NE5" "NE7"
new.cluster.ids <- rep("Luminal", length(current.cluster.ids))
new.cluster.ids[grepl("B", current.cluster.ids)]="Basal"
new.cluster.ids[grepl("C", current.cluster.ids)]="TypeC"
new.cluster.ids[grepl("I", current.cluster.ids)]="Intermediate"
new.cluster.ids[grepl("NE", current.cluster.ids)]="NE"
crpc10xe@meta.data$celltype <- plyr::mapvalues(x = as.character(crpc10xe@meta.data$epitype), from = current.cluster.ids, to = new.cluster.ids)
pdf("dimplot.crpc10xe.celltype.pdf", height=5, width=5)
DimPlot(crpc10xe, label = FALSE, group.by="celltype", label.size=5, raster=TRUE)
dev.off()

crpc10xe <- SetIdent(crpc10xe, value = "celltype")
DefaultAssay(crpc10xe)="SCT"

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
# check.genes=c("SIX3", "MIAT", "SG3", "MIAT", "CHD7", "KIF19", "NEB", "RGS16", "PCP4", "GPX2", "ID4", "CALB1", "SGB2A1")
# check.genes=c("ASCL1", "NEUROD1")
# check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
# check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3")


png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpc10xe, features = check.genes)
dev.off()
DefaultAssay(crpc10xe)="SCT"
png("temp.VlnPlot.png", width=10*100, height=12*100)
temp=VlnPlot(crpc10xe, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(crpc10xe, "crpc10xe.rds")
# crpc10xe=readRDS("crpc10xe.rds")


### barplot the percentage of each samples in cell types
{
  # give sample ID as phenotype + number
  table(crpc10xe@meta.data$pheno)
  # CRPC NEPC
  # 5464 8788
  table(crpc10xe@meta.data$fig.patient);unique(crpc10xe@meta.data$fig.patient)
  crpc10xe@meta.data$pheno=factor(crpc10xe@meta.data$pheno, levels=c("CRPC", "NEPC"))
  pheno.sample=crpc10xe@meta.data[, c("pheno", "fig.sample")]
  pheno.sample$fig.sample=as.character(pheno.sample$fig.sample)
  str(pheno.sample)
  # 'data.frame':   14252 obs. of  2 variables:
  #  $ pheno     : Factor w/ 2 levels "CRPC","NEPC": 1 1 1 1 1 1 1 1 1 1 ...
  #  $ fig.sample: chr  "P1" "P1" "P1" "P1" ...
  pheno.sample=pheno.sample[order(pheno.sample$pheno, pheno.sample$fig.sample, decreasing=c(FALSE, FALSE)), ]
  unique.sample=unique(pheno.sample$fig.sample)
  pheno.sample$fig.sample=factor(pheno.sample$fig.sample, levels=unique.sample)
  pheno.sample$fig.sample.id=as.numeric(pheno.sample$fig.sample)
  pheno.sample$pheno.sample=paste(as.character(pheno.sample$pheno), pheno.sample$fig.sample.id)
  unique.pheno.sample=unique(pheno.sample$pheno.sample)
  pheno.sample$pheno.sample=factor(pheno.sample$pheno.sample, levels=unique.pheno.sample)
  crpc10xe@meta.data$pheno.sample=pheno.sample[rownames(crpc10xe@meta.data), ]$pheno.sample

  celltype.sta=table(crpc10xe@meta.data$pheno.sample, crpc10xe@meta.data$celltype)
  celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
  # celltype.sta=celltype.sta[ c("Clike", "TypeC", "Basal", "NE"), !grepl("incidental|Normal", colnames(celltype.sta), perl=TRUE) ]
  celltype.sta=celltype.sta[ c("Basal",  "NE"),]
  barN=ncol(celltype.sta)
  n_top = nrow(celltype.sta) 
  pdf(paste("crpc10xe.barplot.celltypeInSample.pdf",sep=""))
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


temp.marker=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
temp.marker=c("UPK1A", "CD55", "TGFBR3", "FGFR3")
temp=subset(crpc10xe, downsample=1000)
# temp=NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA" )
# temp=ScaleData(temp, features = rownames(temp), assay="RNA" )
levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
pdf(paste("temp.heatmap.pdf", sep=""), width=(length(unique(crpc10xe@active.ident))/4+5), height=length(temp.marker)/10+3)
temp.heatmap<-DoHeatmap(temp , features = temp.marker, slot="scale.data", assay="SCT", angle = 0, size=3) 
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

# crpc10xe.markers=readRDS("crpc10xe.markers.celltype.rds")
crpc10xe.markers[crpc10xe.markers$cluster=="Clike" & crpc10xe.markers$gene=="TGFBR3", ]
#              p_val avg_log2FC pct.1 pct.2   p_val_adj cluster   gene
# TGFBR3 2.03857e-07  0.2556965   0.3 0.104 0.003816204   Clike TGFBR3
crpc10xe.markers[crpc10xe.markers$cluster=="Clike" & crpc10xe.markers$gene=="FGFR3", ]
             # p_val avg_log2FC pct.1 pct.2    p_val_adj cluster  gene
# FGFR3 5.369109e-49  0.6667128 0.457 0.053 1.005097e-44   Clike FGFR3
 
### check specific marker across all cell types
# crpc10x=readRDS("crpc10x.rds")
check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
clike.cells=rownames(crpc10xe@meta.data)[crpc10xe@meta.data$celltype=="Clike"]
crpc10x@meta.data$celltype[rownames(crpc10x@meta.data)%in%clike.cells]="Clike"
crpc10x <- SetIdent(crpc10x, value = "celltype")
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpc10x, features = check.genes)
dev.off()
DefaultAssay(crpc10x)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(crpc10x, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()



### Not employed --->>>
  ### CCA to identify epithelial cell types WithiN luminal tumor cells
  # norme=readRDS("norme.rds")
  crpc10xl=subset(crpc10xe, subset=celltype=="Luminal")
  anchors <- FindTransferAnchors(reference = norme, query = crpc10xl, 
      dims = 1:30, normalization.method="LogNormalize", reference.assay="SCT")
  predictions <- TransferData(anchorset = anchors, refdata = norme$celltype, 
      dims = 1:30)
  crpc10xl@meta.data$predtype=predictions$predicted.id 
  png("dimplot.predtype.png", height=5*100, width=5*100)
  DimPlot(crpc10xl, label = TRUE, group.by="predtype")
  dev.off()
### Not employed <<<---


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(crpc10xe@meta.data$celltype)
# CellCycle     Clike   Luminal        NE
#       795        70      5422       132
table(crpc10xe@meta.data$fig.patient)
# DGY-PG    HYL    QJZ    XYM
#   3489   1195    169   1566
crpc10xe@meta.data$fig.patient=factor(crpc10xe@meta.data$fig.patient, levels=c("QJZ", "DGY-PG", "HYL", "XYM"))
pct.matrix=table(crpc10xe@meta.data$fig.patient, crpc10xe@meta.data$celltype)
temp=apply(pct.matrix, 1, sum)
pct.matrix=round(apply(pct.matrix, 2, function(x) x/temp), 4)
temp.name=rownames(pct.matrix)
temp.name
#  [1] "CZK-T"      "DHB-T"      "HYQ"        "LPX"        "PXL"
#  [6] "QLX-M"      "R-M"        "RSC1023"    "S6-TUMOR-M" "SHSCR"
# [11] "SJ-M"       "T"          "ZQF-M"
### 9.5: N1.Gleason5
# GS=c( NA, NA, 9.5, 5, 9.5, 9.5, 2, 9.4, 5, 5, 2, NA, 3)
# pct.matrix=pct.matrix[order(GS, decreasing=TRUE), ]
write.csv(pct.matrix, file="part4__crpc10xe__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(crpc10xe@meta.data$epitype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(crpc10xe@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(crpc10xe@meta.data$fig.patient, crpc10xe@meta.data$epitype)
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
write.csv(pct.matrix, file="part4__crpc10xe__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


# saveRDS(crpc10xe, "crpc10xe.rds")
# crpc10xe=readRDS("crpc10xe.rds")


### celltype-specific DEG identification
table(crpc10xe@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   12893     925
DefaultAssay(crpc10xe)="SCT"
### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(crpc10xe@meta.data$celltype)
celltype=sort(celltype)
crpc10xe = SetIdent(crpc10xe, value = "celltype")
for(i in celltype)
{
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        markerlist[[temp.name]]=FindMarkers(crpc10xe, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
      }
  }
}
# saveRDS(markerlist, "crpc10xe.celltype.pairwise.rds")
# markerlist=readRDS("crpc10xe.celltype.pairwise.rds")

sp.markerlist=list()
celltype=unique(crpc10xe@meta.data$celltype)
celltype=sort(celltype)
for(i in celltype)
{
  temp.genes=rownames(crpc10xe@assays$SCT@data)
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

temp=subset(crpc10xe, downsample=300)
sp.markerlist=sp.markerlist[c("Basal", "Intermediate", "Luminal", "NE")]
tocheck=Reduce(union, sp.markerlist)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("Basal", "Intermediate", "Luminal", "NE")
pdf(paste("heatmap.crpc10xe.typeC.pdf", sep=""),  width=10, height=length(tocheck)/3)
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


### DEG identification, across all clusters
crpc10xe <- SetIdent(crpc10xe, value = "celltype")
crpc10xe.markers <- FindAllMarkers(crpc10xe, slot="data", assay="SCT")
# saveRDS(crpc10xe.markers, "crpc10xe.markers.celltype.rds")
# crpc10xe.markers=readRDS("crpc10xe.markers.celltype.rds")

### plot heatmap
# crpc10xe.markers=readRDS("crpc10xe.indivSCT.markers.as.fig.cluster.rds")
# crpc10xe <- NormalizeData(crpc10xe, normalization.method = "LogNormalize", scale.factor = 10000, assay="SCT" )
# crpc10xe <- ScaleData(crpc10xe, features = rownames(crpc10xe), assay="RNA")

# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=crpc10xe.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.crpc10xe.cluster","csv",sep="."), row.name=T)

enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)


# heatmap top 10 DEGs for every cluster
markers=crpc10xe.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(2) & markers$pct.1>0.3 & markers$pct.2<0.1, ]
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
top10$cluster=factor(top10$cluster, levels=c("Clike", "Luminal", "NE", "CellCycle"))
top10=top10[order(top10$cluster), ]
temp=subset(crpc10xe, downsample=300)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
tocheck=c(top10$gene, "KLK3", "AR")
tocheck=c(tocheck, "KRT4", "TACSTD2", "PSCA")
pdf(paste("heatmap.crpc10xe.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=length(tocheck)/8)
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### plot heatmap for special categories of DEGs
setwd("..")
subtype_markers=subtypeMarker(crpc10xe)
setwd("./ ... ")
crpc10xe.markers=readRDS("crpc10xe.markers.celltype.rds")
markers=crpc10xe.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
# for TFs, change log2(1.5) to log(1.2)
# crpc10xe <- ScaleData(crpc10xe, features = rownames(crpc10xe), assay="RNA" )
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
   top10$cluster=factor(top10$cluster, levels=c("Clike", "Luminal", "NE", "CellCycle"))
  top10=top10[order(top10$cluster), ]
  if(nrow(top10)>0)
  {
    DefaultAssay(crpc10xe)="SCT"
    pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
    temp=subset(crpc10xe, downsample=300)
    temp <- SetIdent(temp, value = "celltype")
    levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
    temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()
    DefaultAssay(crpc10xe)="SCT"
  }
}

# for too few cluster DExpressing genes, do vlnplot
tocheck=c("CCL2", "CXCL11", "CXCL1", "CXCL8", "CXCL17", "AR")
pdf("VlnPlot.chemokines_celltype.pdf", width = 10, height =10)
temp=VlnPlot(crpc10xe, features = tocheck, pt.size = 0, group.by="celltype", ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# enrichment of Clike compared with all other cells
markers=crpc10xe.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)

# enrichment of Clike compared with TypeC
markerlist=readRDS("crpc10xe.celltype.pairwise.rds")
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
crpc10xe@meta.data$for_temp_color=as.numeric(factor(crpc10xe@meta.data$celltype, levels=c("Clike", "Luminal", "NE",  "CellCycle")))
run_qusage_heatmap.seurat3(crpc10xe, nm = 'crpc10xe.celltype', kegg, my.seed=100)
run_qusage_heatmap.seurat3(crpc10xe, nm = 'crpc10xe.celltype', hall.list, my.seed=100)


# read CNV result
{
  ### inferCNV by "part5__inferCNV.aggAndIndo.epitype.sct.counts.R"
  ### read inferCNV result (based on $SCT@counts), and classify cancer and normal cells
  # crpc10xe=readRDS("crpc10xe.rds")
  cnv=list()
  cnv[[1]]=read.table("../inferCNV/crpc10x.sct.counts.out/infercnv.observations.txt")
  cnv[[2]]=read.table("../inferCNV/crpc10x.sct.counts.out/infercnv.references.txt")
  for(i in 1:2)
  {
    cnv[[i]]=round((cnv[[i]]-1)^2, 3)
    cnv[[i]]=apply(cnv[[i]], 2, mean)
  }
  str(cnv)
  cnv=c(cnv[[1]], cnv[[2]])

  table(names(cnv)%in%colnames(crpc10xe@assays$SCT@counts))
  #  TRUE
  # 18950
  # str(names(cnv)[!names(cnv)%in%colnames(crpc10xe@assays$SCT@counts)])
  # # chr [1:3354] "GCTGGGTTCGTCGTTC.1_1" "GGCCGATGTGCCTTGG.1_1" ...
  # str(colnames(crpc10xe@assays$SCT@counts)[grepl("^G", colnames(crpc10xe@assays$SCT@counts), perl=TRUE)])
  # # chr [1:781] "GAACCTACACCAGGTC-1_1" "GAACGGAGTCTAAAGA-1_1" ...

  # names.temp=gsub("\\.", "-", names(cnv))
  # str(names.temp[!names.temp%in%colnames(crpc10xe@assays$SCT@counts)])
  # # chr(0)
  # table(names.temp%in%colnames(crpc10xe@assays$SCT@counts))
  # #  TRUE
  # # 14439
  # table(colnames(crpc10xe@assays$SCT@counts)%in%names.temp)
  # #  TRUE
  # # 14439
  # names(cnv)=names.temp

  temp=crpc10xb
  cnv.sub=cnv[names(cnv)%in%rownames(temp@meta.data)]
  temp@meta.data$cnv.counts=rep(0, nrow(temp@meta.data))
  temp@meta.data[names(cnv.sub), "cnv.counts"]=cnv.sub

  # png("temp.VlnPlot.CNV.png", width=12*100, height=8*100)
  pdf("temp.VlnPlot.CNV_COUNTS.pdf", width=12, height=8)
  check.genes="cnv.counts"
  temp=VlnPlot(temp, features = check.genes, pt.size = 0, ncol = 1, group.by = "fig.cluster")
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()
}


### monocle2
{
  # for NE and crpc10x individually
  require(monocle)
  # crpc10xe=readRDS("crpc10xe.rds")
  temp=crpc10xe
  # temp=subset(crpc10xe, subset=(pheno=="crpc10x"))
  gene_metadata=data.frame(gene_short_name=rownames(temp@assays$SCT@counts))
  rownames(gene_metadata)=rownames(temp@assays$SCT@counts)
  temp.m <- newCellDataSet(  temp@assays$SCT@counts,
                            phenoData = new("AnnotatedDataFrame", temp@meta.data),
                            featureData = new("AnnotatedDataFrame", gene_metadata),
                            expressionFamily=negbinomial.size() )

  table(temp.m@phenoData@data$epitype)
  table(temp.m@phenoData@data$celltype)
  temp.m <- estimateSizeFactors(temp.m)
  temp.m <- estimateDispersions(temp.m)

  # temp.m <- detectGenes(temp.m, min_expr = 0.1)
  # ### only keep expressed genes
  # expressed_genes <- row.names(temp.m)[temp.m@featureData@data$num_cells_expressed>= 10]
  # temp.m <- temp.m[expressed_genes,]

  ### use all significant markers of clusters as ordering genes
  tumor.markers=readRDS("crpc10xe.markers.celltype.rds")
  markers=tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[abs(markers$avg_log2FC)>log2(1.5), ]
  markers$foldChange=2^(markers$avg_log2FC)
  # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
  markers=markers[order(markers$foldChange, decreasing=TRUE), ]
  ordering.genes=unique(markers$gene)
  # ordering.genes=unique(rownames(markers))
  str(ordering.genes)
  # chr [1:866] "LCN2" "BPIFB1" "CHGA" "SCGB1A1" "PSCA" "SCGN" "FCGBP" ...

  temp.m <- setOrderingFilter(temp.m,  ordering.genes) # 
  pdf("plot_ordering_genes.pdf")
  plot_ordering_genes(temp.m)
  dev.off()

  temp.m <- reduceDimension(temp.m, max_components = 2,
      method = 'DDRTree')

  temp.m <- orderCells(temp.m)
  pdf("temp_trajectory_epitype.pdf")
  plot_cell_trajectory(temp.m, color_by = "epitype")
  dev.off()
  pdf("temp_trajectory_celltype.pdf")
  plot_cell_trajectory(temp.m, color_by = "celltype")
  dev.off()

  # set a indicator for time, and sort again
  temp.m@phenoData@data$timeIset=rep(1, ncol(temp.m))
  temp.m@phenoData@data$timeIset[temp.m@phenoData@data$celltype=="Clike"]=0
  # temp.m@phenoData@data$timeIset[temp.m@phenoData@data$celltype=="Luminal"]=0
  temp.m<-orderCells(temp.m, root_state = temp.m@phenoData@data$timeIset)

  pdf("plot_cell_trajectory_byPseudotime_temp.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  plot_cell_trajectory(temp.m, color_by = "Pseudotime", cell_size=0.4)
  dev.off()

  pdf("plot_cell_trajectory_fig.cluster_temp.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  plot_cell_trajectory(temp.m, color_by = "epitype", cell_size=0.4)
  dev.off()

  pdf("plot_cell_trajectory_celltype_temp.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  plot_cell_trajectory(temp.m, color_by = "celltype", cell_size=0.05, cell_name_size = 8)
  dev.off()

  pdf("plot_cell_trajectory_fig.zone_temp.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  plot_cell_trajectory(temp.m, color_by = "fig.zone", cell_size=0.05, cell_name_size = 8)
  dev.off()

  pdf("plot_cell_trajectory_pheno_temp.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  plot_cell_trajectory(temp.m, color_by = "pheno", cell_size=0.1, cell_name_size = 8)
  dev.off()

  temp.m.time=temp.m@phenoData@data$Pseudotime
  names(temp.m.time)=rownames(temp.m@phenoData@data)

  png("plot_cell_trajectory_details_celltype_temp.m.png", width=30*100, height=6*100)
  temp=plot_cell_trajectory(temp.m, color_by = "celltype") +
      facet_wrap(~celltype, nrow = 1)
  print(temp)
  dev.off()

  png("plot_cell_trajectory_details_zone_temp.m.png", width=15*100, height=6*100)
  temp=plot_cell_trajectory(temp.m, color_by = "fig.zone") +
      facet_wrap(~fig.zone, nrow = 1)
  print(temp)
  dev.off()

  png("plot_cell_trajectory_details_epitype_temp.m.png", width=15*100, height=6*100)
  temp=plot_cell_trajectory(temp.m, color_by = "epitype") +
      facet_wrap(~epitype, nrow = 1)
  print(temp)
  dev.off()

  png("plot_cell_trajectory_details_cluster_temp.m.png", width=30*100, height=5*100)
  temp=plot_cell_trajectory(temp.m, color_by = "fig.cluster", cell_size = 0.8) +
      facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()
  pdf("plot_cell_trajectory_details_cluster_temp.m.pdf", width=30, height=5)
  temp=plot_cell_trajectory(temp.m, color_by = "fig.cluster", cell_size = 0.4) + facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()

  # saveRDS(temp.m, "monocle.crpc10xe.rds")
  # temp.m=readRDS("monocle.crpc10xe.rds")

  # saveRDS(temp.m, "monocle.nepce.rds")
  # temp.m=readRDS("monocle.nepce.rds")

  # find DEGs across branch 
  BEAM_res <- BEAM(temp.m, branch_point = 1, cores = 20)
  # saveRDS(BEAM_res, "monocle.BEAM_res.temp.m.rds")
  # BEAM_res=readRDS("monocle.BEAM_res.temp.m.rds")
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
    subtype_markers[[temp_subtype]]=intersect(subtype_markers[[temp_subtype]], rownames(temp@assays$SCT@counts))
  }
  rm(subtypeMarks_csv)

  sig_gene_names=intersect(row.names(subset(BEAM_res, qval<1e-4)), subtype_markers$"TF_fromLieBing")
  sig_gene_names=c(sig_gene_names, "CHGA", "SYP", "ENO2", "KLK3", "AR")
  pdf("heatmap.trajectory.temp.m.pdf")
  plot_genes_branched_heatmap(temp.m[sig_gene_names,],
                                          branch_point = 2,
                                          num_clusters = 4,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = T)
  dev.off()

  # plot specified genes
  temp_genes <- c("FEV", "ASCL1", "IRF7", "TP63", "CHGA", "SYP", "KLK3", "AR")
  pdf("dotplot.trajectory.temp.m.pdf")
  plot_genes_branched_pseudotime(temp.m[temp_genes,],
                         branch_point = 2,
                         color_by = "celltype",
                         ncol = 1)
  dev.off()

  temp@meta.data$time=rep(-1, ncol(temp))
  temp@meta.data[names(temp.m.time), ]$time=temp.m.time
  temp@meta.data$time=round(temp@meta.data$time, 1 )

  pdf("fplot.time.pdf", width=5, height=5)
  FeaturePlot(temp, features = c("time"), cols=c("#FFFFCC", "#CC3300"))
  # LightGreen 2 DarkGreen
  dev.off()
  # cols=c("grey", "blue")

  pseudotime=temp@meta.data$time
  names(pseudotime)=rownames(temp@meta.data)
  pseudotime=sort(pseudotime)
  tocheck=c(ko.up.cor.genes, wt.up.cor.genes)
  pdf("heatmap.DEG.as.pseudotime.phenotype.pdf", width=8, height=14)
  DoHeatmap( temp, features = tocheck, cells= names(pseudotime), group.by = "phenotype", slot="scale.data", assay="RNA", angle = 0 ) # 
  dev.off()


  ### find DEGs along Pseudotime 
  diff_test_res.temp.m <- differentialGeneTest(temp.m, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  # saveRDS(diff_test_res.temp.m, "diff_test_res.temp.m.NEPC.rds")
  # diff_test_res.temp.m=readRDS("diff_test_res.temp.m.NEPC.rds")
  # saveRDS(diff_test_res.temp.m, "diff_test_res.temp.m.crpc10x.rds")
  # diff_test_res.temp.m=readRDS("diff_test_res.temp.m.crpc10x.rds")
  diff_test_res =  diff_test_res.temp.m
  str(diff_test_res[,c("gene_short_name", "pval", "qval")])
  diff_test_res=subset(diff_test_res, qval < 0.1 & pval < 0.01)
  sig_gene_names=intersect(diff_test_res[,"gene_short_name"], subtype_markers$"TF_fromLieBing")
  imp.genes=c()
  imp.genes=c(imp.genes, c("CHGA", "SYP", "ENO2", "KLK3", "AR", "ASCL1"))
  imp.genes=c(imp.genes, sig_gene_names)
  imp.genes=intersect(rownames(temp@assays$SCT@counts), imp.genes)
  pdf("heatmap.trajectory.temp.m.pdf")
  plot_pseudotime_heatmap(temp.m[imp.genes,],
                  num_clusters = 3,
                  cores = 1,
                  show_rownames = T)
  dev.off()

}
