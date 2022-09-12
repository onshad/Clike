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

### read CRPC/NEPC samples
### process from raw counts
data.directory="/ ... /counts/crpc"
matrix.directory="/raw_feature_bc_matrix"
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
for (i in 1:length(inte.list))
{
  pdf(paste("temp.VlnPlot_",sample.names[i],"_quality.pdf", sep=""), width = 15, height =4)
  temp=VlnPlot(inte.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size = 0.01, group.by=NULL)
  temp[[1]]=temp[[1]]+geom_hline(yintercept=300)
  print(temp)
  dev.off()
}
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
crpc=merge(x = inte.list[[1]], y = inte.list[2:length(inte.list)] )
crpc@meta.data$fig.sample=as.character(crpc@meta.data$orig.ident)
crpc@meta.data$fig.patient=as.character(crpc@meta.data$orig.ident)
crpc@meta.data$pheno=rep("CRPC", nrow(crpc@meta.data))
crpc@meta.data$pheno[crpc@meta.data$fig.patient%in%c("QJZ")]="NEPC"
# pat2gleason=data.frame(patient=sort(unique(crpc@meta.data$fig.patient)))
# pat2gleason$patient
# pat2gleason$gleason=factor(c(NA, NA, "5", "5", "5", "5", "2", "4", "5", "5", "2", NA, "3"))
# pat2gleason
# crpc@meta.data$gleason=pat2gleason$gleason[match(crpc@meta.data$fig.patient, pat2gleason$patient)]
# crpc@meta.data$gleason=as.character(crpc@meta.data$gleason)
# table(crpc@meta.data$gleason)
str(crpc@meta.data)
# 'data.frame':   12343 obs. of  7 variables:
#  $ orig.ident  : chr  "DGY-PG" "DGY-PG" "DGY-PG" "DGY-PG" ...
#  $ nCount_RNA  : num  31899 36175 8654 4910 2793 ...
#  $ nFeature_RNA: int  5584 6899 2972 2170 1140 4956 2932 2154 1473 2011 ...
#  $ percent.mt  : num  11.87 13.48 7.78 10.75 9.42 ...
#  $ fig.sample  : chr  "DGY-PG" "DGY-PG" "DGY-PG" "DGY-PG" ...
#  $ fig.patient : chr  "DGY-PG" "DGY-PG" "DGY-PG" "DGY-PG" ...
#  $ pheno       : chr  "CRPC" "CRPC" "CRPC" "CRPC" ...
crpc=SCTransform(crpc, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting

# check most variable genes
str(crpc@assays$SCT@var.features)
# chr [1:3000] "RGS5" "COL1A1" "COL3A1" "PLA2G2A" "COL1A2" "SPP1" "SPARC" ...
# top10 <- head(VariableFeatures(crpc, assay="SCT"), 10)
# # plot variable features with and without labels
# toplot<- VariableFeaturePlot(crpc, assay = "SCT", selection.method="sctransform")
# toplot <- LabelPoints(plot = toplot, points = top10, repel = TRUE)
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()
# str(crpc)

DefaultAssay(crpc)="SCT"
crpc <- RunPCA(crpc, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(crpc,  ndims = 100)
dev.off()

### PCA-based clustering
{
  # crpc <- RunUMAP(crpc, reduction = "pca", dims = 1:25, verbose = FALSE)
  # crpc <- FindNeighbors(crpc, reduction = "pca", dims = 1:25, verbose = FALSE)
  # crpc <- FindClusters(crpc, verbose = FALSE, resolution = 0.2)
  # table(crpc@active.ident)

  # crpc.markers <- FindAllMarkers(crpc, assay="SCT")
  # saveRDS(crpc.markers, "crpc.markers.as.fig.cluster.byPCA.rds")
  # crpc.markers=readRDS("crpc.markers.as.fig.cluster.byPCA.rds")
}

### Harmony-based clustering
{
  require(harmony)
  DefaultAssay(crpc)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  crpc <- RunHarmony(crpc, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  dev.off()
  Sys.time()

  # crpc <- RunUMAP(crpc, dims = 1:25, verbose = FALSE)
  crpc <- RunUMAP(crpc, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpc <- FindNeighbors(crpc, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpc <- FindClusters(crpc, verbose = FALSE, resolution = 0.2)
  table(crpc@active.ident)
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 11897 10508  3867  2841  1996  1815  1803  1503  1345   933   861   555   511
#    13    14    15    16
#   322   316   206   198
  crpc.markers <- FindAllMarkers(crpc, assay="SCT")
  # saveRDS(crpc.markers, "crpc.markers.as.fig.cluster.byHarmony.rds")
  # crpc.markers =readRDS("crpc.markers.as.fig.cluster.byHarmony.rds")
}


crpc@meta.data$fig.cluster=crpc@meta.data$seurat_clusters
crpc <- SetIdent(crpc, value = "fig.cluster")
pdf("crpc_cluster_dimplot.pdf", width=5, height=4.5)
DimPlot(crpc, label = TRUE, group.by="fig.cluster", pt.size=0.001, label.size=5, raster=TRUE)# + NoLegend()
dev.off()
png("temp_patient_dimplot.png", width=6*100, height=5*100)
DimPlot(crpc, label = FALSE, group.by="fig.sample", pt.size=0.001)# + NoLegend()
dev.off()
png("temp_pheno_dimplot.png", width=6*100, height=5*100)
DimPlot(crpc, label = FALSE, group.by="pheno", pt.size=0.001)# + NoLegend()
dev.off()
pdf("crpc.dimplot.sample.pdf", height=4.5, width=6)
DimPlot(crpc, label = FALSE, group.by="fig.sample", label.size=5, raster=TRUE)# + NoLegend()
dev.off()


### check celltype markers' expression to make sure every cluster's celltype
markers=crpc.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[markers$avg_log2FC>log2(1.5), ]
setwd("..")
source("functions.seurat.R")
crpc@meta.data$celltype=celltypeAnno(crpc, markers)
setwd("./ ... ")
  # [1] "Luminal cell: 0.56 Epithelial cell: 0.55"
  # [3] "T cell: 1.2"
  # [5] "Luminal cell: 0.35 Epithelial cell: 0.33"
  # [12] "MDSCfunction_profRen: 3.83"
# cluster 11 is epithelial cells (KRT5+, TACSTD2+, KRT19+), and not from bladder tissue ( EMA-, CDH2-,  UGP85-)
crpc@meta.data$celltype[crpc@meta.data$fig.cluster=="11"]="Epithelial cell"

crpc@meta.data$celltype[crpc@meta.data$fig.cluster=="0"]="Epithelial cell"
crpc@meta.data$celltype[crpc@meta.data$fig.cluster=="4"]="Epithelial cell"
crpc@meta.data$celltype[crpc@meta.data$fig.cluster=="7"]="Epithelial cell"
# plot with cluster names
pdf("crpc_celltype_dimplot.pdf", width=6, height=4.5)
DimPlot(crpc, label = TRUE, group.by="celltype", pt.size=0.001, label.size=5, raster=TRUE)# + NoLegend()
dev.off()


# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=crpc.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.crpc.cluster","csv",sep="."), row.name=T)


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
# check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3")
# check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
# check.genes=c("HOXB-AS1", "HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXB9", "HOXB-AS3", "SKAP1")
# check.genes=c("INHBA", "VCAN", "ZFHX4")

# crpc=SetIdent(crpc, value = "fig.cluster")
crpc=SetIdent(crpc, value = "celltype")
DefaultAssay(crpc)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpc, features = check.genes)
dev.off()
DefaultAssay(crpc)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(crpc, features = check.genes, pt.size = 0, group.by="celltype", ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(crpc, "crpc.rds")
# crpc=readRDS("crpc.rds")


imp.genes=c()
imp.genes=c("CD3D", "CD3E", "CD3G",    "CD4", "CD14", "CD163",   "FTL", "CD55", "CD44",    "KRT8", "KRT18", "AR",     "TPSB2", "TPSAB1", "MS4A2",     "ENG", "VWF", "PECAM1",  "ACTA2", "VIM", "CNN1",    "COL1A1", "FBLN1", "FN1",       "IGHG4", "MS4A1", "MZB1" )
# DES: myofibroblast
imp.genes=unique(imp.genes)
imp.genes=imp.genes[length(imp.genes):1]# reverse
# delete celltype with less than 100 cells
crpc=subset(crpc, subset=celltype!="Myocyte")
crpc@meta.data$fig.celltype=factor(crpc@meta.data$celltype, levels=c("T cell", "Macrophage", "Epithelial cell", "Mast cell", "Endothelial cell", "Myofibroblast", "Fibroblast", "B/Plasma cell"))
pdf("crpc_dotmap_celltype_markers.pdf", width = 18, height = 3);
DotPlot(crpc, features = imp.genes, cols = c('white', 'blue'), dot.scale = 6, group.by="fig.celltype");
dev.off();




### further analysis focused on epithelial cells
table(crpc@meta.data$celltype)
   # B/Plasma cell Endothelial cell  Epithelial cell       Fibroblast
   #           145              850              133              292
   #  Luminal cell       Macrophage        Mast cell    Myofibroblast
   #          6482             2135              140              999
   #        T cell
   #          1167
epitypes=c("Epithelial cell", "Luminal cell", "Basal cell")
crpce=subset(crpc, subset=(celltype%in%epitypes))
# str(crpce)
crpce@commands=list()


dim(crpce)
# [1] 19811  6615
DefaultAssay(crpce)="SCT"
str(crpce@assays$SCT@var.features)
# chr [1:3000] "RGS5" "COL1A1" "COL3A1" "PLA2G2A" "COL1A2" "SPP1" "SPARC" ...
# crpce <- FindVariableFeatures(crpce, assay="SCT", nfeatures = 3000)
# str(crpce@assays$SCT@var.features)
# top10 <- head(VariableFeatures(crpce, assay="SCT"), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(crpce, assay = "SCT")
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# toplot=plot1 + plot2
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()

crpce=SCTransform(crpce, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
crpce <- RunPCA(crpce, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(crpce,  ndims = 100)
dev.off()

### clustering basd on PCA
{
  # DefaultAssay(crpce)="SCT"
  # crpce <- RunUMAP(crpce, dims = 1:20, verbose = FALSE)
  # crpce <- FindNeighbors(crpce, dims = 1:20, verbose = FALSE)

  # DefaultAssay(crpce)="SCT" 
  # crpce <- FindClusters(crpce, verbose = FALSE, resolution = 0.2)
  # table(crpce@meta.data$seurat_clusters)
  # crpce.markers <- FindAllMarkers(crpce, assay="SCT")
  # saveRDS(crpce.markers, "crpce.markers.as.fig.cluster.byPCA.rds")
  # crpce.markers =readRDS("crpce.markers.as.fig.cluster.byPCA.rds")

}

### clustering basd on Harmony
{
  require(harmony)
  DefaultAssay(crpce)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  crpce <- RunHarmony(crpce, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
  dev.off()
  Sys.time()

  # crpc <- RunUMAP(crpc, dims = 1:25, verbose = FALSE)
  crpce <- RunUMAP(crpce, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpce <- FindNeighbors(crpce, reduction = "harmony", dims = 1:25, verbose = FALSE)
  crpce <- FindClusters(crpce, verbose = FALSE, resolution = 0.2)
  table(crpce@meta.data$seurat_clusters)
  #    0    1    2    3    4    5    6    7    8    9   10
  # 4069 3853 2325 1659 1446 1289  499  229   88   76   43
  crpce.markers <- FindAllMarkers(crpce, assay="SCT")
  saveRDS(crpce.markers, "crpce.markers.as.fig.cluster.byHarmony.rds")
  # crpce.markers =readRDS("crpce.markers.as.fig.cluster.byHarmony.rds")
}

markers=crpce.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.crpce.cluster","csv",sep="."), row.name=T)

crpce@meta.data$fig.cluster=crpce@meta.data$seurat_clusters
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(crpce, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
pdf("crpce.dimplot.patient.pdf", height=4.5, width=5)
DimPlot(crpce, label = FALSE, group.by="fig.patient", label.size=5, raster=TRUE)
dev.off()
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(crpce, label = FALSE, group.by="gleason")# + NoLegend()
dev.off()
png("dimplot.pheno.png", height=5*100, width=5*100)
DimPlot(crpce, label = FALSE, group.by="pheno")# + NoLegend()
dev.off()


### analysis of crpctasis STOP here, because the small sample number (206 cells)
# saveRDS(crpce, "crpce.rds")
# crpce=readRDS("crpce.rds")

# cluster 4 is CD3E+CD52+ (T cell markers), but it is also AR+ KLK3+
# cluster 7 is ACTA2+CD34+COL1A1+ (Fibroblast markers), but it is also AR+ KLK3+
crpce=subset(crpce, cells = colnames(crpce)[!crpce@meta.data$fig.cluster%in%c("4", "7")])
table(crpce@meta.data$fig.cluster)
#    0    1    2    3    4    5    6    7
# 1972 1755 1695  795    0  132   70    0

# plot previously identified basal/cycle signatures
# all 203 basal markers
# crpce <- ScaleData(crpce, features = rownames(crpce), assay="RNA" )
basal.posi.genes=read.table("../basal.in.all.cell.markers.posi.genes.txt", head=FALSE, stringsAsFactors=FALSE)$V1
# all 126 cycle markers
cycle.posi.genes=read.table("../cycle.in.all.cell.markers.posi.genes.txt", head=FALSE, stringsAsFactors=FALSE)$V1
basal.posi.genes=basal.posi.genes[basal.posi.genes%in%rownames(crpce@assays$SCT@scale.data)]
cycle.posi.genes=cycle.posi.genes[cycle.posi.genes%in%rownames(crpce@assays$SCT@scale.data)]
crpce@meta.data$basal.sig=colMeans(crpce@assays$SCT@scale.data[basal.posi.genes, ])
DefaultAssay(crpce)="SCT"
pdf("temp.basal.sig.featureplot.pdf", width=5, height=5)
FeaturePlot(crpce, features = "basal.sig")
dev.off()
crpce@meta.data$cycle.sig=colMeans(crpce@assays$SCT@scale.data[cycle.posi.genes, ])
pdf("temp.cycle.sig.featureplot.pdf", width=5, height=5)
FeaturePlot(crpce, features = "cycle.sig")
dev.off()
DefaultAssay(crpce)="SCT"
pdf("temp.VlnPlot.pdf", width = 11, height =4)
tocheck=c("basal.sig", "cycle.sig")
temp=VlnPlot(crpce, features = tocheck, ncol = 2,  pt.size = 0, group.by="fig.cluster", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1", "UPK2", "AR")
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # check genes 
# check.genes=c("CCL2", "CD74")

# crpce=subset(crpce, subset=(pheno!="NEPC"))
crpce <- SetIdent(crpce, value = "fig.cluster")
DefaultAssay(crpce)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpce, features = check.genes)
dev.off()
DefaultAssay(crpce)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(crpce, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# indi.plot(crpce, tocheck="fig.patient")
# indi.plot(crpce, tocheck="fig.cluster")

### re-cluster basal & typeC cells
{
  crpcb=subset(crpce, subset=(fig.cluster%in%"2"))
  dim(crpcb)
  # [1] 20066  2325
  DefaultAssay(crpcb)="SCT"
  str(crpcb@assays$SCT@var.features)
  # chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "MMP7" "KRT13" "LCN2" ...
  # crpcb <- FindVariableFeatures(crpcb, assay="SCT", nfeatures = 3000)
  # str(crpcb@assays$SCT@var.features)
  # top10 <- head(VariableFeatures(crpcb, assay="SCT"), 10)
  # # plot variable features with and without labels
  # plot1 <- VariableFeaturePlot(crpcb, assay = "SCT")
  # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # toplot=plot1 + plot2
  # png("temp.png", width=10*100, height=5*100)
  # print(toplot)
  # dev.off()

  crpcb=SCTransform(crpcb, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  crpcb <- RunPCA(crpcb, assay="SCT", verbose = FALSE, npcs = 100)
  ElbowPlot(crpcb,  ndims = 100)
  dev.off()

  ### clustering basd on PCA
  {
    # DefaultAssay(crpcb)="SCT"
    # crpcb <- RunUMAP(crpcb, dims = 1:20, verbose = FALSE)
    # crpcb <- FindNeighbors(crpcb, dims = 1:20, verbose = FALSE)

    # DefaultAssay(crpcb)="SCT" 
    # crpcb <- FindClusters(crpcb, verbose = FALSE, resolution = 0.2)
    # table(crpcb@meta.data$seurat_clusters)
    # crpcb.markers <- FindAllMarkers(crpcb, assay="SCT")
    # saveRDS(crpcb.markers, "crpcb.markers.as.fig.cluster.byPCA.rds")
    # crpcb.markers =readRDS("crpcb.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(crpcb)
    # [1] "SCT"
    Sys.time()
    crpcb <- RunHarmony(crpcb, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    dev.off()
    Sys.time()

    # crpc <- RunUMAP(crpc, dims = 1:25, verbose = FALSE)
    crpcb <- RunUMAP(crpcb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    crpcb <- FindNeighbors(crpcb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    crpcb <- FindClusters(crpcb, verbose = FALSE, resolution = 0.2)
    table(crpcb@meta.data$seurat_clusters)
    #   0   1   2   3   4
    # 894 652 577 171  31
    # crpcb.markers <- FindAllMarkers(crpcb, assay="SCT")
    # saveRDS(crpcb.markers, "crpcb.markers.as.fig.cluster.byHarmony.rds")
    # crpcb.markers =readRDS("crpcb.markers.as.fig.cluster.byHarmony.rds")
  }

  crpcb@meta.data$fig.cluster=as.character(crpcb@meta.data$seurat_clusters)
  png("dimplot.cluster.png", height=5*100, width=5*100)
  DimPlot(crpcb, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  dev.off()
  png("dimplot.patient.png", height=5*100, width=5*100)
  DimPlot(crpcb, label = TRUE, group.by="fig.patient")
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  DimPlot(crpcb, label = TRUE, group.by="gleason")
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  DimPlot(crpcb, label = TRUE, group.by="pheno")
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
  # check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
  # check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3")

  crpcb <- SetIdent(crpcb, value = "fig.cluster")
  DefaultAssay(crpcb)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(crpcb, features = check.genes)
  dev.off()
  DefaultAssay(crpcb)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=8*100)
  temp=VlnPlot(crpcb, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  ### name clusters
  current.cluster.ids <- c(0:4)
  new.cluster.ids <- c("tC20", "tB21", "tL22", "tC23", "tC24")
  # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
  crpcb@meta.data$epitype <- plyr::mapvalues(x = as.character(crpcb@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  pdf("dimplot.epitype.pdf", height=4.5, width=5)
  DimPlot(crpcb, label = TRUE, group.by="epitype", label.size=5, raster=TRUE)
  dev.off()

  btype=crpcb@meta.data$epitype
  names(btype)=rownames(crpcb@meta.data)
  table(btype)
  # btype
  # tB21 tC20 tC23 tC24 tL22
  #  652  894  171   31  577

  # indi.plot(crpcb, tocheck="fig.patient")
  # indi.plot(crpcb, tocheck="fig.cluster")

  # saveRDS(crpcb, "crpcb.rds")
  # crpcb=readRDS("crpcb.rds")
}

### epitype annotation
current.cluster.ids <- c(0:3, 5:6)
new.cluster.ids <- c("crL0", "crL1", "crL2", "crCy3", "crNE5", "crC6")
# tL: tumor luminal, tB: tumor basal, tCy: tumor CellCycle
crpce@meta.data$epitype <- plyr::mapvalues(x = as.character(crpce@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot
# crpce@meta.data[names(btype), "epitype"]=btype
pdf("crpce.dimplot.epitype.pdf", height=4.5, width=5)
DimPlot(crpce, label = TRUE, group.by="epitype", label.size=5, raster=TRUE)
dev.off()

### celltype annotation
current.cluster.ids <- unique(crpce@meta.data$epitype)
current.cluster.ids=sort(current.cluster.ids)
current.cluster.ids
#  [1] "crC6"  "crCy3" "crL0"  "crL1"  "crL2"  "crNE5"
new.cluster.ids <- c("Clike", "CellCycle", "Luminal", "Luminal", "Luminal", "NE")
# tLum: tumor luminal, tBI: tumor basal/intermediate, tCyc: tumor CellCycle
crpce@meta.data$celltype <- plyr::mapvalues(x = as.character(crpce@meta.data$epitype), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot
### f2.g
png("dimplot.celltype.png", height=5*100, width=5*100)
DimPlot(crpce, label = TRUE, group.by="celltype", label.size=5)
dev.off()
pdf("f2.g.crpce.celltype.dimplot.pdf", height=4.5, width=5)
temp.color=color[names(color)%in%unique(crpce@meta.data$celltype)]
temp.oorder=oorder[oorder%in%unique(crpce@meta.data$celltype)]
DimPlot(crpce, label = TRUE, group.by="celltype", label.size=5, raster=TRUE, cols=temp.color, order=temp.oorder)
dev.off()


### f2h
imp.genes = c( "AR", "KRT8",    "KRT5",    "CHGA","SYP","ENO2",    "CYP1A1", "ARL14",    "KRT4", "PSCA" )
# "CD55", "FGFR3",
imp.genes=c(  "KRT4", "PSCA",     "KRT5", "AR", "KRT8",   "CHGA","SYP","ENO2")
imp.genes=c(  "KRT4", "TACSTD2", "PSCA",     "KRT5", "AR", "KRT8",   "CHGA","SYP","ENO2",   "NCAM1", "ACSL4")

DefaultAssay(crpce)="SCT"
temp.color=color[names(color)%in%unique(crpce@meta.data$celltype)]
temp.oorder=oorder[oorder%in%unique(crpce@meta.data$celltype)]
crpce@meta.data$fig.celltype=factor(crpce@meta.data$celltype, levels=temp.oorder)
crpce <- SetIdent(crpce, value = "fig.celltype")
pdf("crpce.VlnPlot.pdf", width=20, height=2)
temp=VlnPlot(crpce, features = imp.genes, pt.size = 0, group.by=NULL, ncol = 11, cols=temp.color)
for(i in 1:length(imp.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 5, colour = "black", shape = 95)
}
print(temp)
dev.off()


crpce <- SetIdent(crpce, value = "celltype")
DefaultAssay(crpce)="SCT"

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

DefaultAssay(crpce)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpce, features = check.genes)
dev.off()
DefaultAssay(crpce)="SCT"
png("temp.VlnPlot.png", width=10*100, height=12*100)
temp=VlnPlot(crpce, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# saveRDS(crpce, "crpce.rds")
# crpce=readRDS("crpce.rds")

temp.marker=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
temp.marker=c("UPK1A", "CD55", "TGFBR3", "FGFR3")
temp.marker=c("TP63", "TRIM29", "ITGB4", "KRT5", "KRT14", "AR", "KLK3", "DPP4", "ESR1", "ESR2", "NR3C1", "PGR", "DUSP6", "SPRY4", "ETV4", "GDF15", "EGR1", "ETV5", "C3orf52", "UPP1", "DUSP5", "TRIB3", "MYEOV", "MAFF", "PHLDA1") # Double-neg marker, Cancer Cell 2017
temp=subset(crpce, downsample=1000)
# temp=NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA" )
# temp=ScaleData(temp, features = rownames(temp), assay="RNA" )
levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
pdf(paste("temp.heatmap.pdf", sep=""), width=(length(unique(crpce@active.ident))/4+5), height=length(temp.marker)/10+3)
temp.heatmap<-DoHeatmap(temp , features = temp.marker, slot="scale.data", assay="SCT", angle = 0, size=3) 
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

# crpce.markers=readRDS("crpce.markers.celltype.rds")
crpce.markers[crpce.markers$cluster=="Clike" & crpce.markers$gene=="TGFBR3", ]
#              p_val avg_log2FC pct.1 pct.2   p_val_adj cluster   gene
# TGFBR3 2.03857e-07  0.2556965   0.3 0.104 0.003816204   Clike TGFBR3
crpce.markers[crpce.markers$cluster=="Clike" & crpce.markers$gene=="FGFR3", ]
             # p_val avg_log2FC pct.1 pct.2    p_val_adj cluster  gene
# FGFR3 5.369109e-49  0.6667128 0.457 0.053 1.005097e-44   Clike FGFR3
 
### check specific marker across all cell types
# crpc=readRDS("crpc.rds")
check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
clike.cells=rownames(crpce@meta.data)[crpce@meta.data$celltype=="Clike"]
crpc@meta.data$celltype[rownames(crpc@meta.data)%in%clike.cells]="Clike"
crpc <- SetIdent(crpc, value = "celltype")
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(crpc, features = check.genes)
dev.off()
DefaultAssay(crpc)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(crpc, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()



### Not employed --->>>
  ### CCA to identify epithelial cell types WithiN luminal tumor cells
  # norme=readRDS("norme.rds")
  crpcl=subset(crpce, subset=celltype=="Luminal")
  anchors <- FindTransferAnchors(reference = norme, query = crpcl, 
      dims = 1:30, normalization.method="LogNormalize", reference.assay="SCT")
  predictions <- TransferData(anchorset = anchors, refdata = norme$celltype, 
      dims = 1:30)
  crpcl@meta.data$predtype=predictions$predicted.id 
  png("dimplot.predtype.png", height=5*100, width=5*100)
  DimPlot(crpcl, label = TRUE, group.by="predtype")
  dev.off()
### Not employed <<<---


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(crpce@meta.data$celltype)
# CellCycle     Clike   Luminal        NE
#       795        70      5422       132
table(crpce@meta.data$fig.patient)
# DGY-PG    HYL    QJZ    XYM
#   3489   1195    169   1566
crpce@meta.data$fig.patient=factor(crpce@meta.data$fig.patient, levels=c("QJZ", "DGY-PG", "HYL", "XYM"))
pct.matrix=table(crpce@meta.data$fig.patient, crpce@meta.data$celltype)
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
write.csv(pct.matrix, file="part4__crpce__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)

### boxplot Clike percentage for CRPC+NEPC vs mHSPC vs Primary and incidental aggressive vs indolent
{
  # give sample ID as phenotype + number
  {
    crpce@meta.data$fig.pheno=factor(crpce@meta.data$pheno, levels=c("CRPC", "NEPC"))
    pheno.sample=crpce@meta.data[, c("fig.pheno", "fig.sample")]
    pheno.sample$fig.sample=as.character(pheno.sample$fig.sample)
    # str( pheno.sample)
    # 'data.frame':   118341 obs. of  2 variables:
    #  $ fig.pheno : Factor w/ 8 levels "Normal prostate",..: 1 1 1 1 1 1  ...
    #  $ fig.sample: chr  "waizhouyouCBH" "waizhouyouCBH" "waizhouyouCBH" ...
    pheno.sample=pheno.sample[order(pheno.sample$fig.pheno, pheno.sample$fig.sample, decreasing=c(FALSE, FALSE)), ]
    unique.sample=unique(pheno.sample$fig.sample)
    pheno.sample$fig.sample=factor(pheno.sample$fig.sample, levels=unique.sample)
    pheno.sample$fig.sample.id=as.numeric(pheno.sample$fig.sample)
    pheno.sample$pheno.sample=paste(as.character(pheno.sample$fig.pheno), pheno.sample$fig.sample.id)
    unique.pheno.sample=unique(pheno.sample$pheno.sample)
    pheno.sample$pheno.sample=factor(pheno.sample$pheno.sample, levels=unique.pheno.sample)
    crpce@meta.data$pheno.sample=pheno.sample[rownames(crpce@meta.data), ]$pheno.sample
  }

  # get percentage table of phenotype+sampleID vs celltype
  celltype.sta=table(crpce@meta.data$pheno.sample, crpce@meta.data$celltype)
  celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
  print(unique(rownames(celltype.sta)))
  celltype.sta=celltype.sta[c("CellCycle", "Clike", "Luminal", "NE" ), ]# 
  barN=ncol(celltype.sta)
  n_top = nrow(celltype.sta) 
  pdf(paste("crpce.barplot.celltypeInSample.pdf",sep=""))
  # the width between paper edge and plot ( down,left,up,right : four directions)
  par(mar = c(14, 4, 1, 8), xpd=T)
  # coll=DiscretePalette(n_top, palette = "polychrome")
  require("scales")
  # get 4 colors, then delete the luminal and cellcycle's color 
  coll=hue_pal()(n_top)
  # coll=coll[c(-3, -1)]
  barplot( celltype.sta, xlab="", ylab="Percentage", col=coll, names.arg=colnames(celltype.sta), las=2, cex.names=1) # , xaxt="n"
  # facts=barN+1
  end_point = barN #0.5 + n_top *facts-1
  # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
  legend( barN+1, 15, rownames(celltype.sta), fill=coll );
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
  dev.off()

}



### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(crpce@meta.data$epitype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(crpce@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(crpce@meta.data$fig.patient, crpce@meta.data$epitype)
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
write.csv(pct.matrix, file="part4__crpce__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


# saveRDS(crpce, "crpce.rds")
# crpce=readRDS("crpce.rds")


### celltype-specific DEG identification
table(crpce@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   12893     925
DefaultAssay(crpce)="SCT"
### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(crpce@meta.data$celltype)
celltype=sort(celltype)
crpce = SetIdent(crpce, value = "celltype")
for(i in celltype)
{
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        markerlist[[temp.name]]=FindMarkers(crpce, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
      }
  }
}
# saveRDS(markerlist, "crpce.celltype.pairwise.rds")
# markerlist=readRDS("crpce.celltype.pairwise.rds")

### cannonical process to find specific markers of Clike cells
{
  markerlist=readRDS("crpce.celltype.pairwise.rds")
  sp.markerlist=list()
  celltype=unique(crpce@meta.data$celltype)
  celltype=sort(celltype)
  for(i in celltype)
  {
    temp.genes=rownames(crpce@assays$SCT@data)
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          cur.list=markerlist[[temp.name]]
          cur.list=cur.list[order(cur.list$avg_log2FC, decreasing=TRUE), ]
          # cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(2) & cur.list$pct.1>0.3 & cur.list$pct.2<0.1,]) # mature version
          # cur.genes=rownames( cur.list[cur.list$avg_log2FC>log2(1.5),] ) # best intersect between prime, crpce, triple
          cur.genes=rownames( cur.list[cur.list$avg_log2FC>log2(1.5) & cur.list$pct.1>0.4,] )
          temp.genes=intersect(temp.genes, cur.genes)
        }
    }
    sp.markerlist[[i]]=temp.genes
  }
  str(sp.markerlist)
  sp.markerlist[["Luminal"]]=c(sp.markerlist[["Luminal"]], "KLK3", "AR")
  sp.markerlist[["TypeC"]]=c(sp.markerlist[["TypeC"]], "KRT4", "TACSTD2", "PSCA")

  temp=subset(crpce, downsample=300)
  sp.markerlist=sp.markerlist[c("Clike", "Luminal", "NE", "CellCycle")]
  tocheck=Reduce(union, sp.markerlist)
  temp <- SetIdent(temp, value = "celltype")
  levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
  pdf(paste("heatmap.crpce.typeC.pdf", sep=""),  width=10, height=length(tocheck)/3)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

  crpce.clike.markers=sp.markerlist$Clike
  str(crpce.clike.markers)
}

intersect(crpce.clike.markers, prime.clike.markers)

best=c("CD55", "IL1RN", "IL1A", "UPK1A", "UPK1B", "PIM1", "PMAIP1")
tocheck=best

### DEG identification, across all clusters
crpce <- SetIdent(crpce, value = "celltype")
crpce.markers <- FindAllMarkers(crpce, slot="data", assay="SCT")
# saveRDS(crpce.markers, "crpce.markers.celltype.rds")
# crpce.markers=readRDS("crpce.markers.celltype.rds")

### plot heatmap
# crpce.markers=readRDS("crpce.indivSCT.markers.as.fig.cluster.rds")
# crpce <- NormalizeData(crpce, normalization.method = "LogNormalize", scale.factor = 10000, assay="SCT" )
# crpce <- ScaleData(crpce, features = rownames(crpce), assay="RNA")

# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=crpce.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.crpce.cluster","csv",sep="."), row.name=T)

enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)


# heatmap top 10 DEGs for every cluster
markers=crpce.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(2) & markers$pct.1>0.3 & markers$pct.2<0.1, ]
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
top10$cluster=factor(top10$cluster, levels=c("Clike", "Luminal", "NE", "CellCycle"))
top10=top10[order(top10$cluster), ]
temp=subset(crpce, downsample=300)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
tocheck=c(top10$gene, "KLK3", "AR")
tocheck=c(tocheck, "KRT4", "TACSTD2", "PSCA")
pdf(paste("heatmap.crpce.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=length(tocheck)/8)
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### plot heatmap for special categories of DEGs
setwd("..")
subtype_markers=subtypeMarker(crpce)
setwd("./ ... ")
crpce.markers=readRDS("crpce.markers.celltype.rds")
markers=crpce.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) , ]#& markers$pct.1>0.2
# for TFs, change log2(1.5) to log(1.2)
# crpce <- ScaleData(crpce, features = rownames(crpce), assay="RNA" )
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
    DefaultAssay(crpce)="SCT"
    pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
    temp=subset(crpce, downsample=300)
    temp <- SetIdent(temp, value = "celltype")
    levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
    temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()
    DefaultAssay(crpce)="SCT"
  }
}


### vlnplot for special categories of DEGs
setwd("..")
subtype_markers=subtypeMarker(crpce)
setwd("./ ... ")
crpce.markers=readRDS("crpce.markers.celltype.rds")
markers=crpce.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.4, ]#

MHCII=c("HLA-DRB3", "HLA-DRA", "HLA-DRB1", "HLA-DPB1", "HLA-DPA1", "HLA-DQA1", "HLA-DQB1")
crpce@meta.data$MHCII=colMeans(crpce@assays$RNA@data[MHCII[MHCII%in%rownames(crpce@assays$RNA@data)], ])


tocheck="chemokine.ligand" # "inhibi.check.point.ligand" "chemokine.ligand"
tocheck=subtype_markers[[tocheck]]
temp.markers=markers[markers$gene%in%tocheck, ]
top10 <- temp.markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
top10=top10[top10$cluster=="Clike", ]
top10$cluster=factor(top10$cluster, levels=c("Clike", "Luminal", "NE", "CellCycle"))
top10=top10[order(top10$cluster), ]
check.genes=top10$gene
check.genes=c(check.genes, "MHCII")
DefaultAssay(crpce)="RNA"
pdf("VlnPlot.chemokines_celltype.pdf", width = 10, height =10)
temp.crpce=crpce
temp.crpce@meta.data$celltype=factor(temp.crpce@meta.data$celltype, levels=c("Clike", "Luminal", "NE", "CellCycle"))
temp=VlnPlot(temp.crpce, features = check.genes, pt.size = 0, group.by="celltype", ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()
DefaultAssay(crpce)="SCT"


### get JAK-STAT, STEMNESS, and EMT genesets
{
  library(msigdbr)
  m_df.kegg = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  m_t2g.kegg = m_df.kegg %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  jakstat.bak=unique(m_t2g.kegg$gene_symbol[m_t2g.kegg$gs_name=="KEGG_JAK_STAT_SIGNALING_PATHWAY"])
  m_df.cgp = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
  m_t2g.cgp = m_df.cgp %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  stem.bak=unique(m_t2g.cgp$gene_symbol[m_t2g.cgp$gs_name=="MALTA_CURATED_STEMNESS_MARKERS"])
  m_df.hall = msigdbr(species = "Homo sapiens", category = "H", subcategory = "")
  m_t2g.hall = m_df.hall %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  emt.bak=unique(m_t2g.hall$gene_symbol[m_t2g.hall$gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])
}

# find differntiall jak-stat genes
tumor.markers=readRDS("crpce.markers.celltype.rds")
markers=tumor.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[abs(markers$avg_log2FC)>log2(1.5), ]
markers$foldChange=2^(markers$avg_log2FC)
# markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
markers=markers[order(markers$foldChange, decreasing=TRUE), ]
ordering.genes=unique(markers$gene)

### check average value of genesets
# crpce <- ScaleData(crpce, features = rownames(crpce), assay="RNA" )
jakstat=jakstat.bak
jakstat=jakstat[jakstat%in%rownames(crpce@assays$SCT@scale.data)]
crpce@meta.data$jakstat.sig=colMeans(crpce@assays$SCT@scale.data[jakstat, ])

stem=stem.bak
stem=stem[stem%in%rownames(crpce@assays$SCT@scale.data)]
crpce@meta.data$stem.sig=colMeans(crpce@assays$SCT@scale.data[stem, ])

emt=emt.bak
emt=emt[emt%in%rownames(crpce@assays$SCT@scale.data)]
crpce@meta.data$emt.sig=colMeans(crpce@assays$SCT@scale.data[emt, ])


check.genes=c("emt.sig", "jakstat.sig", "stem.sig") # jakstat.sig   stem.sig   emt.sig   c("SOX2", "RB1", "TP53", "PTEN", "IL6")
temp.crpce=crpce
temp.crpce@meta.data$celltype=factor(temp.crpce@meta.data$celltype, levels=c("Clike", "Luminal", "NE", "CellCycle"))
temp.crpce <- SetIdent(temp.crpce, value = "celltype")
DefaultAssay(temp.crpce)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(temp.crpce, features = check.genes)
dev.off()
DefaultAssay(temp.crpce)="SCT"
pdf("temp.VlnPlot.pdf", width=12, height=5)
temp=VlnPlot(temp.crpce, features = check.genes, pt.size = 0, group.by="celltype", ncol = 3)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()
DefaultAssay(temp.crpce)="SCT"


#
tocheck=c( "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "FGFR3", "TGFBR3")
crpce<- SetIdent(crpce, value = "celltype")
temp=subset(crpce, downsample=300)
levels(temp)=c("Clike", "Luminal", "NE", "CellCycle")
pdf(paste("heatmap.CRPC.Clike.pdf", sep=""),  width=5, height=5)#height=length(tocheck)/3
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)#+ scale_fill_gradientn(colors = c("green", "red"))
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


# enrichment of Clike compared with all other cells
markers=crpce.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)

# enrichment of Clike compared with TypeC
markerlist=readRDS("crpce.celltype.pairwise.rds")
markers=markerlist[["Clike__TypeC"]]
markers=markers[order(markers$avg_log2FC,decreasing = TRUE), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(rownames(markers), topNterm=20)


# enrichment of Clike compared with Luminal
markerlist=readRDS("crpce.celltype.pairwise.rds")
markers=markerlist[["Clike__Luminal"]]
markers=markers[order(markers$avg_log2FC,decreasing = TRUE), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & (markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(rownames(markers), topNterm=20)



### deconvolution for validation in bulk data
{
  ### add crpce's cell types to crpc's cell types
  crpc=readRDS("crpc.rds")
  crpce=readRDS("crpce.rds")

  png("temp_crpc_celltype_dimplot.png", width=5*100, height=5*100)
  DimPlot(crpc, label = TRUE, group.by="celltype") + NoLegend()
  dev.off()

  png("temp_crpce_celltype_dimplot.png", width=5*100, height=5*100)
  DimPlot(crpce, label = TRUE, group.by="celltype") + NoLegend()
  dev.off()

  crpc@meta.data$celltype[crpc@meta.data$celltype=="B cell Lymph node"]="B/Plasma cell"
  temp.epi.cellnames=colnames(crpc)[crpc@meta.data$celltype%in%c("Basal cell", "Epithelial cell", "Luminal cell")]
  temp.epi.cellnames=temp.epi.cellnames[!temp.epi.cellnames%in%colnames(crpce)]
  temp.epi.cellnames=colnames(crpc)[!colnames(crpc)%in%temp.epi.cellnames]
  crpc=subset(crpc, cells = temp.epi.cellnames)
  epitype=crpce@meta.data$celltype
  names(epitype)=colnames(crpce)
  crpc@meta.data[names(epitype), "celltype"]=epitype
  table(crpc@meta.data$celltype)
   # B/Plasma cell        CellCycle            Clike Endothelial cell
   #           145              795               70              850
   #    Fibroblast          Luminal       Macrophage        Mast cell
   #           292             5422             2135              140
   # Myofibroblast               NE           T cell
   #           999              132             1167


  DefaultAssay(crpc)="SCT"
  ### pairwise customized cluster's DEGs
  markerlist=list()
  celltype=unique(crpc@meta.data$celltype)
  celltype=sort(celltype)
  crpc = SetIdent(crpc, value = "celltype")
  for(i in celltype)
  {
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          markerlist[[temp.name]]=FindMarkers(crpc, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
        }
    }
  }
  saveRDS(markerlist, "crpc.celltype.pairwise.rds")
  

  # @data
  markerlist=readRDS("crpc.celltype.pairwise.rds")
  # @scale.data
  # markerlist=readRDS("markerlist.crpc.pairwise.initiating.types.scaled.data.rds")
  ######
  celltypes=sort(unique(crpc@meta.data$celltype))
  markerlist.filter=list()
  for(j in celltypes)
  {
    print("celltype:")
    print(j)
    print("------------>>>>>>>>>>>")
    sp.markers=rownames(crpc@assays$SCT@data)
    sp.markerlist="1"
    for(i in names(markerlist))
    {
      temp=strsplit(i, "_")[[1]][1]
      if(temp==j)
      {
        markers=markerlist[[i]]
        markers=markers[markers$p_val_adj<0.05, ]
        markers=markers[order(markers$avg_log2FC, decreasing = TRUE), ]
        markers=markers[markers$avg_log2FC>log2(1.5), ] # mature version
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
  write.table(markerlist.filter[["Clike"]], file="marker.Clike.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["TypeC"]], file="marker.TypeC.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["Basal"]], file="marker.Basal.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["Luminal"]], file="marker.Luminal.txt", col.name=F, row.name=F, quote=F, sep='\t')
  for(i in names(markerlist.filter))
  {
    markerlist.filter[[i]]=markerlist.filter[[i]][1:4]
    markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
  }
  markers=Reduce(union, markerlist.filter)
  str(markerlist.filter)
  str(markers)

  pdf(paste("temp.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/9 )
  temp.heatmap<-DoHeatmap( subset(crpc, downsample=300), features = markers, group.by = "celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

  # @data-based markers
  # write.table(markers, file="marker.union.crpc.pairwise.initiating.types.txt", col.name=F, row.name=F, quote=F, sep='\t')
  # @scale.data-based markers:
  write.table(markers, file="markerlist.crpc.pairwise.celltype.txt", col.name=F, row.name=F, quote=F, sep='\t')

  # @data-based markers
  # markers=read.table("marker.union.crpc.pairwise.initiating.types.txt", sep='\t')$V1
  # @scale.data-based markers:
  markers=read.table("markerlist.crpc.pairwise.celltype.txt", sep='\t')$V1
  # make Cibersort's Signature matrix
  celltypes=sort(unique(crpc@meta.data$celltype))
 # [1] "B/Plasma cell"    "Basal"            "Clike"            "Endothelial cell"
 # [5] "Fibroblast"       "Luminal"          "Macrophage"       "Mast cell"
 # [9] "Myofibroblast"    "Neutrophil"       "T cell"           "TypeC"
  ### create signature matrix "sigm"
  sigm=crpc@assays$SCT@scale.data[markers, ]
  # sigm=crpc@assays$SCT@data[markers, ]
  colnames(sigm)=crpc@meta.data$celltype
  sigm=t(sigm)

  sigm=aggregate(sigm, by = list(rownames(sigm)), FUN = mean, na.rm = TRUE)
  sigm[1:3, 1:3]
  rownames(sigm)=sigm$"Group.1"
  sigm=sigm[, c(-1)]
  sigm=as.matrix(sigm)
  sigm=t(sigm)
  sigm[1:3, 1:3]
  # saveRDS(sigm, "crpc.celltype.average.sig.matrix.rds")
  saveRDS(sigm, "crpc.celltype.average.sig.matrix.scale.data.rds")
  # saveRDS(sigm, "crpc.celltype.average.sig.matrix.2fc.dataslot.rds")

  ### use FARDEEP to deconvolution
  ### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522071/
  # Rscript stage.R ../ ... /crpc.celltype.average.sig.matrix.scale.data.rds fardeep no no Clike CRPC.Clike;Rscript survival2.R ../ ... /crpc.celltype.average.sig.matrix.scale.data.rds fardeep no no Clike CRPC.Clike
}


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
crpce@meta.data$for_temp_color=as.numeric(factor(crpce@meta.data$celltype, levels=c("Clike", "Luminal", "NE",  "CellCycle")))
run_qusage_heatmap.seurat3(crpce, nm = 'crpce.celltype', kegg, my.seed=100)
run_qusage_heatmap.seurat3(crpce, nm = 'crpce.celltype', hall.list, my.seed=100)


### monocle2
{
  # for NE and CRPC individually
  require(monocle)
  # crpce=readRDS("crpce.rds")
  temp=crpce
  # temp=subset(crpce, subset=(pheno=="CRPC"))
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
  tumor.markers=readRDS("crpce.markers.celltype.rds")
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

  pdf("crpce.plot_cell_trajectory_byState_temp.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "State", cell_size=0.5)
  print(temp)
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
  crpce@meta.data[names(temp.m.time), "pseudotime"]=temp.m.time

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
    temp.m.time=temp.m@phenoData@data$Pseudotime
    names(temp.m.time)=rownames(temp.m@phenoData@data)
    crpce@meta.data[names(temp.m.time), "pseudotime"]=temp.m.time

    temp.df=crpce@meta.data[, c("celltype", "pseudotime")]
    # get colors
    require("scales")
    uni.color=c("CellCycle", "Clike", "Luminal", "NE")
    color=hue_pal()(length(uni.color))
    # color=color[length(oorder):1]# 
    # specify Clike's colour
    # color[length(oorder)]="#FF6600"
    pdf("f.3m.crpce.CellTypeDensityAlongPseudotime.pdf", width=10, height=5)
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
      p1<- ggplot(temp.df, aes(x=pseudotime, color=celltype)) +
      geom_density() +
      labs(title="Cell type density curve",x="Pseudotime", y = "Density") +
      scale_color_manual(values=color)+
      theme_classic()
    }
    print(p1)
    dev.off()
  }


  pdf("crpc_plot_cell_trajectory_details_celltype_temp.m.pdf", width=21, height=6)
  temp=plot_cell_trajectory(temp.m, color_by = "celltype", cell_size=0.4) +
      facet_wrap(~celltype, nrow = 1)
  print(temp)
  dev.off()

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

  # saveRDS(temp.m, "monocle.crpce.rds")
  # temp.m=readRDS("monocle.crpce.rds")

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
  # saveRDS(diff_test_res.temp.m, "diff_test_res.temp.m.CRPC.rds")
  # diff_test_res.temp.m=readRDS("diff_test_res.temp.m.CRPC.rds")
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
