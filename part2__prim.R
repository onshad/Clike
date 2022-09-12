### Seurat V4
require(ggplot2)
require(RColorBrewer)
require(dplyr)
source("../functions.seurat.R")
require(Seurat)
require(GEOquery)
# [1] '1.0.0'
packageVersion("Seurat")
# [1] '4.0.1'


# working.directory:
setwd("/ ... ")

### read primary PCa samples
### process from raw counts
data.directory="/ ... /counts/primary"
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
  plot1=plot1+geom_vline(xintercept=1000)+geom_hline(yintercept=30)
  plot2 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2=plot2+geom_vline(xintercept=mean(inte.list[[i]]@meta.data$nCount_RNA))
  plot3 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2 + plot3)
  dev.off()
  print("Sample ID:")
  print(sample.names[i])
  print(dim(inte.list[[i]]))
  print(mean(inte.list[[i]]@meta.data$nCount_RNA))
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
prim=merge(x = inte.list[[1]], y = inte.list[2:length(inte.list)] )
prim@meta.data$fig.sample=as.character(prim@meta.data$orig.ident)
prim@meta.data$fig.patient=as.character(prim@meta.data$orig.ident)
prim@meta.data$pheno=rep("Prim", nrow(prim@meta.data))
prim@meta.data$pheno[prim@meta.data$fig.patient%in%c("HYQ", "PXL", "QLX-M", "RSC1023")]="mHSPC"
prim@meta.data$pheno[prim@meta.data$fig.patient%in%c("T", "CZK-T")]="ADT"
pat2gleason=data.frame(patient=sort(unique(prim@meta.data$fig.patient)))
pat2gleason$patient
 # old:
 # [1] "HYQ"        "LPX"        "PXL"        "QLX-M"      "R-M"
 # [6] "RSC1023"    "S6-TUMOR-M" "SHSCR"      "SJ-M"       "ZQF-M"
 # pat2gleason$gleason=factor(c("5", "5", "5", "5", "2", "4", "5", "5", "2", "3"))
#  [1] "CZK-T"      "DHB-T"      "HYQ"        "LPX"        "PXL"
#  [6] "QLX-M"      "R-M"        "RSC1023"    "S6-TUMOR-M" "SHSCR"
# [11] "SJ-M"       "T"          "ZQF-M"
pat2gleason$gleason=factor(c(NA, NA, "5", "5", "5", "5", "2", "4", "5", "5", "2", NA, "3"))
pat2gleason
#       patient gleason
# 1       CZK-T    <NA>
# 2       DHB-T    <NA>
# 3         HYQ       5
# 4         LPX       5
# 5         PXL       5
# 6       QLX-M       5
# 7         R-M       2
# 8     RSC1023       4
# 9  S6-TUMOR-M       5
# 10      SHSCR       5
# 11       SJ-M       2
# 12          T    <NA>
# 13      ZQF-M       3
prim@meta.data$gleason=pat2gleason$gleason[match(prim@meta.data$fig.patient, pat2gleason$patient)]
prim@meta.data$gleason=as.character(prim@meta.data$gleason)
table(prim@meta.data$gleason)
 #    2     3     4     5
 # 4310  4675  4162 29066
# check ingredients in meta.data
str(prim@meta.data)
# 'data.frame':   68339 obs. of  8 variables:
#  $ orig.ident  : chr  "CZK-T" "CZK-T" "CZK-T" "CZK-T" ...
#  $ nCount_RNA  : num  1994 3653 2746 2389 4268 ...
#  $ nFeature_RNA: int  1071 1471 1198 1227 1801 2061 1298 3889 1243 2011 ...
#  $ percent.mt  : num  4.86 9.77 21.56 16.07 8.22 ...
#  $ fig.sample  : chr  "CZK-T" "CZK-T" "CZK-T" "CZK-T" ...
#  $ fig.patient : chr  "CZK-T" "CZK-T" "CZK-T" "CZK-T" ...
#  $ pheno       : chr  "ADT" "ADT" "ADT" "ADT" ...
#  $ gleason     : chr  NA NA NA NA ...
prim=SCTransform(prim, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting

# check most variable genes
str(prim@assays$SCT@var.features)
# chr [1:3000] "TPSB2" "MSMB" "NPY" "RGS5" "TPSAB1" "SPP1" "ACPP" "PLA2G2A" ...
# top10 <- head(VariableFeatures(prim, assay="SCT"), 10)
# # plot variable features with and without labels
# toplot<- VariableFeaturePlot(prim, assay = "SCT", selection.method="sctransform")
# toplot <- LabelPoints(plot = toplot, points = top10, repel = TRUE)
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()
# str(prim)

DefaultAssay(prim)="SCT"
prim <- RunPCA(prim, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(prim,  ndims = 100)
dev.off()

### PCA-based clustering
{
  # prim <- RunUMAP(prim, reduction = "pca", dims = 1:25, verbose = FALSE)
  # prim <- FindNeighbors(prim, reduction = "pca", dims = 1:25, verbose = FALSE)
  # prim <- FindClusters(prim, verbose = FALSE, resolution = 0.2)
  # table(prim@active.ident)

  # prim.markers <- FindAllMarkers(prim, assay="SCT")
  # saveRDS(prim.markers, "prim.markers.as.fig.cluster.byPCA.rds")
  # prim.markers=readRDS("prim.markers.as.fig.cluster.byPCA.rds")
}

### Harmony-based clustering
{
  require(harmony)
  DefaultAssay(prim)
  # [1] "SCT"
  Sys.time()
  prim <- RunHarmony(prim, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  dev.off()
  Sys.time()

  # prim <- RunUMAP(prim, dims = 1:25, verbose = FALSE)
  prim <- RunUMAP(prim, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prim <- FindNeighbors(prim, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prim <- FindClusters(prim, verbose = FALSE, resolution = 0.2)
  table(prim@active.ident)
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 18145 11552 11430  7288  5183  2234  2226  2223  1745  1695  1645  1493   686
#    13    14    15    16
#   325   235   131   103
  prim.markers <- FindAllMarkers(prim, assay="SCT")
  # saveRDS(prim.markers, "prim.markers.as.fig.cluster.byHarmony.rds")
  # prim.markers =readRDS("prim.markers.as.fig.cluster.byHarmony.rds")
}


prim@meta.data$fig.cluster=prim@meta.data$seurat_clusters
prim <- SetIdent(prim, value = "fig.cluster")
pdf("prim_cluster_dimplot.pdf", width=5, height=4.5)
DimPlot(prim, label = TRUE, group.by="fig.cluster", pt.size=0.001, label.size=5, raster=TRUE)# + NoLegend()
dev.off()
png("temp_patient_dimplot.png", width=6*100, height=5*100)
DimPlot(prim, label = FALSE, group.by="fig.sample", pt.size=0.001)# + NoLegend()
dev.off()
png("temp_pheno_dimplot.png", width=6*100, height=5*100)
DimPlot(prim, label = FALSE, group.by="pheno", pt.size=0.001)# + NoLegend()
dev.off()
pdf("prim.dimplot.sample.pdf", height=4.5, width=6)
DimPlot(prim, label = FALSE, group.by="fig.sample", label.size=5, raster=TRUE)# + NoLegend()
dev.off()


### check celltype markers' expression to make sure every cluster's celltype
markers=prim.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[markers$avg_log2FC>log2(1.5), ]
setwd("..")
source("functions.seurat.R")
prim@meta.data$celltype=celltypeAnno(prim, markers)
 # [10] "BladderEpithelial_nc202010: 0.9
 # [11] "Mast cell: 0.19 T cell: 0.19"
prim@meta.data$celltype[prim@meta.data$fig.cluster=="9"]="Epithelial cell"
prim@meta.data$celltype[prim@meta.data$fig.cluster=="10"]="T cell"
prim@meta.data$celltype[prim@meta.data$fig.cluster=="11"]="B/Plasma cell"
prim@meta.data$celltype[prim@meta.data$fig.cluster=="12"]="B/Plasma cell"

prim@meta.data$celltype[prim@meta.data$fig.cluster=="6"]="Epithelial cell"
prim@meta.data$celltype[prim@meta.data$fig.cluster=="16"]="Epithelial cell"
setwd("./ ... ")
# no obscure identification
# plot with cluster names
pdf("prim_celltype_dimplot.pdf", width=6, height=4.5)
DimPlot(prim, label = TRUE, group.by="celltype", pt.size=0.001, label.size=5, raster=TRUE)# + NoLegend()
dev.off()


# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=prim.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.prim.cluster","csv",sep="."), row.name=T)


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
# check.genes=c("INHBA", "VCAN", "ZFHX4")


DefaultAssay(prim)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prim, features = check.genes)
dev.off()
DefaultAssay(prim)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(prim, features = check.genes, pt.size = 0, group.by="celltype", ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(prim, "prim.rds")
# prim=readRDS("prim.rds")


# draw important genes' dotmap
imp.genes=c()
imp.genes=c("CD3D", "CD3E", "CD3G",    "CD4", "CD14", "CD163",    "KRT8", "KRT18", "AR",     "TPSB2", "TPSAB1", "MS4A2",     "ENG", "VWF", "PECAM1",  "ACTA2", "VIM", "CNN1",    "COL1A1", "FBLN1", "FN1",       "IGHG4", "MS4A1", "MZB1",   "FTL", "CD55", "CD44")
# DES: myofibroblast
imp.genes=unique(imp.genes)
imp.genes=imp.genes[length(imp.genes):1]# reverse
# delete celltype with less than 100 cells
prim=subset(prim, subset=celltype!="Myocyte")
prim@meta.data$fig.celltype=factor(prim@meta.data$celltype, levels=c("T cell", "Macrophage", "Epithelial cell", "Mast cell", "Endothelial cell", "Myofibroblast", "Fibroblast", "B/Plasma cell", "Neutrophil"))
pdf("prim_dotmap_celltype_markers.pdf", width = 18, height = 3);
DotPlot(prim, features = imp.genes, cols = c('white', 'blue'), dot.scale = 6, group.by="fig.celltype");
dev.off();


### further analysis focused on epithelial cells
table(prim@meta.data$celltype)
   # B/Plasma cell       Basal cell Endothelial cell  Epithelial cell
   #         14092            16515            17314            49178
   #    Fibroblast     Luminal cell       Macrophage        Mast cell
   #          4508             3933            25145             4160
   # Myofibroblast       Neutrophil           T cell
   #         10326             8175            43894
epitypes=c("Epithelial cell", "Luminal cell", "Basal cell")
prime=subset(prim, subset=(celltype%in%epitypes))
str(prime)
prime@commands=list()


dim(prime)
# [1] 22146 15576
DefaultAssay(prime)="SCT"
str(prime@assays$SCT@var.features)
# chr [1:3000] "TPSB2" "MSMB" "NPY" "RGS5" "TPSAB1" "SPP1" "ACPP" "PLA2G2A" ...
# prime <- FindVariableFeatures(prime, assay="SCT", nfeatures = 3000)
# str(prime@assays$SCT@var.features)
# top10 <- head(VariableFeatures(prime, assay="SCT"), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(prime, assay = "SCT")
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# toplot=plot1 + plot2
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()

prime=SCTransform(prime, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
prime <- RunPCA(prime, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(prime,  ndims = 100)
dev.off()

### clustering basd on PCA
{
  # DefaultAssay(prime)="SCT"
  # prime <- RunUMAP(prime, dims = 1:20, verbose = FALSE)
  # prime <- FindNeighbors(prime, dims = 1:20, verbose = FALSE)

  # DefaultAssay(prime)="SCT" 
  # prime <- FindClusters(prime, verbose = FALSE, resolution = 0.2)
  # table(prime@meta.data$seurat_clusters)
  # prime.markers <- FindAllMarkers(prime, assay="SCT")
  # saveRDS(prime.markers, "prime.markers.as.fig.cluster.byPCA.rds")
  # prime.markers =readRDS("prime.markers.as.fig.cluster.byPCA.rds")

}

### clustering basd on Harmony
{
  require(harmony)
  DefaultAssay(prime)
  # [1] "SCT"
  Sys.time()
  prime <- RunHarmony(prime, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
  dev.off()
  Sys.time()

  # prim <- RunUMAP(prim, dims = 1:25, verbose = FALSE)
  prime <- RunUMAP(prime, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prime <- FindNeighbors(prime, reduction = "harmony", dims = 1:25, verbose = FALSE)
  prime <- FindClusters(prime, verbose = FALSE, resolution = 0.2)
  table(prime@meta.data$seurat_clusters)
  #    0    1    2    3    4    5    6    7    8    9   10
  # 4069 3853 2325 1659 1446 1289  499  229   88   76   43
  prime.markers <- FindAllMarkers(prime, assay="SCT")
  # saveRDS(prime.markers, "prime.markers.as.fig.cluster.byHarmony.rds")
  # prime.markers =readRDS("prime.markers.as.fig.cluster.byHarmony.rds")
}

markers=prime.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.prime.cluster","csv",sep="."), row.name=T)

prime@meta.data$pheno[prime@meta.data$fig.patient%in%c("DHB-T", "HYQ", "PXL", "QLX-M", "RSC1023")]="mHSPC"# DHB-T wasn't annotated as HSPC before, but we later knew this tumor had stage M1 
prime@meta.data$fig.cluster=as.character(prime@meta.data$seurat_clusters)
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(prime, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
pdf("prime.dimplot.patient.pdf", height=4.5, width=5.5)
DimPlot(prime, label = FALSE, group.by="fig.patient", label.size=5, raster=TRUE)
dev.off()
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(prime, label = FALSE, group.by="gleason")# + NoLegend()
dev.off()
png("dimplot.pheno.png", height=5*100, width=5*100)
DimPlot(prime, label = FALSE, group.by="pheno")# + NoLegend()
dev.off()
prime@meta.data$temp="Other cells"
prime@meta.data$temp[prime@meta.data$pheno=="ADT"]="Cells from post-ADT"
png("dimplot.png", height=5*100, width=6*100)
temp=DimPlot(prime, label = FALSE, group.by="temp", label.size=8, cells.highlight=colnames(prime)[prime@meta.data$pheno=="ADT"], sizes.highlight=0.5) # + NoLegend()
print(temp)
dev.off()


# cluster 6 is T cell(CD3E+), remove it
# cluster 7 is B cell(MS4A1+, CD79A+), remove it
# remove clusters 8, 9 and 10 with too few cells (<100 cells)
prime=subset(prime, cells = colnames(prime)[!prime@meta.data$fig.cluster%in%c("6", "7", "8", "9", "10")])
table(prime@meta.data$fig.cluster)
#    0    1    2    3    4    5
# 4069 3853 2325 1659 1446 1289

# plot previously identified basal/cycle signatures
# all 203 basal markers
# prime <- ScaleData(prime, features = rownames(prime), assay="RNA" )
basal.posi.genes=read.table("../basal.in.all.cell.markers.posi.genes.txt", head=FALSE, stringsAsFactors=FALSE)$V1
# all 126 cycle markers
cycle.posi.genes=read.table("../cycle.in.all.cell.markers.posi.genes.txt", head=FALSE, stringsAsFactors=FALSE)$V1
basal.posi.genes=basal.posi.genes[basal.posi.genes%in%rownames(prime@assays$SCT@scale.data)]
cycle.posi.genes=cycle.posi.genes[cycle.posi.genes%in%rownames(prime@assays$SCT@scale.data)]
prime@meta.data$basal.sig=colMeans(prime@assays$SCT@scale.data[basal.posi.genes, ])
DefaultAssay(prime)="SCT"
pdf("temp.basal.sig.featureplot.pdf", width=5, height=5)
FeaturePlot(prime, features = "basal.sig")
dev.off()
prime@meta.data$cycle.sig=colMeans(prime@assays$SCT@scale.data[cycle.posi.genes, ])
pdf("temp.cycle.sig.featureplot.pdf", width=5, height=5)
FeaturePlot(prime, features = "cycle.sig")
dev.off()
DefaultAssay(prime)="SCT"
pdf("temp.VlnPlot.pdf", width = 11, height =4)
tocheck=c("basal.sig", "cycle.sig")
temp=VlnPlot(prime, features = tocheck, ncol = 2,  pt.size = 0, group.by="fig.cluster", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
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
# check.genes=c("SIX3", "MIAT", "SG3", "MIAT", "CHD7", "KIF19", "NEB", "RGS16", "PCP4", "GPX2", "ID4", "CALB1", "SGB2A1")

prime <- SetIdent(prime, value = "fig.cluster")
DefaultAssay(prime)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prime, features = check.genes)
dev.off()
DefaultAssay(prime)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(prime, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# indi.plot(prime, tocheck="fig.patient")
# indi.plot(prime, tocheck="fig.cluster")

### re-cluster basal & typeC cells
{
  primb=subset(prime, subset=(fig.cluster%in%"2"))
  dim(primb)
  # [1] 20066  2325
  DefaultAssay(primb)="SCT"
  str(primb@assays$SCT@var.features)
  # chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "MMP7" "KRT13" "LCN2" ...
  # primb <- FindVariableFeatures(primb, assay="SCT", nfeatures = 3000)
  # str(primb@assays$SCT@var.features)
  # top10 <- head(VariableFeatures(primb, assay="SCT"), 10)
  # # plot variable features with and without labels
  # plot1 <- VariableFeaturePlot(primb, assay = "SCT")
  # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # toplot=plot1 + plot2
  # png("temp.png", width=10*100, height=5*100)
  # print(toplot)
  # dev.off()

  primb=SCTransform(primb, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  primb <- RunPCA(primb, assay="SCT", verbose = FALSE, npcs = 100)
  ElbowPlot(primb,  ndims = 100)
  dev.off()

  ### clustering basd on PCA
  {
    # DefaultAssay(primb)="SCT"
    # primb <- RunUMAP(primb, dims = 1:20, verbose = FALSE)
    # primb <- FindNeighbors(primb, dims = 1:20, verbose = FALSE)

    # DefaultAssay(primb)="SCT" 
    # primb <- FindClusters(primb, verbose = FALSE, resolution = 0.2)
    # table(primb@meta.data$seurat_clusters)
    # primb.markers <- FindAllMarkers(primb, assay="SCT")
    # saveRDS(primb.markers, "primb.markers.as.fig.cluster.byPCA.rds")
    # primb.markers =readRDS("primb.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(primb)
    # [1] "SCT"
    Sys.time()
    primb <- RunHarmony(primb, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    dev.off()
    Sys.time()

    # prim <- RunUMAP(prim, dims = 1:25, verbose = FALSE)
    primb <- RunUMAP(primb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    primb <- FindNeighbors(primb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    primb <- FindClusters(primb, verbose = FALSE, resolution = 0.2)
    table(primb@meta.data$seurat_clusters)
    #   0   1   2   3   4
    # 894 652 577 171  31
    # primb.markers <- FindAllMarkers(primb, assay="SCT")
    # saveRDS(primb.markers, "primb.markers.as.fig.cluster.byHarmony.rds")
    # primb.markers =readRDS("primb.markers.as.fig.cluster.byHarmony.rds")
  }

  primb@meta.data$fig.cluster=as.character(primb@meta.data$seurat_clusters)
  png("dimplot.cluster.png", height=5*100, width=5*100)
  DimPlot(primb, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  dev.off()
  pdf("dimplot.patient.pdf", height=4.5, width=5)
  DimPlot(primb, label = TRUE, group.by="fig.patient", label.size=5, raster=TRUE)
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  DimPlot(primb, label = TRUE, group.by="gleason")
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  DimPlot(primb, label = TRUE, group.by="pheno")
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

  primb <- SetIdent(primb, value = "fig.cluster")
  DefaultAssay(primb)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(primb, features = check.genes)
  dev.off()
  DefaultAssay(primb)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=8*100)
  temp=VlnPlot(primb, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
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
  primb@meta.data$epitype <- plyr::mapvalues(x = as.character(primb@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  png("dimplot.epitype.pdf", height=4.5, width=5)
  DimPlot(primb, label = TRUE, group.by="epitype", label.size=8, raster=TRUE)
  dev.off()

  btype=primb@meta.data$epitype
  names(btype)=rownames(primb@meta.data)
  table(btype)
  # btype
  # tB21 tC20 tC23 tC24 tL22
  #  652  894  171   31  577

  # indi.plot(primb, tocheck="fig.patient")
  # indi.plot(primb, tocheck="fig.cluster")

  # saveRDS(primb, "primb.rds")
  # primb=readRDS("primb.rds")
}

### epitype annotation
current.cluster.ids <- c(0:5, 7)
new.cluster.ids <- c("tL0", "tL1", "tB2", "tL3", "tL4", "tL5", "tCy7")
# tL: tumor luminal, tB: tumor basal, tCy: tumor CellCycle
prime@meta.data$epitype <- plyr::mapvalues(x = as.character(prime@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot
prime@meta.data[names(btype), "epitype"]=btype
pdf("prime.dimplot.epitype.pdf", height=4.5, width=5)
DimPlot(prime, label = TRUE, group.by="epitype", label.size=5, raster=TRUE)
dev.off()

### celltype annotation
current.cluster.ids <- unique(prime@meta.data$epitype)
current.cluster.ids=sort(current.cluster.ids)
current.cluster.ids
#  [1] "tB21" "tC20" "tC23" "tC24" "tCy7" "tL0"  "tL1"  "tL22" "tL3"  "tL4"
# [11] "tL5"
new.cluster.ids <- c("Basal", "TypeC", "Clike", "TypeC", "CellCycle", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal")
# tLum: tumor luminal, tBI: tumor basal/intermediate, tCyc: tumor CellCycle
prime@meta.data$celltype <- plyr::mapvalues(x = as.character(prime@meta.data$epitype), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot

### f2.e
png("f2.e.prim.celltype.dimplot.png", height=5*100, width=5*100)
DimPlot(prime, label = TRUE, group.by="celltype", label.size=5)
dev.off()
pdf("f2.e.prim.celltype.dimplot.pdf", height=4.5, width=5)
temp.color=color[names(color)%in%unique(prime@meta.data$celltype)]
temp.oorder=oorder[oorder%in%unique(prime@meta.data$celltype)]
DimPlot(prime, label = TRUE, group.by="celltype", label.size=5, raster=TRUE, cols=temp.color, order=temp.oorder)
dev.off()

### f2f
imp.genes=c(  "AR", "KRT8",    "KRT5",    "CHGA","SYP","ENO2",    "CYP1A1", "ARL14",    "KRT4", "PSCA")
imp.genes=c("IL1A", "ARL14", "FABP4", "FGFR3", "TGFBR3", "TACSTD2", "GPRC5A", "UPK1A", "CD55" )
imp.genes=c( "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "FGFR3", "TGFBR3")
# "CD55", "FGFR3",
imp.genes=c( "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "FGFR3", "TGFBR3")

### f2b
imp.genes=c(  "AR", "KRT8",    "KRT5",    "CHGA","SYP","ENO2",    "CYP1A1", "ARL14",    "KRT4", "PSCA")
# "CD55", "FGFR3",
imp.genes=c(  "KRT4", "PSCA",     "KRT5", "AR", "KRT8",   "CHGA","SYP","ENO2")
imp.genes=c(  "KRT4", "TACSTD2", "PSCA",     "KRT5", "AR", "KRT8",   "CHGA","SYP","ENO2",   "NCAM1", "ACSL4")

DefaultAssay(prime)="SCT"
temp.color=color[names(color)%in%unique(prime@meta.data$celltype)]
temp.oorder=oorder[oorder%in%unique(prime@meta.data$celltype)]
prime@meta.data$fig.celltype=factor(prime@meta.data$celltype, levels=temp.oorder)
prime <- SetIdent(prime, value = "fig.celltype")
pdf("prime.VlnPlot.pdf", width=20, height=2)
temp=VlnPlot(prime, features = imp.genes, pt.size = 0, group.by=NULL, ncol = 11, cols=temp.color)
for(i in 1:length(imp.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 5, colour = "black", shape = 95)
}
print(temp)
dev.off()



{
  tocheck=c( "IL1A", "ARL14", "FGFR3",  "PSCA", "KRT4" ,  "KRT5", "KRT19",  "KLK3", "AR" )
  # tocheck=c("TACSTD2", "GPRC5A", "UPK1A", "CD55", "FGFR3", "TGFBR3", "IL1A", "ARL14", "FABP4")
  # tocheck=c( "KRT5", "KRT19", "PSCA", "KRT4" ,   "KLK3", "ACPP",    "AMACR", "AR",     "TOP2A", "STMN1",      "ENO2", "FOXA1",  "cnv.counts", "nFeature_RNA") # MMP7, AMACR
  # pdf("f1_epithelial.dotmap_celltype_markers.pdf", width = 12, height = 2);
  # DotPlot(temp, features = tocheck, cols = c('white', 'blue'), dot.scale = 6, group.by="fig.celltype");
  # dev.off();
  # tocheck=c("KRT4", "PSCA", "KRT5", "IL1A", "ARL14", "UPK1A", "CD55",  "FGFR3"  )

  pdf("VlnPlot_prime_celltype.pdf", width = 10, height =20)
  # keep colors consistent with all-epithelial-celltypes image
  oorder=c("Clike", "TypeC",  "Basal", "Luminal")
  # oorder=oorder[length(oorder):1]
  # get colour
  require("scales")
  color=hue_pal()(length(oorder))
  # color=color[length(oorder):1]# 
  # specify Clike's colour
  # color[length(oorder)]="#FF6600"
  # oorder=oorder[length(oorder):1]

  prime@meta.data$factor.fig.celltype=factor(prime@meta.data$celltype, levels=oorder)
  temp=VlnPlot(prime, features = tocheck, pt.size = 0, group.by="factor.fig.celltype", ncol = 1, cols = NULL)
  for(i in 1:length(tocheck)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
  }
  # adjust the border line's width, not succesfully
  temp$layers[[1]]$aes_params$size = 0.1
  # adjust = 0.1 # refer to ggplot2.tidyverse.org/reference/geom_violin.html
  print(temp)
  dev.off()
  DefaultAssay(prime)="SCT"
}


check.genes=c("SIX3", "MIAT", "SG3", "MIAT", "CHD7", "KIF19", "NEB", "RGS16", "PCP4", "GPX2", "ID4", "CALB1", "SGB2A1")
check.genes=c("ASCL1", "NEUROD1")
check.genes=c("PSCA", "LY6D", "ATXN1")
check.genes=c( "LY6D", "CD24A", "KRT8", "KRT18", "PROM1", "EPCAM", "KRT7", "KRT19", "PSCA", "ELF3")
prime <- SetIdent(prime, value = "celltype")
DefaultAssay(prime)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prime, features = check.genes)
dev.off()
DefaultAssay(prime)="SCT"
png("temp.VlnPlot.png", width=10*100, height=12*100)
temp=VlnPlot(prime, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(prime, "prime.rds")
# prime=readRDS("prime.rds")


### Not employed --->>>
  ### CCA to identify epithelial cell types WithiN luminal tumor cells
  # norme=readRDS("norme.rds")
  priml=subset(prime, subset=celltype=="Luminal")
  anchors <- FindTransferAnchors(reference = norme, query = priml, 
      dims = 1:30, normalization.method="LogNormalize", reference.assay="SCT")
  predictions <- TransferData(anchorset = anchors, refdata = norme$celltype, 
      dims = 1:30)
  priml@meta.data$predtype=predictions$predicted.id 
  png("dimplot.predtype.png", height=5*100, width=5*100)
  DimPlot(priml, label = TRUE, group.by="predtype")
  dev.off()
### Not employed <<<---


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(prime@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(prime@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(prime@meta.data$fig.patient, prime@meta.data$celltype)
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
write.csv(pct.matrix, file="part4__prime__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


table(prime@meta.data$pheno)
  # ADT mHSPC  Prim
  # 327  6316  7998
prime@meta.data$temp <- prime@meta.data$pheno
mHSPC=c("DHB-T", "HYQ", "PXL", "QLX-M", "RSC1023")# DHB-T wasn't annotated as HSPC before, but we later knew this tumor had stage M1 
prime@meta.data$temp[prime@meta.data$fig.patient%in%mHSPC]="mHSPC"
prime@meta.data$temp[prime@meta.data$temp=="Prim"]="Primary PCa"
prime@meta.data$temp[prime@meta.data$temp=="Norm"]="Normal prostate"
prime@meta.data$temp[prime@meta.data$temp=="ADT"]="post-ADT"
prime@meta.data$temp=factor(prime@meta.data$temp, levels=c("Primary PCa", "mHSPC", "post-ADT"))
prime@meta.data$fig.pheno=prime@meta.data$temp
table(prime@meta.data$fig.pheno)
# Primary PCa       mHSPC    post-ADT
#        6970        7344         327

### barplot Clike percentage for CRPC+NEPC vs mHSPC vs Primary and incidental aggressive vs indolent
{
  # give sample ID as phenotype + number
  {
    pheno.sample=prime@meta.data[, c("fig.pheno", "fig.sample")]
    pheno.sample$fig.sample=as.character(pheno.sample$fig.sample)
    # str( pheno.sample)
    # 'data.frame':   14641 obs. of  4 variables:
    #  $ fig.pheno    : Factor w/ 3 levels "Primary PCa",..: 1 1 1 1 1 1 1 1 1 1 ...
    #  $ fig.sample   : Factor w/ 13 levels "LPX","R-M","S6-TUMOR-M",..: 1 1 1 1 1 1 1 1 1 1 ...
    #  $ fig.sample.id: num  1 1 1 1 1 1 1 1 1 1 ...
    #  $ pheno.sample : Factor w/ 13 levels "Primary PCa 1",..: 1 1 1 1 1 1 1 1 1 1 ...
    pheno.sample=pheno.sample[order(pheno.sample$fig.pheno, pheno.sample$fig.sample, decreasing=c(FALSE, FALSE)), ]
    unique.sample=unique(pheno.sample$fig.sample)
    pheno.sample$fig.sample=factor(pheno.sample$fig.sample, levels=unique.sample)
    pheno.sample$fig.sample.id=as.numeric(pheno.sample$fig.sample)
    pheno.sample$pheno.sample=paste(as.character(pheno.sample$fig.pheno), pheno.sample$fig.sample.id)
    unique.pheno.sample=unique(pheno.sample$pheno.sample)
    pheno.sample$pheno.sample=factor(pheno.sample$pheno.sample, levels=unique.pheno.sample)
    prime@meta.data$pheno.sample=pheno.sample[rownames(prime@meta.data), ]$pheno.sample
  }

  # get percentage table of phenotype+sampleID vs celltype
  celltype.sta=table(prime@meta.data$pheno.sample, prime@meta.data$celltype)
  celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
  print(unique(rownames(celltype.sta)))
  celltype.sta=celltype.sta[c("Clike", "TypeC", "Basal" ), ]# 
  barN=ncol(celltype.sta)
  n_top = nrow(celltype.sta) 
  pdf(paste("prime.barplot.celltypeInSample.pdf",sep=""))
  # the width between paper edge and plot ( down,left,up,right : four directions)
  par(mar = c(14, 4, 1, 8), xpd=T)
  # coll=DiscretePalette(n_top, palette = "polychrome")
  require("scales")
  # get 4 colors, then delete the luminal's color 
  coll=hue_pal()(n_top+1)
  # coll=coll[c(-3)]
  barplot( celltype.sta, xlab="", ylab="Percentage", col=coll, names.arg=colnames(celltype.sta), las=2, cex.names=1) # , xaxt="n"
  # facts=barN+1
  end_point = barN #0.5 + n_top *facts-1
  # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
  legend( barN+3, 15, rownames(celltype.sta), fill=coll );
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
  dev.off()

}


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(prime@meta.data$epitype)
# tB21 tC20 tC23 tC24  tL0  tL1 tL22  tL3  tL4  tL5
 # 652  894  171   31 4069 3853  577 1659 1446 1289
 table(prime@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(prime@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(prime@meta.data$fig.patient, prime@meta.data$epitype)
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
write.csv(pct.matrix, file="part4__prime__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


### check correlation between NE cell number and Clike cell number
### get cell type percentage in each sample, 
### identify NE marker+ cells as NE cells
prime@meta.data$celltype.ne=prime@meta.data$celltype
tocheck=c("CHGA", "SYP", "ENO2")
temp=colMeans(prime@assays$SCT@counts[tocheck, ])>0
prime@meta.data$celltype.ne[temp]="NE"
celltype.sta=table(prime@meta.data$fig.sample, prime@meta.data$celltype.ne)
celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
cor(celltype.sta["Clike", ], celltype.sta["NE", ])
# [1] -0.244116


# saveRDS(prime, "prime.rds")
# prime=readRDS("prime.rds")


### celltype-specific DEG identification
table(prime@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   12893     925
DefaultAssay(prime)="SCT"
### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(prime@meta.data$celltype)
celltype=sort(celltype)
prime = SetIdent(prime, value = "celltype")
for(i in celltype)
{
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        markerlist[[temp.name]]=FindMarkers(prime, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
      }
  }
}
# saveRDS(markerlist, "prime.celltype.pairwise.rds")
# markerlist=readRDS("prime.celltype.pairwise.rds")

### cannonical process to find specific markers of Clike cells
{
  markerlist=readRDS("prime.celltype.pairwise.rds")
  sp.markerlist=list()
  celltype=unique(prime@meta.data$celltype)
  celltype=sort(celltype)
  for(i in celltype)
  {
    temp.genes=rownames(prime@assays$SCT@data)
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          cur.list=markerlist[[temp.name]]
          cur.list=cur.list[order(cur.list$avg_log2FC, decreasing=TRUE), ]
          # cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(1.5),]) # best intersect between prime, crpce, triple
          cur.genes=rownames( cur.list[cur.list$avg_log2FC>log2(1.5) & cur.list$pct.1>0.4,] )
          # cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(2) & cur.list$pct.1>0.3 & cur.list$pct.2<0.1,]) # canonical markers used in current manuscript
          temp.genes=intersect(temp.genes, cur.genes)
        }
    }
    sp.markerlist[[i]]=temp.genes
  }
  str(sp.markerlist)
  sp.markerlist[["Luminal"]]=c(sp.markerlist[["Luminal"]], "KLK3", "AR")
  sp.markerlist[["TypeC"]]=c(sp.markerlist[["TypeC"]], "KRT4", "TACSTD2", "PSCA")

  temp=subset(prime, downsample=300)
  sp.markerlist=sp.markerlist[c("Clike", "TypeC", "Basal", "Luminal")]
  tocheck=Reduce(union, sp.markerlist)
  temp <- SetIdent(temp, value = "celltype")
  levels(temp)=c("Clike", "TypeC", "Basal", "Luminal")
  pdf(paste("heatmap.prime.typeC.pdf", sep=""),  width=10, height=length(tocheck)/3)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

  prime.clike.markers=sp.markerlist$Clike
  str(prime.clike.markers)
}

temp=Reduce(intersect, list( convert_mouse_to_human_symbols(triple.clike.markers), prime.clike.markers, crpce.clike.markers, incidente.clike.markers))
temp=Reduce(intersect, list( convert_mouse_to_human_symbols(triple.clike.markers), prime.clike.markers, crpce.clike.markers ))
temp=Reduce(intersect, list( prime.clike.markers, crpce.clike.markers ))
best=c("CD55", "IL1RN", "IL1A", "UPK1A", "UPK1B", "PIM1", "PMAIP1")


## trial process to find specific/surface markers of Type C tumor cells
markerlist=readRDS("prime.celltype.pairwise.rds")
sp.markerlist=list()
celltype=unique(prime@meta.data$celltype)
celltype=sort(celltype)
for(i in celltype)
{
  temp.genes=rownames(prime@assays$SCT@data)
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        cur.list=markerlist[[temp.name]]
        cur.list=cur.list[order(cur.list$avg_log2FC, decreasing=TRUE), ]
        # cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(1.5) & cur.list$pct.1>0.2 , ])  # surface marker parameters
        # cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(2) & cur.list$pct.1>0.5 , ])  # specific marker parameters
        cur.genes=rownames(cur.list[cur.list$avg_log2FC>log2(1.5) & cur.list$pct.1>0.1 , ])  # function enrichment's marker parameters
        temp.genes=intersect(temp.genes, cur.genes)
      }
  }
  sp.markerlist[[i]]=temp.genes
}
str(sp.markerlist)
# sp.markerlist[["Luminal"]]=c(sp.markerlist[["Luminal"]], "KLK3", "AR")
# sp.markerlist[["TypeC"]]=c(sp.markerlist[["TypeC"]], "KRT4", "TACSTD2", "PSCA")

temp=subset(prime, downsample=300)
sp.markerlist=sp.markerlist[c("Clike", "TypeC", "Basal", "Luminal")]
tocheck=Reduce(union, sp.markerlist)

# require(msigdbr)
# cc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
# cc = cc %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
# cdgenes=sort(unique(cc$gene_symbol[cc$gs_name=="GO_CELL_SURFACE"]))
# tocheck=tocheck[tocheck%in%cdgenes]

temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("Clike", "TypeC", "Basal", "Luminal")
pdf(paste("heatmap.prime.typeC.pdf", sep=""),  width=10, height=length(tocheck)/3)
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

enrichment.gene.set.plot(sp.markerlist[["Clike"]], topNterm=20)
enrichment.gene.set.plot(sp.markerlist[["TypeC"]], topNterm=20)

temp.genes=markerlist[["Clike__TypeC"]]
temp.genes=temp.genes[temp.genes$avg_log2FC>log2(1.5) & temp.genes$pct.1>0.1 , ]
str(temp.genes)
enrichment.gene.set.plot(rownames(temp.genes), topNterm=20)
### heatmap intersection of DEGs and genes of Pathway
{
  require(clusterProfiler)
  require(msigdbr)
  m_df.hallmarker = msigdbr(species = "Homo sapiens", category = "H", subcategory = "")
  m_t2g.hallmarker = m_df.hallmarker %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  str(m_t2g.hallmarker)
  # pathway.name="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
  # pathway.name="HALLMARK_HYPOXIA"
  pathway.name="HALLMARK_ESTROGEN_RESPONSE_EARLY"
  pathway.name2="HALLMARK_ESTROGEN_RESPONSE_LATE"
  tocheck=m_t2g.hallmarker$gene_symbol[m_t2g.hallmarker$gs_name==pathway.name]
  tocheck2=m_t2g.hallmarker$gene_symbol[m_t2g.hallmarker$gs_name==pathway.name2]
  tocheck=unique(c(tocheck, tocheck2))
  tocheck=intersect(tocheck, rownames(temp.genes))

  prime<- SetIdent(prime, value = "celltype")
  temp=subset(prime, downsample=300)
  levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
  pdf(paste("temp.heatmap.pdf", sep=""),  width=10, height=length(tocheck)/3)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
}


markerlist=readRDS("prime.celltype.pairwise.rds")

temp=markerlist[["Clike__TypeC"]]
temp=temp[temp$avg_log2FC>log2(1.5) & temp$pct.1>0.1 & temp$p_val_adj<0.05, ]
str(temp)
enrichment.gene.set.plot(rownames(temp), topNterm=20)

temp=markerlist[["Clike__Luminal"]]
temp=temp[temp$avg_log2FC>log2(1.5)  & temp$pct.1>0.1 & temp$p_val_adj<0.05, ] #  
str(temp)
enrichment.gene.set.plot(rownames(temp), topNterm=15)


### manually pick up some markers witih good specificity
tocheck=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "CLCN3", "SEMA7A", "DSG2", "VTCN1", "LTF", "CD74", "CFTR", "KRT4", "C3", "INSR", "CEACAM6", "CP", "TM4SF1", "RARRES1", "LCN2", "PPP1R1B", "EGFR", "SFRP1", "NGFR", "PRNP", "FOLH1", "ACPP")
tocheck=c("TP63", "TRIM29", "ITGB4", "KRT5", "KRT14", "AR", "KLK3", "DPP4", "ESR1", "ESR2", "NR3C1", "PGR", "DUSP6", "SPRY4", "ETV4", "GDF15", "EGR1", "ETV5", "C3orf52", "UPP1", "DUSP5", "TRIB3", "MYEOV", "MAFF", "PHLDA1")
tocheck=c("EPCAM", "KRT8", "KRT5", "TP63", "TACSTD2", "CD55", "FGFR3", "UPK1A", "TGFBR3")
tocheck=c( "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "FGFR3", "TGFBR3")

prime<- SetIdent(prime, value = "celltype")
temp=subset(prime, downsample=300)
levels(temp)=c("Clike", "TypeC", "Basal", "Luminal")
pdf(paste("heatmap.prime.Clike.pdf", sep=""),  width=5, height=5)#height=length(tocheck)/3
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)#+ scale_fill_gradientn(colors = c("green", "red"))
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()

### check whether KRT8+ FGFR+ basal cells are actually Clike cells
{
  tocheck=c("EPCAM", "KRT8", "KRT5", "TP63", "TACSTD2", "CD55", "FGFR3", "UPK1A", "TGFBR3")
  png("temp.VlnPlot.png.png", width=10*100, height=8*100)
  temp=VlnPlot(prime, features = tocheck, pt.size = 0, group.by="celltype", ncol = 2)
  for(i in 1:length(tocheck)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  # find basal cell with FGFR3 and KRT8 expression > median value of Clike cells
  prime@meta.data$basal.plus=prime@meta.data$celltype
  fgfr3median=median(prime@assays$SCT@data["FGFR3", prime@meta.data$celltype=="Clike"])
  krt8median=median(prime@assays$SCT@data["KRT8", prime@meta.data$celltype=="Clike"]) 
  prime@meta.data$basal.plus[prime@meta.data$celltype=="Basal" & prime@assays$SCT@data["FGFR3", ]>fgfr3median & prime@assays$SCT@data["KRT8", ]>krt8median ]="Basal.FGFR3+KRT8+"
  table(prime@meta.data$basal.plus)
            # Basal Basal.FGFR3+KRT8+             Clike           Luminal
            #   636                16               171             12893
            # TypeC
            #   925
  tocheck=c("EPCAM", "KRT8", "KRT5", "TP63", "TACSTD2", "CD55", "FGFR3", "UPK1A", "TGFBR3")
  png("temp.VlnPlot.png.png", width=10*100, height=12*100)
  temp=VlnPlot(prime, features = tocheck, pt.size = 0, group.by="basal.plus", ncol = 2)
  for(i in 1:length(tocheck)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

}

### check specific marker across all cell types
# prim=readRDS("prim.rds")
check.genes=c("IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1")
clike.cells=rownames(prime@meta.data)[prime@meta.data$celltype=="Clike"]
prim@meta.data$celltype[rownames(prim@meta.data)%in%clike.cells]="Clike"
prim <- SetIdent(prim, value = "celltype")
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prim, features = check.genes)
dev.off()
DefaultAssay(prim)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(prim, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

prim <- SetIdent(prim, value = "celltype")
clike.allcelltype=FindMarkers(prim, ident.1 = "Clike", assay="SCT")
# saveRDS(clike.allcelltype, "markers.clikeVSallcelltypes.rds")
# clike.allcelltype=readRDS("markers.clikeVSallcelltypes.rds")
markers=clike.allcelltype
markers=markers[order(markers$avg_log2FC,decreasing = TRUE), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("foldchange","pct.1", "pct.2", "p_val_adj")], file=paste("markers.clike.in.prim","csv",sep="."), row.name=T)



### DEG identification, across all clusters
prime <- SetIdent(prime, value = "celltype")
prime.markers <- FindAllMarkers(prime, slot="data", assay="SCT")
# saveRDS(prime.markers, "prime.markers.celltype.rds")
# prime.markers=readRDS("prime.markers.celltype.rds")


prime <- SetIdent(prime, value = "celltype")
prime.markers <- FindAllMarkers(prime, slot="data", assay="SCT")
# saveRDS(prime.markers, "prime.markers.celltype.rds")
# prime.markers=readRDS("prime.markers.celltype.rds")


### plot heatmap
# prime.markers=readRDS("prime.indivSCT.markers.as.fig.cluster.rds")
# prime <- NormalizeData(prime, normalization.method = "LogNormalize", scale.factor = 10000, assay="SCT" )
# prime <- ScaleData(prime, features = rownames(prime), assay="RNA")

# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=prime.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.prime.celltype","csv",sep="."), row.name=T)

markers=prime.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.prime.celltype.NOfilter","csv",sep="."), row.name=T)



# heatmap top 10 DEGs for every cluster
markers=prime.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(2) & markers$pct.1>0.3 & markers$pct.2<0.1, ]
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
top10$cluster=factor(top10$cluster, levels=c("Clike", "TypeC", "Luminal", "Basal", "CellCycle"))
top10=top10[order(top10$cluster), ]
pdf(paste("heatmap.prime.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=nrow(top10)/8)
temp=subset(prime, downsample=300)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
tocheck=c(top10$gene, "KLK3", "AR")
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


# check cell-surface markers of each celltype
markers=prime.markers
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
  DefaultAssay(prime)="SCT"
  pdf(paste("temp.heatmap.CD_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
  temp=subset(prime, downsample=1000)
  levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
  temp.heatmap<-DoHeatmap(temp , features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
  DefaultAssay(prime)="SCT"
}


### pairwise DEGs to find three surface marker for Clike cells
markerlist=readRDS("prime.celltype.pairwise.rds")
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
temp.marker=c("UPK1A", "CD55", "TGFBR3", "FGFR3")
temp=subset(prime, downsample=1000)
# temp=NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA" )
# temp=ScaleData(temp, features = rownames(temp), assay="RNA" )
levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
pdf(paste("temp.heatmap.pdf", sep=""), width=(length(unique(prime@active.ident))/4+5), height=length(temp.marker)/10+3)
temp.heatmap<-DoHeatmap(temp , features = temp.marker, slot="scale.data", assay="SCT", angle = 0, size=3) 
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### plot heatmap for special categories of DEGs
setwd("..")
subtype_markers=subtypeMarker(prime)
setwd("./ ... ")
prime.markers=readRDS("prime.markers.celltype.rds")
markers=prime.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
# for TFs, change log2(1.5) to log(1.2)
# prime <- ScaleData(prime, features = rownames(prime), assay="RNA" )
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
    DefaultAssay(prime)="SCT"
    pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
    temp=subset(prime, downsample=300)
    temp <- SetIdent(temp, value = "celltype")
    levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
    temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()
    DefaultAssay(prime)="SCT"
  }
}


# for too few cluster DExpressing genes, do vlnplot
MHCII=c("HLA-DRB3", "HLA-DRA", "HLA-DRB1", "HLA-DPB1", "HLA-DPA1", "HLA-DQA1", "HLA-DQB1")
prime@meta.data$MHCII=colMeans(prime@assays$SCT@data[MHCII[MHCII%in%rownames(prime@assays$SCT@data)], ])
tocheck=c("CCL4", "CCL20", "CXCL2", "CXCL8", "CXCL16", "CXCL17", "CX3CL1",  "CEACAM1", "LGALS9", "NECTIN2", "MHCII")
pdf("prime.vlnPlot.chemokines_celltype.pdf", width = 10, height =20)
temp.prime=subset(prime, celltype%in%c("Clike", "TypeC", "Basal", "Luminal"))
temp.prime@meta.data$celltype=factor(temp.prime@meta.data$celltype, levels=c( "Clike", "TypeC", "Basal", "Luminal" ))
temp=VlnPlot(temp.prime, features = tocheck, pt.size = 0, group.by="celltype", ncol = 1)
# for(i in 1:length(tocheck)) 
# {
#   temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
# }
print(temp)
dev.off()


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
tumor.markers=readRDS("prime.markers.celltype.rds")
markers=tumor.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[abs(markers$avg_log2FC)>log2(1.5), ]
markers$foldChange=2^(markers$avg_log2FC)
# markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
markers=markers[order(markers$foldChange, decreasing=TRUE), ]
ordering.genes=unique(markers$gene)

### check average value of genesets
jakstat=jakstat.bak
jakstat=jakstat[jakstat%in%rownames(prime@assays$RNA@scale.data)]
prime@meta.data$jakstat.sig=colMeans(prime@assays$RNA@scale.data[jakstat, ])

stem=stem.bak
stem=stem[stem%in%rownames(prime@assays$RNA@scale.data)]
prime@meta.data$stem.sig=colMeans(prime@assays$RNA@scale.data[stem, ])

emt=emt.bak
emt=emt[emt%in%rownames(prime@assays$RNA@scale.data)]
prime@meta.data$emt.sig=colMeans(prime@assays$RNA@scale.data[emt, ])

check.genes=c("SOX2", "RB1", "TP53", "PTEN", "IL6") # jakstat.sig   stem.sig   emt.sig   c("SOX2", "RB1", "TP53", "PTEN", "IL6")
prime <- SetIdent(prime, value = "celltype")
DefaultAssay(prime)="RNA"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(prime, features = check.genes)
dev.off()
DefaultAssay(prime)="RNA"
png("temp.VlnPlot.png", width=10*100, height=5*100)
temp=VlnPlot(prime, features = check.genes, pt.size = 0, group.by="celltype", ncol = 5)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()
DefaultAssay(prime)="SCT"

# enrichment of Clike compared with all other cells
markers=prime.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)

# enrichment of Clike compared with TypeC
markerlist=readRDS("prime.celltype.pairwise.rds")
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
prime@meta.data$for_temp_color=as.numeric(factor(prime@meta.data$celltype, levels=c("Clike", "TypeC", "Luminal", "Basal")))
run_qusage_heatmap.seurat3(prime, nm = 'prime.celltype', kegg, my.seed=100)
run_qusage_heatmap.seurat3(prime, nm = 'prime.celltype', hall.list, my.seed=100)


### monocle2
{
  # for NE and CRPC individually
  require(monocle)
  # prime=readRDS("prime.rds")
  temp=prime
  # temp=subset(prime, subset=(pheno=="CRPC"))
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
  # tumor.markers=readRDS("prime.markers.celltype.rds")
  tumor.markers=readRDS("prime.markers.as.fig.cluster.byHarmony.rds")
  markers=tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[abs(markers$avg_log2FC)>log2(1.5), ]
  markers$foldChange=2^(markers$avg_log2FC)
  # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
  markers=markers[order(markers$foldChange, decreasing=TRUE), ]
  ordering.genes=unique(markers$gene)
  # ordering.genes=unique(rownames(markers))
  str(ordering.genes)
  # celltype DEGs
  # chr [1:866] "LCN2" "BPIFB1" "CHGA" "SCGB1A1" "PSCA" "SCGN" "FCGBP" ... 
  # cluster DEGs
  # chr [1:2476] "OLFM4" "PLP1" "SPARC" "VIM" "MMP7" "PMP22" "GPM6B" "TIMP3" ...

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
  str(pData(temp.m)$State)
  pData(temp.m)$State=rep(0, ncol(temp.m))
  temp.m@phenoData@data$State[temp.m@phenoData@data$celltype=="Clike"]=1
  temp.m <- orderCells(temp.m, root_state = temp.m@phenoData@data$State)

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

  saveRDS(temp.m, "monocle.prime.rds")
  # temp.m=readRDS("monocle.prime.rds")

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
  temp_genes <- c("PSCA", "CHGA", "SYP", "ENO2", "KLK3", "AR", "KRT4")
  temp_genes <- c("FEV", "ASCL1", "IRF7", "TP63", "CHGA", "SYP", "KLK3", "AR")
  pdf("dotplot.trajectory.temp.m.pdf")
  plot_genes_branched_pseudotime(temp.m[temp_genes,],
                         branch_point = 1,
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
  saveRDS(diff_test_res.temp.m, "diff_test_res.temp.m.prime.rds")
  # diff_test_res.temp.m=readRDS("diff_test_res.temp.m.prime.rds")
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


### deconvolution for validation in bulk data
{
  ### add prime's cell types to prim's cell types
  prime=readRDS("prime.rds")
  prim=readRDS("prim.rds")

  png("temp_prim_celltype_dimplot.png", width=5*100, height=5*100)
  DimPlot(prim, label = TRUE, group.by="celltype") + NoLegend()
  dev.off()

  png("temp_prime_celltype_dimplot.png", width=5*100, height=5*100)
  DimPlot(prime, label = TRUE, group.by="celltype") + NoLegend()
  dev.off()

  prim@meta.data$celltype[prim@meta.data$celltype=="B cell Lymph node"]="B/Plasma cell"
  temp.epi.cellnames=colnames(prim)[prim@meta.data$celltype%in%c("Basal cell", "Epithelial cell", "Luminal cell")]
  temp.epi.cellnames=temp.epi.cellnames[!temp.epi.cellnames%in%colnames(prime)]
  temp.epi.cellnames=colnames(prim)[!colnames(prim)%in%temp.epi.cellnames]
  prim=subset(prim, cells = temp.epi.cellnames)
  epitype=prime@meta.data$celltype
  # epitype=paste("Epi.", epitype, sep="") 
  names(epitype)=colnames(prime)
  prim@meta.data[names(epitype), "celltype"]=epitype
  table(prim@meta.data$celltype)
     # B/Plasma cell            Basal            Clike Endothelial cell
     #          4637              652              171             7288
     #    Fibroblast          Luminal       Macrophage        Mast cell
     #          2234            12893            11430             1745
     # Myofibroblast       Neutrophil           T cell            TypeC
     #          5183              325            19921              925

  # write percentage matrix to .csv
  pct.matrix=table(prim@meta.data$fig.patient, prim@meta.data$celltype)
  temp=apply(pct.matrix, 1, sum)
  pct.matrix=round(apply(pct.matrix, 2, function(x) x/temp), 4)
  temp.name=rownames(pct.matrix)
  temp.name
  #  [1] "CZK-T"      "DHB-T"      "HYQ"        "LPX"        "PXL"
  #  [6] "QLX-M"      "R-M"        "RSC1023"    "S6-TUMOR-M" "SHSCR"
  # [11] "SJ-M"       "T"          "ZQF-M"
  ### 9.5: N1.Gleason5
  GS=c( NA, NA, 9.5, 5, 9.5, 9.5, 2, 9.4, 5, 5, 2, NA, 3)
  pct.matrix=pct.matrix[order(GS, decreasing=TRUE), ]
  write.csv(pct.matrix, file="part4__prim__pct.matrix.csv", row.name=T)
  res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
  str(res)


  DefaultAssay(prim)="SCT"
  ### pairwise customized cluster's DEGs
  markerlist=list()
  celltype=unique(prim@meta.data$celltype)
  celltype=sort(celltype)
  prim = SetIdent(prim, value = "celltype")
  for(i in celltype)
  {
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          markerlist[[temp.name]]=FindMarkers(prim, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
        }
    }
  }
  # saveRDS(markerlist, "prim.celltype.pairwise.rds")
   

  # @data
  markerlist=readRDS("prim.celltype.pairwise.rds")
  # @scale.data
  # markerlist=readRDS("markerlist.prim.pairwise.initiating.types.scaled.data.rds")
  ######
  ### only consider epithelial cluster to find surface marker
  {
    # prim.bak=prim
    # prim=prim.bak
    # prim=subset(prim, celltype%in%c("TypeC", "Clike", "Luminal", "Basal"))
  }
  celltypes=sort(unique(prim@meta.data$celltype))
  markerlist.filter=list()
  for(j in celltypes)
  {
    print("celltype:")
    print(j)
    print("------------>>>>>>>>>>>")
    sp.markers=rownames(prim@assays$SCT@data)
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
  write.table(markerlist.filter[["Clike"]], file="marker.Clike.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["TypeC"]], file="marker.TypeC.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["Basal"]], file="marker.Basal.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["Luminal"]], file="marker.Luminal.txt", col.name=F, row.name=F, quote=F, sep='\t')
  clike.markers=markerlist.filter[["Clike"]]#[1:20]
  for(i in names(markerlist.filter))
  {
    markerlist.filter[[i]]=markerlist.filter[[i]][1:4]
    markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
  }
  markers=Reduce(union, markerlist.filter)
  str(markerlist.filter)
  str(markers)
  str(clike.markers)

  pdf(paste("temp.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/9 )
  temp.heatmap<-DoHeatmap( subset(prim, downsample=300), features = markers, group.by = "celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()


  ### only keep surface marker
  require(msigdbr)
  cc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
  cc = cc %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  cdgenes=sort(unique(cc$gene_symbol[cc$gs_name=="GO_CELL_SURFACE"]))
  table(clike.markers%in%cdgenes)
  clike.markers=clike.markers[clike.markers%in%cdgenes]

  prim@meta.data$celltype=factor(prim@meta.data$celltype, levels=c("T cell", "Myofibroblast", "Fibroblast", "Endothelial cell", "Mast cell", "Macrophage", "TypeC", "Clike", "Luminal", "Basal", "Neutrophil", "B/Plasma cell"))
  prim = SetIdent(prim, value = "celltype")
  pdf(paste("temp.prim.Clike.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(clike.markers)/7+3 )
  temp.heatmap<-DoHeatmap( subset(prim, downsample=300), features = clike.markers, group.by = "celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
  write.table(clike.markers, file="markerlist.prim.pairwise.celltype.txt", col.name=F, row.name=F, quote=F, sep='\t')

  tocheck=c("CBR1", "CBR3", "GSDMC", "MMP7", "F3", "CEACAM1", "LY6D", "FGFR3")
  pdf(paste("temp.prim.Clike.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(clike.markers)/7+3 )
  temp.heatmap<-DoHeatmap( subset(prim, downsample=300), features = tocheck, group.by = "celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()


  # @data-based markers
  # write.table(markers, file="marker.union.prim.pairwise.initiating.types.txt", col.name=F, row.name=F, quote=F, sep='\t')
  # @scale.data-based markers:
  write.table(markers, file="markerlist.prim.pairwise.celltype.txt", col.name=F, row.name=F, quote=F, sep='\t')

  # @data-based markers
  # markers=read.table("marker.union.prim.pairwise.initiating.types.txt", sep='\t')$V1
  # @scale.data-based markers:
  markers=read.table("markerlist.prim.pairwise.celltype.txt", sep='\t')$V1
  # make Cibersort's Signature matrix
  celltypes=sort(unique(prim@meta.data$celltype))
 # [1] "B/Plasma cell"    "Basal"            "Clike"            "Endothelial cell"
 # [5] "Fibroblast"       "Luminal"          "Macrophage"       "Mast cell"
 # [9] "Myofibroblast"    "Neutrophil"       "T cell"           "TypeC"
  ### create signature matrix "sigm"
  sigm=prim@assays$SCT@scale.data[markers, ]
  # sigm=prim@assays$SCT@data[markers, ]
  colnames(sigm)=prim@meta.data$celltype
  sigm=t(sigm)

  sigm=aggregate(sigm, by = list(rownames(sigm)), FUN = mean, na.rm = TRUE)
  sigm[1:3, 1:3]
  rownames(sigm)=sigm$"Group.1"
  sigm=sigm[, c(-1)]
  sigm=as.matrix(sigm)
  sigm=t(sigm)
  sigm[1:3, 1:3]
  # saveRDS(sigm, "prim.celltype.average.sig.matrix.rds")
  saveRDS(sigm, "prim.celltype.average.sig.matrix.scale.data.rds")
  # saveRDS(sigm, "prim.celltype.average.sig.matrix.2fc.dataslot.rds")

  ### use FARDEEP to deconvolution
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522071/
  # Rscript stage.R ../ ... /prim.celltype.average.sig.matrix.rds fardeep no no Clike Clike
  # Rscript survival2.R ../ ... /prim.celltype.average.sig.matrix.rds fardeep no no Clike Clike
  Rscript stage.R ../ ... /prim.celltype.average.sig.matrix.scale.data.rds fardeep no no Clike Clike;Rscript survival2.R ../ ... /prim.celltype.average.sig.matrix.scale.data.rds fardeep no no Clike Clike
  # Rscript stage.R ../ ... /prim.celltype.average.sig.matrix.2fc.dataslot.rds fardeep no no Clike Clike;Rscript survival2.R ../ ... /prim.celltype.average.sig.matrix.2fc.dataslot.rds fardeep no no Clike Clike

  # Rscript survival2.R ../ ... /marker.Clike.txt mean no no Clike ClikeMean
  # Rscript survival2.R ../ ... /marker.TypeC.txt mean no no TypeC TypeCMean
  # Rscript survival2.R ../ ... /marker.Luminal.txt mean no no Luminal LumiMean
  # Rscript survival2.R ../ ... /marker.Basal.txt mean no no Basal BasalMean

  # Rscript stage.R ../ ... /marker.Clike.txt mean no no Clike ClikeMean
  # Rscript stage.R ../ ... /marker.TypeC.txt mean no no TypeC TypeCMean
  # Rscript stage.R ../ ... /marker.Basal.txt mean no no Basal BasalMean
}


### deconvolution signature / pairwise DEGs of Clike in epithelial cell types
prime=readRDS("prime.rds")


png("temp_prime_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(prime, label = TRUE, group.by="celltype") + NoLegend()
dev.off()

table(prime@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   12893     925

### deconvolution for validation in bulk data
{
  DefaultAssay(prime)="SCT"
  ### pairwise customized cluster's DEGs
  markerlist=list()
  celltype=unique(prime@meta.data$celltype)
  celltype=sort(celltype)
  prime = SetIdent(prime, value = "celltype")
  for(i in celltype)
  {
    for(j in celltype)
    {
        if(i!=j)
        {
          temp.name=paste(i, j, sep="__")
          markerlist[[temp.name]]=FindMarkers(prime, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
        }
    }
  }
  saveRDS(markerlist, "prime.celltype.pairwise.rds")
   

  # @data
  markerlist=readRDS("prime.celltype.pairwise.rds")
  # @scale.data
  # markerlist=readRDS("markerlist.prime.pairwise.initiating.types.scaled.data.rds")
  ######
  celltypes=sort(unique(prime@meta.data$celltype))
  markerlist.filter=list()
  for(j in celltypes)
  {
    print("celltype:")
    print(j)
    print("------------>>>>>>>>>>>")
    sp.markers=rownames(prime@assays$SCT@data)
    for(i in names(markerlist))
    {
      temp=strsplit(i, "_")[[1]][1]
      if(temp==j)
      {
        markers=markerlist[[i]]
        markers=markers[markers$p_val_adj<0.05, ]
        markers=markers[order(markers$avg_log2FC, decreasing = TRUE), ]
        # markers=markers[markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
        markers=markers[markers$avg_log2FC>log2(2) , ]
        sp.markers=intersect(sp.markers, rownames(markers))
        print("pairwise:")
        print(i)
        print("length(sp.markers):")
        print(length(sp.markers))
      }
    }
    markerlist.filter[[j]]=sp.markers
  }
  str(markerlist.filter)
  write.table(markerlist.filter[["Clike"]], file="marker.Clike.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["TypeC"]], file="marker.TypeC.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["Basal"]], file="marker.Basal.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(markerlist.filter[["Luminal"]], file="marker.Luminal.txt", col.name=F, row.name=F, quote=F, sep='\t')
  for(i in names(markerlist.filter))
  {
    markerlist.filter[[i]]=markerlist.filter[[i]][1:5]
    markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
  }
  markers=Reduce(union, markerlist.filter)
  str(markers)

  pdf(paste("temp.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/9 )
  temp.heatmap<-DoHeatmap( subset(prime, downsample=300), features = markers, group.by = "celltype", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()



  # @data-based markers
  # write.table(markers, file="marker.union.prime.pairwise.initiating.types.txt", col.name=F, row.name=F, quote=F, sep='\t')
  # @scale.data-based markers:
  write.table(markers, file="markerlist.prime.pairwise.celltype.txt", col.name=F, row.name=F, quote=F, sep='\t')

  # @data-based markers
  # markers=read.table("marker.union.prime.pairwise.initiating.types.txt", sep='\t')$V1
  # @scale.data-based markers:
  markers=read.table("markerlist.prime.pairwise.celltype.txt", sep='\t')$V1
  # make Cibersort's Signature matrix
  celltypes=sort(unique(prime@meta.data$celltype))
 # [1] "B/Plasma cell"    "Basal"            "Clike"            "Endothelial cell"
 # [5] "Fibroblast"       "Luminal"          "Macrophage"       "Mast cell"
 # [9] "Myofibroblast"    "Neutrophil"       "T cell"           "TypeC"
  # create signature matrix "sigm"
  sigm=prime@assays$SCT@scale.data[markers, ]
  colnames(sigm)=prime@meta.data$celltype
  sigm=t(sigm)

  sigm=aggregate(sigm, by = list(rownames(sigm)), FUN = mean, na.rm = TRUE)
  sigm[1:3, 1:3]
  rownames(sigm)=sigm$"Group.1"
  sigm=sigm[, c(-1)]
  sigm=as.matrix(sigm)
  sigm=t(sigm)
  sigm[1:3, 1:3]
  # saveRDS(sigm, "prime.celltype.average.sig.matrix.rds")
  # saveRDS(sigm, "prime.celltype.average.sig.matrix.scale.data.rds")

  ### use FARDEEP to deconvolution
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522071/
  # Rscript stage.R ../ ... /prime.celltype.average.sig.matrix.rds fardeep no no Clike Clike
  # Rscript survival2.R ../ ... /prime.celltype.average.sig.matrix.rds fardeep no no Clike Clike

  # Rscript survival2.R ../ ... /marker.Clike.txt mean no no Clike ClikeMean
  # Rscript survival2.R ../ ... /marker.TypeC.txt mean no no TypeC TypeCMean
  # Rscript survival2.R ../ ... /marker.Luminal.txt mean no no Luminal LumiMean
  # Rscript survival2.R ../ ... /marker.Basal.txt mean no no Basal BasalMean

  # Rscript stage.R ../ ... /marker.Clike.txt mean no no Clike ClikeMean
  # Rscript stage.R ../ ... /marker.TypeC.txt mean no no TypeC TypeCMean
  # Rscript stage.R ../ ... /marker.Basal.txt mean no no Basal BasalMean
}
