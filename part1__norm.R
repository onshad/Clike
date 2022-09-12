### Seurat V3
source("../functions.seurat.R")
require(ggplot2)
require(RColorBrewer)
require(BoutrosLab.plotting.general)
require(dplyr)
require(canprot)
require(Seurat)
require(canprot
packageVersion("Seurat")
# [1] '4.0.1'

# working.directory:
setwd("/ ... ")

### read normalary PCa samples
### process from raw counts
data.directory="/ ... /counts/normal"
matrix.directory="/raw_feature_bc_matrix"
dir.raw=list.dirs(path = data.directory, full.names = TRUE, recursive = F)
dir.raw=paste(dir.raw, matrix.directory, sep="")
pos.1=regexpr(data.directory, dir.raw, perl=TRUE)
pos.1=attr(pos.1, "match.length")+2
pos.2=regexpr(matrix.directory, dir.raw, perl=TRUE)-1
sample.names=substring(dir.raw, pos.1, pos.2)
keep1=grepl("CBH", sample.names, fixed=TRUE)
keep2=grepl("GHL", sample.names, fixed=TRUE)
keep=keep1|keep2
dir.raw=dir.raw[keep]
sample.names=sample.names[keep]
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
normal=merge(x = inte.list[[1]], y = inte.list[2:length(inte.list)] )
normal@meta.data$fig.sample=as.character(normal@meta.data$orig.ident)
normal@meta.data$fig.patient=as.character(normal@meta.data$orig.ident)
normal@meta.data$fig.patient=substr(normal@meta.data$fig.patient, (nchar(normal@meta.data$fig.patient)-2), nchar(normal@meta.data$fig.patient))
normal@meta.data$fig.patient=as.factor(normal@meta.data$fig.patient)
normal@meta.data$fig.patient=as.character(normal@meta.data$fig.patient)
normal@meta.data$pheno=rep("Norm", nrow(normal@meta.data))
# normal@meta.data$pheno[normal@meta.data$fig.patient%in%c("QJZ")]="NEPC"
# pat2gleason=data.frame(patient=sort(unique(normal@meta.data$fig.patient)))
# pat2gleason$patient
# pat2gleason$gleason=factor(c(NA, NA, "5", "5", "5", "5", "2", "4", "5", "5", "2", NA, "3"))
# pat2gleason
# normal@meta.data$gleason=pat2gleason$gleason[match(normal@meta.data$fig.patient, pat2gleason$patient)]
# normal@meta.data$gleason=as.character(normal@meta.data$gleason)
# table(normal@meta.data$gleason)
table(normal@meta.data$pheno)
#  Norm
# 18627
# get prostate zone annotation
normal@meta.data$fig.zone=as.character(normal@meta.data$orig.ident)
normal@meta.data$fig.zone=substr(normal@meta.data$fig.zone, 1, (nchar(normal@meta.data$fig.zone)-6))
current.cluster.ids <- c("waizhou", "yixing", "zhongyang")
new.cluster.ids <- c("Peri", "Tran", "Cent") # c("Peripheral", "Transition", "Central")
normal@meta.data$fig.zone <- plyr::mapvalues(x = as.character(normal@meta.data$fig.zone), from = current.cluster.ids, to = new.cluster.ids)
table(normal@meta.data$fig.zone)
# Cent Peri Tran
# 6840 6488 5299
str(normal@meta.data)
# 'data.frame':   18627 obs. of  8 variables:
#  $ orig.ident  : chr  "waizhouyouCBH" "waizhouyouCBH" "waizhouyouCBH" "waizhouyouCBH" ...
#  $ nCount_RNA  : num  3183 2546 6801 6661 3346 ...
#  $ nFeature_RNA: int  1193 1378 2541 2386 1196 1965 1489 2168 1987 1007 ...
#  $ percent.mt  : num  11.2 10.4 14.9 14.7 29.5 ...
#  $ fig.sample  : chr  "waizhouyouCBH" "waizhouyouCBH" "waizhouyouCBH" "waizhouyouCBH" ...
#  $ fig.patient : chr  "CBH" "CBH" "CBH" "CBH" ...
#  $ pheno       : chr  "Norm" "Norm" "Norm" "Norm" ...
#  $ fig.zone    : chr  "Peri" "Peri" "Peri" "Peri" ...
normal=SCTransform(normal, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting

# check most variable genes
str(normal@assays$SCT@var.features)
# chr [1:3000] "RGS5" "COL1A1" "COL3A1" "PLA2G2A" "COL1A2" "SPP1" "SPARC" ...
# top10 <- head(VariableFeatures(normal, assay="SCT"), 10)
# # plot variable features with and without labels
# toplot<- VariableFeaturePlot(normal, assay = "SCT", selection.method="sctransform")
# toplot <- LabelPoints(plot = toplot, points = top10, repel = TRUE)
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()
# str(normal)

DefaultAssay(normal)="SCT"
normal <- RunPCA(normal, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(normal,  ndims = 100)
dev.off()

### PCA-based clustering
{
  # normal <- RunUMAP(normal, reduction = "pca", dims = 1:25, verbose = FALSE)
  # normal <- FindNeighbors(normal, reduction = "pca", dims = 1:25, verbose = FALSE)
  # normal <- FindClusters(normal, verbose = FALSE, resolution = 0.2)
  # table(normal@active.ident)

  # normal.markers <- FindAllMarkers(normal, assay="SCT")
  # saveRDS(normal.markers, "normal.markers.as.fig.cluster.byPCA.rds")
  # normal.markers=readRDS("normal.markers.as.fig.cluster.byPCA.rds")
}

### Harmony-based clustering
{
  require(harmony)
  DefaultAssay(normal)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  normal <- RunHarmony(normal, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  dev.off()
  Sys.time()

  # normal <- RunUMAP(normal, dims = 1:25, verbose = FALSE)
  normal <- RunUMAP(normal, reduction = "harmony", dims = 1:25, verbose = FALSE)
  normal <- FindNeighbors(normal, reduction = "harmony", dims = 1:25, verbose = FALSE)
  normal <- FindClusters(normal, verbose = FALSE, resolution = 0.2)
  table(normal@active.ident)
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 11897 10508  3867  2841  1996  1815  1803  1503  1345   933   861   555   511
#    13    14    15    16
#   322   316   206   198
  normal.markers <- FindAllMarkers(normal, assay="SCT")
  saveRDS(normal.markers, "normal.markers.as.fig.cluster.byHarmony.rds")
  # normal.markers =readRDS("normal.markers.as.fig.cluster.byHarmony.rds")
}


normal@meta.data$fig.cluster=normal@meta.data$seurat_clusters
normal <- SetIdent(normal, value = "fig.cluster")
pdf("normal_cluster_dimplot.pdf", width=5, height=4.5)
DimPlot(normal, label = TRUE, group.by="fig.cluster", pt.size=0.001, label.size=5, raster=TRUE)# + NoLegend()
dev.off()
png("temp_patient_dimplot.png", width=6*100, height=5*100)
DimPlot(normal, label = FALSE, group.by="fig.patient", pt.size=0.001)# + NoLegend()
dev.off()
png("temp_pheno_dimplot.png", width=6*100, height=5*100)
DimPlot(normal, label = FALSE, group.by="pheno", pt.size=0.001)# + NoLegend()
dev.off()
pdf("normal.dimplot.sample.pdf", height=4.5, width=6)
DimPlot(normal, label = FALSE, group.by="fig.sample", label.size=5, raster=TRUE)# + NoLegend()
dev.off()


### check celltype markers' expression to make sure every cluster's celltype
markers=normal.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[markers$avg_log2FC>log2(1.5), ]
setwd("..")
source("functions.seurat.R")
normal@meta.data$celltype=celltypeAnno(normal, markers)
setwd("./ ... ")
  # [1] "Luminal cell: 0.39"
  # [10] "Luminal cell: 0.24">NKX3-1+>epithelial???
  # [13] "Mast cell_cell201808: 1.83">TACSTD2+PSCA+KRT7/8/18/20+GATA3+>epithelial
  # [14] "Luminal cell: 0.51" >TSPAN1+TPPP3+EPCAM+>Ciliated epithelial cell
# epithelial cells (KRT5+, TACSTD2+, KRT19+), and not from bladder tissue ( EMA-, CDH2-,  UGP85-)
normal@meta.data$celltype[normal@meta.data$fig.cluster=="9"]="Epithelial cell"
normal@meta.data$celltype[normal@meta.data$fig.cluster=="12"]="Epithelial cell"
normal@meta.data$celltype[normal@meta.data$fig.cluster=="13"]="Epithelial cell"

normal@meta.data$celltype[normal@meta.data$fig.cluster=="0"]="Epithelial cell"
normal@meta.data$celltype[normal@meta.data$fig.cluster=="1"]="Epithelial cell"
normal@meta.data$celltype[normal@meta.data$fig.cluster=="4"]="Epithelial cell"
normal@meta.data$celltype[normal@meta.data$fig.cluster=="7"]="Epithelial cell"
# plot with cluster names
pdf("normal_celltype_dimplot.pdf", width=5, height=4)
DimPlot(normal, label = TRUE, group.by="celltype", pt.size=0.001, label.size=5, raster=TRUE) #+ NoLegend()
dev.off()


# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=normal.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.normal.cluster","csv",sep="."), row.name=T)


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
# check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "CLCN3", "SEMA7A", "DSG2", "VTCN1", "LTF", "CD74", "CFTR", "KRT4", "C3", "INSR", "CEACAM6", "CP", "TM4SF1", "RARRES1", "LCN2", "PPP1R1B", "EGFR", "SFRP1", "NGFR", "PRNP", "FOLH1", "ACPP")

DefaultAssay(normal)="RNA"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(normal, features = check.genes)
dev.off()
DefaultAssay(normal)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(normal, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(normal, "normal.rds")
# normal=readRDS("normal.rds")


# draw important genes' dotmap
imp.genes=c()
imp.genes=c("CD3D", "CD3E", "CD3G",    "CD4", "CD14", "CD163",    "KRT8", "KRT18", "AR",     "TPSB2", "TPSAB1", "MS4A2",     "ENG", "VWF", "PECAM1",  "ACTA2", "VIM", "CNN1",    "COL1A1", "FBLN1", "FN1"       )
imp.genes=c("CD3D", "CD3E", "CD3G",    "CD4", "CD14", "CD163",    "FTL", "CD55", "CD44",    "KRT8", "KRT18", "AR",     "TPSB2", "TPSAB1", "MS4A2",     "ENG", "VWF", "PECAM1",  "ACTA2", "VIM", "CNN1",    "COL1A1", "FBLN1", "FN1"  )
# DES: myofibroblast
imp.genes=unique(imp.genes)
imp.genes=imp.genes[length(imp.genes):1]# reverse
# delete celltype with less than 100 cells
normal=subset(normal, subset=celltype!="Myocyte")
normal@meta.data$fig.celltype=factor(normal@meta.data$celltype, levels=c("T cell", "Macrophage", "Epithelial cell", "Mast cell", "Endothelial cell", "Myofibroblast", "Fibroblast"))
pdf("normal_dotmap_celltype_markers.pdf", width = 18, height = 3);
DotPlot(normal, features = imp.genes, cols = c('white', 'blue'), dot.scale = 6, group.by="fig.celltype");
dev.off();


### further analysis focused on epithelial cells
table(normal@meta.data$celltype)
# Endothelial cell  Epithelial cell       Fibroblast       Macrophage
#             2464            12463              337             1983
#        Mast cell    Myofibroblast           T cell
#               42              838              500
table(is.na(normal@meta.data$celltype))
epitypes=c("Epithelial cell", "Luminal cell", "Basal cell")
normale=subset(normal, subset=(celltype%in%epitypes))
# str(normale)
normale@commands=list()


dim(normale)
# [1] 19599 12463
DefaultAssay(normale)="SCT"
str(normale@assays$SCT@var.features)
# chr [1:3000] "MSMB" "OLFM4" "LTF" "SCGB1A1" "MMP7" "KRT13" "SCGB3A1" ...
# normale <- FindVariableFeatures(normale, assay="SCT", nfeatures = 3000)
# str(normale@assays$SCT@var.features)
# top10 <- head(VariableFeatures(normale, assay="SCT"), 10)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(normale, assay = "SCT")
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# toplot=plot1 + plot2
# png("temp.png", width=10*100, height=5*100)
# print(toplot)
# dev.off()

normale=SCTransform(normale, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
normale <- RunPCA(normale, assay="SCT", verbose = FALSE, npcs = 100)
ElbowPlot(normale,  ndims = 100)
dev.off()

### clustering basd on PCA
{
  # DefaultAssay(normale)="SCT"
  # normale <- RunUMAP(normale, dims = 1:20, verbose = FALSE)
  # normale <- FindNeighbors(normale, dims = 1:20, verbose = FALSE)

  # DefaultAssay(normale)="SCT" 
  # normale <- FindClusters(normale, verbose = FALSE, resolution = 0.2)
  # table(normale@meta.data$seurat_clusters)
  # normale.markers <- FindAllMarkers(normale, assay="SCT")
  # saveRDS(normale.markers, "normale.markers.as.fig.cluster.byPCA.rds")
  # normale.markers =readRDS("normale.markers.as.fig.cluster.byPCA.rds")

}

### clustering basd on Harmony
{
  require(harmony)
  DefaultAssay(normale)
  # [1] "SCT"
  Sys.time()
  pdf("temp.runHarmony.pdf")
  normale <- RunHarmony(normale, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
  dev.off()
  Sys.time()

  # normal <- RunUMAP(normal, dims = 1:25, verbose = FALSE)
  normale <- RunUMAP(normale, reduction = "harmony", dims = 1:25, verbose = FALSE)
  normale <- FindNeighbors(normale, reduction = "harmony", dims = 1:25, verbose = FALSE)
  normale <- FindClusters(normale, verbose = FALSE, resolution = 0.2)
  table(normale@meta.data$seurat_clusters)
  #    0    1    2    3    4    5    6    7    8
  # 3705 3561 2925 1299  435  280  149   90   19
  normale.markers <- FindAllMarkers(normale, assay="SCT")
  # saveRDS(normale.markers, "normale.markers.as.fig.cluster.byHarmony.rds")
  # normale.markers =readRDS("normale.markers.as.fig.cluster.byHarmony.rds")
}

markers=normale.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.normale.cluster","csv",sep="."), row.name=T)

normale@meta.data$fig.cluster=normale@meta.data$seurat_clusters
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(normale, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
pdf("normale.dimplot.patient.pdf", height=4.5, width=5)
DimPlot(normale, label = FALSE, group.by="fig.patient", raster=TRUE)
dev.off()
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(normale, label = FALSE, group.by="gleason")# + NoLegend()
dev.off()
png("dimplot.pheno.png", height=5*100, width=5*100)
DimPlot(normale, label = FALSE, group.by="pheno")# + NoLegend()
dev.off()
png("dimplot.zone.png", height=5*100, width=5*100)
DimPlot(normale, label = FALSE, group.by="fig.zone")# + NoLegend()
dev.off()
pdf("normale.dimplot.sample.pdf", height=4.5, width=6)
DimPlot(normale, label = FALSE, group.by="fig.sample", label.size=5, raster=TRUE)# + NoLegend()
dev.off()

### analysis of normaltasis STOP here, because the small sample number (206 cells)
# saveRDS(normale, "normale.rds")
# normale=readRDS("normale.rds")

# cluster 4 ribosome(RPL, RPS, )+    such as RPL12 RPS14
# cluster 5 KRT5+
# cluster 6 > UPK2+PSCA+KRT20+exHighGATA3+ > bladder epithelial
# cluster 7 > PECAM1+ VIM+
# cluster 8 > KRT8+ EPCAM+ NCAM1+ ACSL4+ > neuroendocrine
# Delete clusters with < 100 cells
table(normale@meta.data$fig.cluster)
#    0    1    2    3    4    5    6    7    8
# 3705 3561 2925 1299  435  280  149   90   19
normale=subset(normale, cells = colnames(normale)[!normale@meta.data$fig.cluster%in%c("4", "6", "7")])
table(normale@meta.data$fig.cluster)
#    0    1    2    3    4    5    6    7    8
# 3705 3561 2925 1299    0  280    0    0   19
normale@meta.data$fig.cluster=as.character(normale@meta.data$fig.cluster)
table(normale@meta.data$fig.cluster)
#    0    1    2    3    5    8
# 3705 3561 2925 1299  280   19
table(is.na(normale@meta.data$fig.cluster))
normale@meta.data$fig.cluster=factor(normale@meta.data$fig.cluster, levels=sort(unique(normale@meta.data$fig.cluster)))
str(normale@meta.data$fig.cluster)

# check known epithelial cell types
check.genes=c("KLK3", "KRT14","KRT13", "KRT8", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2","NCAM1", "TP63", "AR", "NCAM1", "ACSL4", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1", "UPK2", "AR")#
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")
# # check genes 
# check.genes=c("CCL2", "CD74")

normale <- SetIdent(normale, value = "fig.cluster")
DefaultAssay(normale)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(normale, features = check.genes)
dev.off()
DefaultAssay(normale)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
temp=VlnPlot(normale, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# indi.plot(normale, tocheck="fig.patient")
# indi.plot(normale, tocheck="fig.cluster")

### re-cluster typeC cells to find Clike cluster
{
  normalb=subset(normale, subset=(fig.cluster%in%"3"))
  dim(normalb)
  # [1] 18521  1377
  DefaultAssay(normalb)="SCT"
  str(normalb@assays$SCT@var.features)
  # chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "MMP7" "KRT13" "LCN2" ...
  # normalb <- FindVariableFeatures(normalb, assay="SCT", nfeatures = 3000)
  # str(normalb@assays$SCT@var.features)
  # top10 <- head(VariableFeatures(normalb, assay="SCT"), 10)
  # # plot variable features with and without labels
  # plot1 <- VariableFeaturePlot(normalb, assay = "SCT")
  # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # toplot=plot1 + plot2
  # png("temp.png", width=10*100, height=5*100)
  # print(toplot)
  # dev.off()

  normalb=SCTransform(normalb, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  normalb <- RunPCA(normalb, assay="SCT", verbose = FALSE, npcs = 100)
  ElbowPlot(normalb,  ndims = 100)
  dev.off()

  ### clustering basd on PCA
  {
    # DefaultAssay(normalb)="SCT"
    # normalb <- RunUMAP(normalb, dims = 1:20, verbose = FALSE)
    # normalb <- FindNeighbors(normalb, dims = 1:20, verbose = FALSE)

    # DefaultAssay(normalb)="SCT" 
    # normalb <- FindClusters(normalb, verbose = FALSE, resolution = 0.2)
    # table(normalb@meta.data$seurat_clusters)
    # normalb.markers <- FindAllMarkers(normalb, assay="SCT")
    # saveRDS(normalb.markers, "normalb.markers.as.fig.cluster.byPCA.rds")
    # normalb.markers =readRDS("normalb.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(normalb)
    # [1] "SCT"
    Sys.time()
    normalb <- RunHarmony(normalb, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    dev.off()
    Sys.time()

    # normal <- RunUMAP(normal, dims = 1:25, verbose = FALSE)
    normalb <- RunUMAP(normalb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    normalb <- FindNeighbors(normalb, reduction = "harmony", dims = 1:25, verbose = FALSE)
    normalb <- FindClusters(normalb, verbose = FALSE, resolution = 0.2)
    table(normalb@meta.data$seurat_clusters)
      # 0   1   2   3
    # 445 414 306 212
    # normalb.markers <- FindAllMarkers(normalb, assay="SCT")
    # saveRDS(normalb.markers, "normalb.markers.as.fig.cluster.byHarmony.rds")
    # normalb.markers =readRDS("normalb.markers.as.fig.cluster.byHarmony.rds")
  }

  normalb@meta.data$fig.cluster=normalb@meta.data$seurat_clusters
  png("dimplot.cluster.png", height=5*100, width=5*100)
  DimPlot(normalb, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  dev.off()
  png("dimplot.patient.png", height=5*100, width=5*100)
  DimPlot(normalb, label = TRUE, group.by="fig.patient")
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  DimPlot(normalb, label = TRUE, group.by="gleason")
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  DimPlot(normalb, label = TRUE, group.by="pheno")
  dev.off()
  png("dimplot.zone.png", height=5*100, width=5*100)
  DimPlot(normalb, label = TRUE, group.by="fig.zone")
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
  # check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "CLCN3", "SEMA7A", "DSG2", "VTCN1", "LTF", "CD74", "CFTR", "KRT4", "C3", "INSR", "CEACAM6", "CP", "TM4SF1", "RARRES1", "LCN2", "PPP1R1B", "EGFR", "SFRP1", "NGFR", "PRNP", "FOLH1", "ACPP")


  normalb <- SetIdent(normalb, value = "fig.cluster")
  DefaultAssay(normalb)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(normalb, features = check.genes)
  dev.off()
  DefaultAssay(normalb)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=8*100)
  temp=VlnPlot(normalb, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  ### name clusters
  # current.cluster.ids <- c(0:4)
  # new.cluster.ids <- c("tC20", "tB21", "tL22", "tC23", "tC24")
  # # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
  # normalb@meta.data$epitype <- plyr::mapvalues(x = as.character(normalb@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  # png("dimplot.epitype.png", height=5*100, width=5*100)
  # DimPlot(normalb, label = TRUE, group.by="epitype", label.size=8)
  # dev.off()

  # btype=normalb@meta.data$epitype
  # names(btype)=rownames(normalb@meta.data)
  # table(btype)

  # indi.plot(normalb, tocheck="fig.patient")
  # indi.plot(normalb, tocheck="fig.cluster")

  ### DEG identification, across all clusters
  normalb <- SetIdent(normalb, value = "fig.cluster")
  normalb.markers <- FindAllMarkers(normalb, slot="data", assay="SCT")
  saveRDS(normalb.markers, "normalb.markers.cluster.rds")
  # normalb.markers=readRDS("normalb.markers.celltype.rds")

  # filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
  markers=normalb.markers
  markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
  markers$foldchange=2^(markers$avg_log2FC)
  markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
  write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.normalb.cluster","csv",sep="."), row.name=T)

  enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)



  # heatmap top 10 DEGs for every cluster
  markers=normalb.markers
  markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3 & markers$pct.2<0.1, ]
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
  # top10$cluster=factor(top10$cluster, levels=c("TypeC", "Luminal", "Basal"))
  top10=top10[order(top10$cluster), ]
  temp=subset(normalb, downsample=300)
  temp <- SetIdent(temp, value = "fig.cluster")
  # levels(temp)=c("TypeC", "Luminal", "Basal")
  tocheck=c(top10$gene, "KLK3", "AR")
  tocheck=c(tocheck, "KRT4", "TACSTD2", "PSCA")
  pdf(paste("heatmap.normalb.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=length(tocheck)/8+2)
  temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

  ### plot heatmap for special categories of DEGs
  setwd("..")
  subtype_markers=subtypeMarker(normalb)
  setwd("./ ... ")
  normalb.markers=readRDS("normalb.markers.celltype.rds")
  markers=normalb.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
  markers$foldchange=exp(markers$avg_log2FC)
  markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
  # for TFs, change log2(1.5) to log(1.2)
  # normalb <- ScaleData(normalb, features = rownames(normalb), assay="RNA" )
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
    # top10$cluster=factor(top10$cluster, levels=c("TypeC", "Luminal", "Basal"))
    top10=top10[order(top10$cluster), ]
    if(nrow(top10)>0)
    {
      DefaultAssay(normalb)="SCT"
      pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
      temp=subset(normalb, downsample=300)
      temp <- SetIdent(temp, value = "fig.cluster")
      # levels(temp)=c("TypeC", "Luminal", "Basal")
      temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
      # , subset=(fig.cluster%in%temp.clusters)
      print(temp.heatmap)
      dev.off()
      DefaultAssay(normalb)="SCT"
    }
  }


  saveRDS(normalb, "normalb.rds")
  # normalb=readRDS("normalb.rds")
}

### epitype annotation
current.cluster.ids <- unique(normale@meta.data$fig.cluster)
current.cluster.ids
# [1] 0 1 2 3 5 8
# Levels: 0 1 2 3 5 8
new.cluster.ids <- c("nL0", "nL1", "nB2", "nC3", "nB5", "nN8")
# tL: tumor luminal, tB: tumor basal, tCy: tumor CellCycle
normale@meta.data$epitype <- plyr::mapvalues(x = as.character(normale@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot
# normale@meta.data[names(btype), "epitype"]=btype
pdf("normale.dimplot.epitype.pdf", height=4.5, width=5)
DimPlot(normale, label = TRUE, group.by="epitype", label.size=8, raster=TRUE)
dev.off()


### celltype annotation
current.cluster.ids <- unique(normale@meta.data$epitype)
current.cluster.ids=sort(current.cluster.ids)
current.cluster.ids
#  [1] "nB2" "nB5" "nC3" "nL0" "nL1" "nN8"
new.cluster.ids <- c("Basal", "Basal", "TypeC", "Luminal", "Luminal", "NE")
# tLum: tumor luminal, tBI: tumor basal/intermediate, tCyc: tumor CellCycle
normale@meta.data$celltype <- plyr::mapvalues(x = as.character(normale@meta.data$epitype), from = current.cluster.ids, to = new.cluster.ids)
# add refined basal/typeC types into epitype slot
table( normale@meta.data$celltype)
  # Basal Luminal      NE   TypeC
  #  3205    7266      19    1299

### f2a
png("normale.celltype.dimplot.png", height=5*100, width=5*100)
DimPlot(normale, label = TRUE, group.by="celltype", label.size=5)
dev.off()
pdf("normale.celltype.dimplot.pdf", height=4.5, width=5)
temp.color=color[names(color)%in%unique(normale@meta.data$celltype)]
temp.oorder=oorder[oorder%in%unique(normale@meta.data$celltype)]
DimPlot(normale, label = TRUE, group.by="celltype", label.size=5, raster=TRUE, cols=temp.color, order=temp.oorder)
dev.off()
pdf("normale.cluster.dimplot.pdf", height=5, width=5)
DimPlot(normale, label = TRUE, group.by="fig.cluster", label.size=5, raster=TRUE)
dev.off()


### f2b
imp.genes=c(  "AR", "KRT8",    "KRT5",    "CHGA","SYP","ENO2",    "CYP1A1", "ARL14",    "KRT4", "PSCA")
# "CD55", "FGFR3",
imp.genes=c(  "KRT4", "TACSTD2", "PSCA",     "KRT5", "AR", "KRT8",   "CHGA","SYP","ENO2",   "NCAM1", "ACSL4")

DefaultAssay(normale)="SCT"
temp.color=color[names(color)%in%unique(normale@meta.data$celltype)]
temp.oorder=oorder[oorder%in%unique(normale@meta.data$celltype)]
normale@meta.data$fig.celltype=factor(normale@meta.data$celltype, levels=temp.oorder)
normale <- SetIdent(normale, value = "fig.celltype")
pdf("norme.VlnPlot.pdf", width=20, height=2)
temp=VlnPlot(normale, features = imp.genes, pt.size = 0, group.by=NULL, ncol = 11, cols=temp.color)
for(i in 1:length(imp.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 5, colour = "black", shape = 95)
}
print(temp)
dev.off()


### for percentage of cell types in different prostate zones, plot a barplot and oncoplot (package maftools), then combine them in adobe illustrator
## get annotation : zone + left/right + patientID (1/2)
# get prostate side (left/right) annotation
normale@meta.data$side=as.character(normale@meta.data$orig.ident)
normale@meta.data$side=substr(normale@meta.data$side, (nchar(normale@meta.data$side)-5), (nchar(normale@meta.data$side)-3))
current.cluster.ids <- c("you", "zuo")
new.cluster.ids <- c("right", "left") 
normale@meta.data$side <- plyr::mapvalues(x = as.character(normale@meta.data$side), from = current.cluster.ids, to = new.cluster.ids)
table(normale@meta.data$side)
# get prostate zone annotation
normale@meta.data$zone=as.character(normale@meta.data$orig.ident)
normale@meta.data$zone=substr(normale@meta.data$zone, 1, (nchar(normale@meta.data$zone)-6))
current.cluster.ids <- c("waizhou", "yixing", "zhongyang")
new.cluster.ids <- c("Peripheral zone", "Transition zone", "Central zone") # c("Peripheral", "Transition", "Central")
normale@meta.data$zone <- plyr::mapvalues(x = as.character(normale@meta.data$zone), from = current.cluster.ids, to = new.cluster.ids)
table(normale@meta.data$zone)
# get patientID (1/2) annotation
normale@meta.data$pID=as.character(normale@meta.data$orig.ident)
normale@meta.data$pID=substr(normale@meta.data$pID, (nchar(normale@meta.data$pID)-2), nchar(normale@meta.data$pID))
current.cluster.ids <- c("CBH", "GHL")
new.cluster.ids <- c("1", "2") # c("Peripheral", "Transition", "Central")
normale@meta.data$pID <- plyr::mapvalues(x = as.character(normale@meta.data$pID), from = current.cluster.ids, to = new.cluster.ids)
table(normale@meta.data$pID)
# zone + left/right + patientID (1/2)
normale@meta.data$fig.zone.side.p=paste(normale@meta.data$zone, normale@meta.data$pID, normale@meta.data$side, sep=" ")
## get side.zone.pID-celltype.percentage matrix
table( normale@meta.data$fig.zone.side.p )
#     Central zone 1 left    Central zone 1 right     Central zone 2 left
#                     947                     874                     598
#    Central zone 2 right  Peripheral zone 1 left Peripheral zone 1 right
#                    1516                     675                    1002
#  Peripheral zone 2 left Peripheral zone 2 right  Transition zone 1 left
#                     940                    1649                     468
# Transition zone 1 right  Transition zone 2 left Transition zone 2 right
#                     815                    1213                    1248
celltype.sta=table(normale@meta.data$fig.zone.side.p, normale@meta.data$celltype)
celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
barN=ncol(celltype.sta)
n_top = nrow(celltype.sta) 
pdf(paste("normale.barplot.celltypeInSample.pdf",sep=""))
# the width between paper edge and plot ( down,left,up,right : four directions)
par(mar = c(10, 4, 1, 8), xpd=T)
# coll=DiscretePalette(n_top, palette = "polychrome")
require("scales")
coll=hue_pal()(n_top)
barplot( celltype.sta, xlab="", ylab="Percentage", col=coll, names.arg=colnames(celltype.sta), las=2, cex.names=1) # , xaxt="n"
# facts=barN+1
end_point = barN #0.5 + n_top *facts-1
# text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
legend( barN+3,60, rownames(celltype.sta), fill=coll );
# text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
# text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
dev.off()


normale <- SetIdent(normale, value = "celltype")
DefaultAssay(normale)="SCT"

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
# check.genes=c("UPK1A", "CD55", "TGFBR3", "FGFR3", "IL1A", "ARL14", "FABP4", "SNX31", "CYP1A1", "CLCN3", "SEMA7A", "DSG2", "VTCN1", "LTF", "CD74", "CFTR", "KRT4", "C3", "INSR", "CEACAM6", "CP", "TM4SF1", "RARRES1", "LCN2", "PPP1R1B", "EGFR", "SFRP1", "NGFR", "PRNP", "FOLH1", "ACPP")

png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(normale, features = check.genes)
dev.off()
DefaultAssay(normale)="SCT"
png("temp.VlnPlot.png", width=10*100, height=12*100)
temp=VlnPlot(normale, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# saveRDS(normale, "normale.rds")
# normale=readRDS("normale.rds")


### Not employed --->>>
  ### CCA to identify epithelial cell types WithiN luminal tumor cells
  # normale=readRDS("normale.rds")
  normale=subset(normale, subset=celltype=="Luminal")
  anchors <- FindTransferAnchors(reference = normale, query = normale, 
      dims = 1:30, normalization.method="LogNormalize", reference.assay="SCT")
  predictions <- TransferData(anchorset = anchors, refdata = normale$celltype, 
      dims = 1:30)
  normale@meta.data$predtype=predictions$predicted.id 
  png("dimplot.predtype.png", height=5*100, width=5*100)
  DimPlot(normale, label = TRUE, group.by="predtype")
  dev.off()
### Not employed <<<---


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(normale@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(normale@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(normale@meta.data$fig.patient, normale@meta.data$celltype)
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
write.csv(pct.matrix, file="part4__normale__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(normale@meta.data$epitype)
  # Basal   Clike Luminal   TypeC
  #   652     171   13122     925
table(normale@meta.data$gleason)
#    2    3    4    5
# 2044  604  819 9906
pct.matrix=table(normale@meta.data$fig.patient, normale@meta.data$epitype)
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
write.csv(pct.matrix, file="part4__normale__pct.matrix.csv", row.name=T)
res <- t.test(as.numeric(pct.matrix[1:4, "TypeC"]), as.numeric(pct.matrix[5:10, "TypeC"]))
str(res)


# saveRDS(normale, "normale.rds")
# normale=readRDS("normale.rds")


### celltype-specific DEG identification
table(normale@meta.data$celltype)
  # Basal   Clike Luminal   TypeC
  #   652     171   12893     925
DefaultAssay(normale)="SCT"
### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(normale@meta.data$celltype)
celltype=sort(celltype)
normale = SetIdent(normale, value = "celltype")
for(i in celltype)
{
  for(j in celltype)
  {
      if(i!=j)
      {
        temp.name=paste(i, j, sep="__")
        markerlist[[temp.name]]=FindMarkers(normale, ident.1 = i, ident.2 = j, slot="data", assay="SCT")
      }
  }
}
# saveRDS(markerlist, "normale.celltype.pairwise.rds")
# markerlist=readRDS("normale.celltype.pairwise.rds")

sp.markerlist=list()
celltype=unique(normale@meta.data$celltype)
celltype=sort(celltype)
for(i in celltype)
{
  temp.genes=rownames(normale@assays$SCT@data)
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

temp=subset(normale, downsample=300)
sp.markerlist=sp.markerlist[c("Clike", "TypeC", "Basal", "Luminal")]
tocheck=Reduce(union, sp.markerlist)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("TypeC", "Basal", "Luminal", "NE")
pdf(paste("heatmap.normale.typeC.pdf", sep=""),  width=10, height=length(tocheck)/3)
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
normale <- SetIdent(normale, value = "celltype")
normale.markers <- FindAllMarkers(normale, slot="data", assay="SCT")
# saveRDS(normale.markers, "normale.markers.celltype.rds")
# normale.markers=readRDS("normale.markers.celltype.rds")

### plot heatmap
# normale.markers=readRDS("normale.indivSCT.markers.as.fig.cluster.rds")
# normale <- NormalizeData(normale, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA" )
# normale <- ScaleData(normale, features = rownames(normale), assay="RNA")

# filter DEGs, merge clusters, plot heatmap, and do enrichment analysis
markers=normale.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=2^(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.normale.cluster","csv",sep="."), row.name=T)

enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)


# heatmap top 10 DEGs for every cluster
markers=normale.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(2) & markers$pct.1>0.3 & markers$pct.2<0.1, ]
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
top10$cluster=factor(top10$cluster, levels=c("TypeC", "Luminal", "Basal"))
top10=top10[order(top10$cluster), ]
temp=subset(normale, downsample=300)
temp <- SetIdent(temp, value = "celltype")
levels(temp)=c("TypeC", "Luminal", "Basal")
tocheck=c(top10$gene, "KLK3", "AR")
tocheck=c(tocheck, "KRT4", "TACSTD2", "PSCA")
pdf(paste("heatmap.normale.cluster.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=length(tocheck)/8)
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### plot heatmap for special categories of DEGs
setwd("..")
subtype_markers=subtypeMarker(normale)
setwd("./ ... ")
normale.markers=readRDS("normale.markers.celltype.rds")
markers=normale.markers
markers=markers[markers$p_val_adj<0.05, ]
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
# for TFs, change log2(1.5) to log(1.2)
# normale <- ScaleData(normale, features = rownames(normale), assay="RNA" )
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
   top10$cluster=factor(top10$cluster, levels=c("TypeC", "Luminal", "Basal"))
  top10=top10[order(top10$cluster), ]
  if(nrow(top10)>0)
  {
    DefaultAssay(normale)="SCT"
    pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
    temp=subset(normale, downsample=300)
    temp <- SetIdent(temp, value = "celltype")
    levels(temp)=c("TypeC", "Luminal", "Basal")
    temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()
    DefaultAssay(normale)="SCT"
  }
}

# for too few cluster DExpressing genes, do vlnplot
tocheck=c("CCL2", "CXCL11", "CXCL1", "CXCL8", "CXCL17", "AR")
pdf("VlnPlot.chemokines_celltype.pdf", width = 10, height =10)
temp=VlnPlot(normale, features = tocheck, pt.size = 0, group.by="celltype", ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


# enrichment of Clike compared with all other cells
markers=normale.markers
markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
markers$foldchange=exp(markers$avg_log2FC)
markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
enrichment.gene.set.plot(markers$gene[markers$cluster=="Clike"], topNterm=20)

# enrichment of Clike compared with TypeC
markerlist=readRDS("normale.celltype.pairwise.rds")
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
normale@meta.data$for_temp_color=as.numeric(factor(normale@meta.data$celltype, levels=c("Luminal", "TypeC",  "Basal")))
# run_qusage_heatmap.seurat3(normale, nm = 'normale.celltype', kegg, my.seed=100)
run_qusage_heatmap.seurat3(normale, nm = 'normale.celltype', hall.list, my.seed=100)



### NOT RUN yet --->>>
### individually analysis of basal cells
bas13=subset(normale, subset=(celltype=="tBI"))
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
lum13=subset(normale, subset=(celltype%in%c("tLum", "tInt")))
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
table(normale@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
# to detect pairwise DEGs
normale@meta.data$celltype2compare=normale@meta.data$celltype
normale@meta.data[names(basal.anno), ]$celltype2compare=basal.anno
normale@meta.data[names(luminal.anno), ]$celltype2compare=luminal.anno
table(normale@meta.data$celltype2compare)
  # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
  # 213   165   109   115   219    87   632 11931  6763  4929   446


### get gleason annotation
table(normale@meta.data$fig.patient)
 #  53   54   55   56   59   62   71   73   74   75   76   77
 # 990  267 1335 5820 1677 4246  994  340 2617 1000  932 5391
current.cluster.ids <- sort(unique(normale@meta.data$fig.patient))
new.cluster.ids <- c("9", "10", "9", "9", "7", "7", "9", "9", "3+4", "9", "7", "9")
normale@meta.data$gleason <- plyr::mapvalues(x = as.character(normale@meta.data$fig.patient), from = current.cluster.ids, to = new.cluster.ids)
table(normale@meta.data$gleason)
  #  10   3+4     7     9
  # 267  2617  6855 15870
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(normale, label = FALSE, group.by="gleason", cols=c("green", "blue", "purple", "red"), order=c( "10", "9", "7", "3+4"), pt.size=0.1)
dev.off()

# saveRDS(normale, "normale.seurat3.indivSCTed.celltyped.rds")
# normale=readRDS("normale.seurat3.indivSCTed.celltyped.rds")
# normale=readRDS("normale.rds")

### pairwise customized cluster's DEGs
markerlist=list()
celltype=unique(normale@meta.data$celltype2compare)
normale = SetIdent(normale, value = "celltype2compare")
for(i in celltype)
{
  for(j in celltype)
  {
    if(i!=j)
    {
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(normale, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "normale.markerlist.pairwise.rds")
# markerlist=readRDS("normale.markerlist.pairwise.rds")


### pairwise fig.cluster's DEGs
markerlist=list()
celltype=unique(normale@meta.data$fig.cluster)
normale = SetIdent(normale, value = "fig.cluster")
for(i in celltype)
{
  for(j in celltype)
  {
    if(i!=j)
    {
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(normale, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "normale.markerlist.pairwise.fig.cluster.rds")
# markerlist=readRDS("normale.markerlist.pairwise.fig.cluster.rds")


### go to "_parts1and2_201216" to combine normal and PCa samples


### check basal/intermediate and luminal/intermediate signature in different-grade bulk samples
# identify signatures
table(normale@meta.data$celltype2compare)
  # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
  # 213   165   109   115   219    87   632 11931  6763  4929   446
png("temp_cluster_dimplot.png", width=5*100, height=5*100)
DimPlot(normale, label = TRUE, group.by="celltype2compare", pt.size=0.01)# + NoLegend()
dev.off()


### get gleason annotation
table(normale@meta.data$fig.patient)
 #  53   54   55   56   59   62   71   73   74   75   76   77
 # 990  267 1335 5820 1677 4246  994  340 2617 1000  932 5391
current.cluster.ids <- sort(unique(normale@meta.data$fig.patient))
new.cluster.ids <- c("9", "10", "9", "9", "7", "7", "9", "9", "3+4", "9", "7", "9")
normale@meta.data$gleason <- plyr::mapvalues(x = as.character(normale@meta.data$fig.patient), from = current.cluster.ids, to = new.cluster.ids)
table(normale@meta.data$gleason)
  #  10   3+4     7     9
  # 267  2617  6855 15870
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(normale, label = FALSE, group.by="gleason", cols=c("green", "blue", "purple", "red"), order=c( "10", "9", "7", "3+4"), pt.size=0.1)
dev.off()

###--->>>originating cell type based gene signature
# get originating-cell-type annotations
table(normale@meta.data$celltype2compare)
  # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
  # 213   165   109   115   219    87   632 11931  6763  4929   446
current.cluster.ids <- sort(unique(normale@meta.data$celltype2compare))
new.cluster.ids <- c("TypeC", "TypeC", "Basal", "TypeC", "Basal", "TypeC", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal")
normale@meta.data$ori.type <- plyr::mapvalues(x = as.character(normale@meta.data$celltype2compare), from = current.cluster.ids, to = new.cluster.ids)
table(normale@meta.data$ori.type)
  # Basal Luminal   TypeC
  #   328   24701     580
png("dimplot.type.originating.png", height=5*100, width=5*100)
DimPlot(normale, label = FALSE, group.by="ori.type", pt.size=0.1)
dev.off()

### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(normale@meta.data$ori.type)
  # Basal  Luminal   TypeC
  #   328    24701     580
pct.matrix=table(normale@meta.data$fig.patient, normale@meta.data$ori.type)
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
table(normale@meta.data$ori.type)
  # Basal Luminal   TypeC
  #   328   24701     580
ori.type=normale@meta.data$ori.type
names(ori.type)=rownames(normale@meta.data)
table(normal13@meta.data$ori.type)
  #      Basal Endothelial cell       Fibroblast          Luminal
  #        328             3212              408            24701
  # Macrophage        Mast cell    Myofibroblast           T cell
  #       2142             1443             1798             3135
  #      TypeC
  #        580
normal13@meta.data$ori.type=normal13@meta.data$celltype
normal13@meta.data[names(ori.type), ]$ori.type=ori.type
table(normal13@meta.data$ori.type)
cellnum=table(normal13@meta.data$ori.type)

markerlist=list()
celltype=unique(normal13@meta.data$ori.type)
normal13 = SetIdent(normal13, value = "ori.type")
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
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(normal13, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "normal13.markerlist.pairwise.ori.type.rds")
# markerlist=readRDS("normal13.markerlist.pairwise.ori.type.rds")

### get originating-cell-type signature
# markerlist=readRDS("normal13.markerlist.pairwise.ori.type.rds")
celltype=unique(normal13@meta.data$ori.type)
# 
signature=list()
for(j in celltype)
{
  geneset=c(rownames(normal13@assays$RNA))
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

# normal13 <- ScaleData(normal13, features = rownames(normal13), assay="RNA" )
# for( i in names(signature) )
normal13 = SetIdent(normal13, value = "ori.type")
for( i in c("Luminal", "TypeC") )
{
  temp = colMeans( normal13@assays$RNA@scale.data[signature[[i]] ,] )
  normal13@meta.data$for_temp_color=temp
  check.genes=c("for_temp_color")
  DefaultAssay(normal13)="RNA"
  png(paste("temp.featureplot", i, "png", sep="."), height=10*100, width=14*100)
  temp=FeaturePlot(normal13, features = check.genes)
  print(temp)
  dev.off()
  DefaultAssay(normal13)="SCT"
  png(paste("temp.VlnPlot", i, "png", sep="."), width=10*100, height=8*100)
  temp=VlnPlot(normal13, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2, assay = "RNA")
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
table(normale@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
png("dimplot.celltype.png", height=5*100, width=5*100)
DimPlot(normale, label = FALSE, group.by="celltype", pt.size=0.1)
dev.off()

### get percentage matrix of basal/int and luminal cells across patients 
### then go to Excel to plot the barplot
table(normale@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
pct.matrix=table(normale@meta.data$fig.patient, normale@meta.data$celltype)
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
table(normale@meta.data$celltype)
  # tBI  tCyc  tInt  tLum
  # 908   632   487 23582
celltype=normale@meta.data$celltype
names(celltype)=rownames(normale@meta.data)
table(normal13@meta.data$celltype)
    #   Basal cell Endothelial cell  Epithelial cell       Fibroblast
    #          964             3212              651              408
    # Luminal cell       Macrophage        Mast cell    Myofibroblast
    #        23994             2142             1443             1798
    #       T cell
    #         3135
normal13@meta.data$finer.type=normal13@meta.data$celltype
normal13@meta.data[names(celltype), ]$finer.type=celltype
table(normal13@meta.data$finer.type)


markerlist=list()
celltype=unique(normal13@meta.data$finer.type)
normal13 = SetIdent(normal13, value = "finer.type")
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
      markerlist[[paste(i, j, sep="__")]]=FindMarkers(normal13, ident.1 = i, ident.2 = j, assay="SCT")
    }
  }
}
# saveRDS(markerlist, "normal13.markerlist.pairwise.finer.type.rds")
# markerlist=readRDS("normal13.markerlist.pairwise.finer.type.rds")

### get originating-cell-type signature
# markerlist=readRDS("normal13.markerlist.pairwise.finer.type.rds")
celltype=unique(normal13@meta.data$finer.type)
# 
signature=list()
for(j in celltype)
{
  geneset=c(rownames(normal13@assays$RNA))
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

# normal13 <- ScaleData(normal13, features = rownames(normal13), assay="RNA" )
# for( i in names(signature) )
normal13 = SetIdent(normal13, value = "ori.type")
for( i in c("tBI", "tInt", "tLum") )
{
  temp = colMeans( normal13@assays$RNA@scale.data[signature[[i]] ,] )
  normal13@meta.data$for_temp_color=temp
  check.genes=c("for_temp_color")
  DefaultAssay(normal13)="RNA"
  png(paste("temp.featureplot", i, "png", sep="."), height=10*100, width=14*100)
  temp=FeaturePlot(normal13, features = check.genes)
  print(temp)
  dev.off()
  DefaultAssay(normal13)="SCT"
  png(paste("temp.VlnPlot", i, "png", sep="."), width=10*100, height=8*100)
  temp=VlnPlot(normal13, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2, assay = "RNA")
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
gleason.genes.up=gleason.genes.up[gleason.genes.up%in%rownames(normale@assays$SCT@scale.data)]
# gleason.genes.down=readRDS("gleason.genes.down.rds")
# gleason.genes.down=gleason.genes.down[gleason.genes.down%in%rownames(normale@assays$SCT@scale.data)]
# normale=ScaleData(normale)
normale@meta.data$bulk.gleason=colMeans(normale@assays$SCT@scale.data[gleason.genes.up, ])#-colMeans(normale@assays$SCT@scale.data[gleason.genes.down, ])
DefaultAssay(normale)="RNA"
png("normale.bulk.gleason.featureplot.png", width=4*100, height=5*100)
FeaturePlot(normale, features = "bulk.gleason", pt.size=0.01)
dev.off()
DefaultAssay(normale)="SCT"
# as celltype2compare
pdf("temp.VlnPlot.fig.cluster.pdf", width = 11, height =4)
tocheck=c("bulk.gleason")
temp=VlnPlot(normale, features = tocheck, ncol = 2,  pt.size = 0, group.by="fig.cluster", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()
# as fig.cluster
pdf("temp.VlnPlot.celltype2compare.pdf", width = 11, height =4)
tocheck=c("bulk.gleason")
temp=VlnPlot(normale, features = tocheck, ncol = 2,  pt.size = 0, group.by="celltype2compare", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()
# as celltype
pdf("temp.VlnPlot.celltype.pdf", width = 11, height =4)
tocheck=c("bulk.gleason")
temp=VlnPlot(normale, features = tocheck, ncol = 2,  pt.size = 0, group.by="celltype", assay = "RNA")
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


### survival and bulk-different-Gleason comparision

normale@assays$RNA@scale.data=matrix(0, 0, 0)
normale@commands=list()
normale@meta.data$basal.sig=NULL
normale@meta.data$cycle.sig=NULL
normale@meta.data$SCT_snn_res.0.5=NULL
normale@meta.data$SCT_snn_res.0.3=NULL
saveRDS(normale, "normale.seurat3.indivSCTed.celltyped.rds")
# normale=readRDS("normale.seurat3.indivSCTed.celltyped.rds")


### check indolent signature from Hansen

### NOT RUN yet <<<---

