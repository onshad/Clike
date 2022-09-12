
### Seurat V4
setwd("/ ... ")
source("../functions.seurat.R")
require(ggplot2)
require(RColorBrewer)
require(BoutrosLab.plotting.general)
require(dplyr)
require(canprot)
require(Seurat)
require(canprot)
require(GEOquery)
# [1] '1.0.0'
packageVersion("Seurat")
# [1] '4.0.1'
# options(future.globals.maxSize=300*1024^3)# memory 300GB
options(future.globals.maxSize=100*1024^3)# memory 100GB
# working.directory: 
setwd("/ ... ")

# saveRDS(epi.tumor, "epi.tumor.highCNV.harmony.adjed.rds")
# epi.tumor=readRDS("epi.tumor.highCNV.harmony.adjed.rds")

### read normal prostate, incidental, primary PCa, CRPC, and CRPC.metastasis samples
normale=readRDS("normale.rds")
incie=readRDS("incidente.rds")
prime=readRDS("prime.rds")
crpce=readRDS("crpce.rds")
rownames1=rownames(normale@assays$RNA@data)
rownames2=rownames(incie@assays$RNA@data)
rownames3=rownames(prime@assays$RNA@data)
rownames4=rownames(crpce@assays$RNA@data)
str(rownames1)
# chr [1:27984] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" "AL627309.3" ...
str(rownames4)
# chr [1:27984] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" "AL627309.3" ...
epi=merge(x = normale, y = list( incidente=incie, prime=prime, crpce=crpce))
str(epi)
# $ : chr [1:27984] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
# $ : chr [1:39463] "AAATCCTGTAAACGTACCGCAACATTA_1_1" "AAATCCTGTAAGAGTAG
table(epi@meta.data$epitype)
 #  crC6 crCy3  crL0  crL1  crL2 crNE5 iaC10 iaC11 iaC12   iB2   iL0   iL3   iL4
 #   70   795  1972  1755  1695   132   737   675   276   989  2235   920   270
 #  iL6  iNE5   nB2   nC3   nL0   nL1   nL4   nL5  tB21  tC20  tC23  tC24   tL0
 #  154   202  2879  1377  3475  3420   681   113   652   894   171    31  4069
 #  tL1  tL22   tL3   tL4   tL5
 # 3853   577  1659  1446  1289
table(is.na(epi@meta.data$epitype))
# FALSE
# 39463
table(epi@meta.data$celltype)
    # Basal CellCycle     Clike    Clike2   Luminal        NE     TypeC
    #  4520       795       517       675     29583       334      3039
table(is.na(epi@meta.data$celltype))
table(is.na(epi@meta.data$percent.mt))
epi=SCTransform(epi, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
# saveRDS(epi, "epi.fullsamples.rds")
# epi=readRDS("epi.fullsamples.rds")


# get cross-dataset pariwise DEGs
epi = SetIdent(epi, value = "epitype")
# epi=subset(epi, downsample=500)

# epi=subset(epi, subset=epitype!="tC24")
table(epi@meta.data$epitype)
# crC6 crCy3  crL0  crL1  crL2 crNE5 iaC10 iaC11 iaC12   iB2   iL0   iL3   iL4
#   70   795  1972  1755  1695   132   737   675   276   989  2235   920   270
#  iL6  iNE5   nB2   nC3   nL0   nL1   nL4   nL5  tB21  tC20  tC23  tC24   tL0
#  154   202  2879  1377  3475  3420   681   113   652   894   171    31  4069
#  tL1  tL22   tL3   tL4   tL5
# 3853   577  1659  1446  1289
dim(epi)
# [1] 21539 39463

markerlist.epi=list()
celltype=sort(unique(epi@meta.data$epitype))
length.celltype=length(celltype)
epi = SetIdent(epi, value = "epitype")
for(i in 1:(length.celltype-1))
{
  for(j in ((i+1):length.celltype))
  {
    print(paste("Comparing", i, "VS.", j, "(total", length.celltype, "clusters)", ":", sep=" "))
    ci=celltype[i]
    cj=celltype[j]
    print(paste("(", ci, "VS.", cj, ")", sep=" "))
    markerlist.epi[[paste(ci, cj, sep="__")]]=FindMarkers(epi, ident.1 = ci, ident.2 = cj, assay="SCT")
  }
}
# saveRDS(markerlist.epi, "epi.markerlist.pairwise.rds")
# markerlist.epi=readRDS("epi.markerlist.pairwise.rds")


### UMAP and find clusters based on 3 markersets 
### (2 intra-dataSet marker lists, and 1 inter-dataSet marker list)
markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
str(epi@assays$SCT@var.features)
# chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "LCN2" "SCGB1A1" "MMP7" ...
geneset=c()
for(i in markerlist.epi)
{
  geneset=union( geneset, rownames(i)[ (i$p_val_adj<0.05) & (abs(i$avg_log2FC)>log2(2))] )
}
str(geneset)
# chr [1:519] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "UPK1A" ... log2(4)
# chr [1:907] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(3)
# chr [1:2127] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(2)
# chr [1:4782] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(1.5)
DefaultAssay(epi)="SCT"
# var.features.bak=epi@assays$SCT@var.features
# epi@assays$SCT@var.features=var.features.bak
epi@assays$SCT@var.features=geneset
epi = SetIdent(epi, value = "epitype")
table(epi@active.ident)
# epi=subset(epi, downsample=500)
epi <- RunPCA(epi, assay="SCT", verbose = FALSE, npcs = 100)
pdf("temp.PCA.pdf")
ElbowPlot(epi,  ndims = 100)
dev.off()
str(epi@assays$SCT@var.features)
 # chr [1:4510] "PIGR" "KRT7" "S100A6" "LCN2" "CX3CL1" "SNCG" "UPK1A" ...
str(epi@reductions$pca@feature.loadings)
 # num [1:4510, 1:100] -0.0561 -0.0545 -0.0721 -0.0603 -0.0277 ...
 # - attr(*, "dimnames")=List of 2
 #  ..$ : chr [1:4510] "PIGR" "KRT7" "S100A6" "LCN2" ...
 #  ..$ : chr [1:100] "PC_1" "PC_2" "PC_3" "PC_4" ...

### PCA-based clustering
{
  epi <- RunUMAP(epi, reduction = "pca", dims = 1:25, verbose = FALSE)
  epi <- FindNeighbors(epi, reduction = "pca", dims = 1:25, verbose = FALSE)
  epi <- FindClusters(epi, verbose = FALSE, resolution = 0.2)
  table(epi@meta.data$seurat_clusters)

  # epi.markers <- FindAllMarkers(epi, assay="SCT")
  # saveRDS(epi.markers, "epi.markers.as.fig.cluster.byPCA.rds")
  # epi.markers=readRDS("epi.markers.as.fig.cluster.byPCA.rds")

  # saveRDS(epi, "epi.noHarmony.rds")
  # epi=readRDS("epi.noHarmony.rds")
}


### Harmony-based clustering
{
  # epi=readRDS("epi.noHarmony.rds")
  epi@assays$SCT@var.features=var.features.bak
  epi <- RunPCA(epi, assay="SCT", verbose = FALSE, npcs = 100)
  pdf("temp.PCA.pdf")
  ElbowPlot(epi,  ndims = 100)
  dev.off()

  require(harmony)
  DefaultAssay(epi)
  # [1] "SCT"
  Sys.time()
  table(is.na(epi@meta.data$fig.sample))
  pdf("temp.runHarmony.pdf")
  epi <- RunHarmony(epi, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  dev.off()
  Sys.time()

  epi <- RunUMAP(epi, reduction = "harmony", dims = 1:25, verbose = FALSE)
  epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:25, verbose = FALSE)
  epi <- FindClusters(epi, verbose = FALSE, resolution = 0.2)
  table(epi@meta.data$seurat_clusters)
  #    0    1    2    3    4    5    6    7    8    9   10   11   12
  # 9884 8891 8197 5822 2194 1228  836  833  683  295  224  198  178

  # epi.markers <- FindAllMarkers(epi, assay="SCT")
  # saveRDS(epi.markers, "epi.markers.as.fig.cluster.byHarmony.rds")
  # epi.markers =readRDS("epi.markers.as.fig.cluster.byHarmony.rds")

  # saveRDS(epi, "epi.harmony.rds")
  # epi=readRDS("epi.harmony.rds")
}


### Harmony-based clustering, separately consider batches within each phenotype
{
  # epi=readRDS("epi.noHarmony.rds")
  markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
  str(epi@assays$SCT@var.features)
  # chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "LCN2" "SCGB1A1" "MMP7" ...
  geneset=c()
  for(i in markerlist.epi)
  {
    geneset=union( geneset, rownames(i)[ (i$p_val_adj<0.05) & (abs(i$avg_log2FC)>log2(2))] )
  }
  str(geneset)
  # chr [1:519] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "UPK1A" ... log2(4)
  # chr [1:907] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(3)
  # chr [1:2127] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(2)
  # chr [1:4782] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(1.5)
  DefaultAssay(epi)="SCT"
  # var.features.bak=epi@assays$SCT@var.features
  # epi@assays$SCT@var.features=var.features.bak
  epi@assays$SCT@var.features=geneset

  epi <- RunPCA(epi, assay="SCT", verbose = FALSE, npcs = 100)
  pdf("temp.PCA.pdf")
  ElbowPlot(epi,  ndims = 100)
  dev.off()

  pheno=levels(epi@meta.data$fig.pheno)
  for(i in 1:length(pheno))
  {
    epi@meta.data[, paste("batch.", i, sep="")]=epi@meta.data$fig.sample
    temp=epi@meta.data[epi@meta.data$fig.pheno==pheno[i], paste("batch.", i, sep="")]
    temp=unique(temp)[1]
    epi@meta.data[epi@meta.data$fig.pheno!=pheno[i], paste("batch.", i, sep="")]=temp
  }

  require(harmony)
  DefaultAssay(epi)
  # [1] "SCT"
  Sys.time()
  table(is.na(epi@meta.data$fig.sample))
  pdf("temp.runHarmony.pdf")
  epi <- RunHarmony(epi, group.by.vars = paste("batch.", 1:length(pheno), sep=""), assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  dev.off()
  Sys.time()

  epi <- RunUMAP(epi, reduction = "harmony", dims = 1:25, verbose = FALSE)
  epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:25, verbose = FALSE)
  epi <- FindClusters(epi, verbose = FALSE, resolution = 0.2)
  table(epi@meta.data$seurat_clusters)

  # epi.markers <- FindAllMarkers(epi, assay="SCT")
  # saveRDS(epi.markers, "epi.markers.as.fig.cluster.byHarmony.rds")
  # epi.markers =readRDS("epi.markers.as.fig.cluster.byHarmony.rds")

  # saveRDS(epi, "epi.harmony.finerBatchCorrect.rds")
  # epi=readRDS("epi.harmony.finerBatchCorrect.rds")
}


### Harmony-based clustering (parameters adjusted), pair-wise DEGs as var genes 
{
  epi = readRDS("epi.fullsamples.rds")
  epi = SetIdent(epi, value = "epitype")
  table(epi@meta.data$epitype)
  # epi=subset(epi, downsample=150)
  # epi=subset(epi, downsample=300)
  markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
  str(epi@assays$SCT@var.features)
  # chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "LCN2" "SCGB1A1" "MMP7" ...
  geneset=c()
  for(i in markerlist.epi)
  {
    geneset=union( geneset, rownames(i)[ (i$p_val_adj<0.05) & (abs(i$avg_log2FC)>log2(1.5))] )
  }
  str(geneset)
  # chr [1:519] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "UPK1A" ... log2(4)
  # chr [1:907] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(3)
  # chr [1:2127] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(2)
  # chr [1:4782] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(1.5)
  DefaultAssay(epi)="SCT"
  # var.features.bak=epi@assays$SCT@var.features
  # epi@assays$SCT@var.features=var.features.bak
  epi@assays$SCT@var.features=geneset

  epi <- RunPCA(epi, assay="SCT", verbose = FALSE, npcs = 100)
  pdf("temp.PCA.pdf")
  temp=ElbowPlot(epi,  ndims = 100)
  print(temp)
  dev.off()

  require(harmony)
  DefaultAssay(epi)
  # [1] "SCT"
  # Sys.time()
  table(is.na(epi@meta.data$fig.sample))
  pdf("temp.runHarmony.pdf")
  # epi <- RunHarmony(epi, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25, theta=0, lambda=1, sigma=0.05, nclust=10, block.size=0.002) # for fullsample cells
  epi <- RunHarmony(epi, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25, theta=2, lambda=1, sigma=0.1, nclust=10, block.size=0.05) # for cnv-high cells
  dev.off()
  # Sys.time()

  epi <- RunUMAP(epi, reduction = "harmony", dims = 1:25, verbose = FALSE)
  epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:25, verbose = FALSE)
  epi <- FindClusters(epi, verbose = FALSE, resolution = 0.2)
  # table(epi@meta.data$seurat_clusters)

  epi@meta.data$fig.cluster=epi@meta.data$seurat_clusters
  png("dimplot.epitype.png", height=5*140, width=5*140)
  temp=DimPlot(epi, label = TRUE, group.by="epitype", label.size=4) # + NoLegend()
  print(temp)
  dev.off()
  png("dimplot.phenotype.png", height=5*100, width=5*100)
  temp=DimPlot(epi, label = FALSE, group.by="pheno")#, cols=c("blue", "red"))
  print(temp)
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  temp=DimPlot(epi, label = FALSE, group.by="gleason")#, cols=c("blue", "red"))
  print(temp)
  dev.off()
  png("dimplot.cluster.png", height=5*100, width=5*100)
  temp=DimPlot(epi, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  print(temp)
  dev.off()

  # epi.markers <- FindAllMarkers(epi, assay="SCT")
  # saveRDS(epi.markers, "epi.markers.as.fig.cluster.byHarmony.rds")
  # epi.markers =readRDS("epi.markers.as.fig.cluster.byHarmony.rds")

  # saveRDS(epi, "epi.harmony.adjed.rds")
  # epi=readRDS("epi.harmony.adjed.rds")
}

epi@meta.data$temp=NA
epi@meta.data$temp[epi@meta.data$pheno=="ADT"]="post-ADT"
temp.epi=subset(epi, subset=pheno=="ADT")
png("dimplot.png", height=5*100, width=6*100)
temp=DimPlot(temp.epi, label = FALSE, group.by="fig.celltype", label.size=8) # + NoLegend()
print(temp)
dev.off()

table(epi@meta.data$seurat_clusters)
epi@meta.data$fig.cluster=epi@meta.data$seurat_clusters
png("dimplot.epitype.png", height=5*140, width=5*140)
DimPlot(epi, label = FALSE, group.by="epitype", label.size=4) # + NoLegend()
dev.off()
pdf("dimplot.epitype.pdf", height=5, width=7)
DimPlot(epi, label = FALSE, group.by="epitype", label.size=4) # + NoLegend()
epi@meta.data$fig.celltype=epi@meta.data$celltype
epi@meta.data$fig.celltype[epi@meta.data$fig.celltype=="Clike2"]="Intermediate"
pdf("dimplot.celltype.pdf", height=5, width=6.5)
DimPlot(epi, label = FALSE, group.by="fig.celltype", label.size=4) # + NoLegend()
dev.off()
epi@meta.data$temp=NA
epi@meta.data$temp[epi@meta.data$epitype=="nL5"]="nL5"
epi@meta.data$temp[epi@meta.data$epitype=="iNE5"]="iNE5"
epi@meta.data$temp[epi@meta.data$epitype=="crNE5"]="crNE5"
epi@meta.data$temp[epi@meta.data$epitype=="iL3"]="iL3"
epi@meta.data$temp[epi@meta.data$epitype=="tC23"]="tC23"
epi@meta.data$temp[epi@meta.data$epitype=="iaC12"]="iaC12"
epi@meta.data$temp[epi@meta.data$epitype=="crC6"]="crC6"
png("dimplot.png", height=5*100, width=5*100)
DimPlot(epi, label = FALSE, group.by="temp")#, cols=c("blue", "red"))
dev.off()
epi@meta.data$temp=epi@meta.data$epitype
epi@meta.data$temp[!epi@meta.data$fig.pheno%in%c("Aggressive incidental PCa", "Indolent incidental PCa")]=NA

table(epi@meta.data$pheno)
  # ADT  Aggr  CRPC  Indo mHSPC  NEPC  Norm  Prim
  # 327  2869  6250  3589  6316   169 11945  7998
mHSPC=c("DHB-T", "HYQ", "PXL", "QLX-M", "RSC1023")# DHB-T wasn't annotated as HSPC before, but we later knew this tumor had stage M1 
epi@meta.data$temp <- epi@meta.data$pheno
epi@meta.data$temp[epi@meta.data$fig.patient%in%mHSPC]="mHSPC"
epi@meta.data$temp[epi@meta.data$temp=="Indo"]="Indolent incidental PCa"
epi@meta.data$temp[epi@meta.data$temp=="Aggr"]="Aggressive incidental PCa"
epi@meta.data$temp[epi@meta.data$temp=="Prim"]="Primary PCa"
epi@meta.data$temp[epi@meta.data$temp=="Norm"]="Normal prostate"
epi@meta.data$temp[epi@meta.data$temp=="ADT"]="post-ADT"
epi@meta.data$temp=factor(epi@meta.data$temp, levels=c("Normal prostate", "Indolent incidental PCa", "Aggressive incidental PCa", "Primary PCa", "mHSPC", "post-ADT", "CRPC", "NEPC"))
epi@meta.data$fig.pheno=epi@meta.data$temp
table(epi@meta.data$fig.pheno)
# Normal prostate   Indolent incidental PCa Aggressive incidental PCa
#           11945                      3589                      2869
#     Primary PCa                     mHSPC                  post-ADT
#            6970                      7344                       327
#            CRPC                      NEPC
#            6250                       169




png("dimplot.png", height=5*100, width=5*100)
DimPlot(epi, label = FALSE, group.by="temp")#, cols=c("blue", "red"))
dev.off()
# f1c
# png("dimplot.celltype.png", height=5*140, width=5*140)

pdf("f1c_dimplot.celltype.pdf")
DimPlot(epi, label = TRUE, group.by="celltype", label.size=4,  pt.size=0.01, raster = TRUE) # + NoLegend()
dev.off()
pdf("f1c_dimplot.celltype.pdf")
DimPlot(epi, label = TRUE, group.by="celltype", label.size=4,  pt.size=0.01, raster = TRUE) # + NoLegend()
dev.off()
png("dimplot.cluster.png", height=5*100, width=5*100)
DimPlot(epi, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
dev.off()
png("dimplot.patient.png", height=5*100, width=5*100)
DimPlot(epi, label = TRUE, group.by="fig.patient")
dev.off()
png("dimplot.zone.png", height=5*100, width=5*100)
DimPlot(epi, label = TRUE, group.by="fig.zone")
dev.off()
png("dimplot.phenotype.png", height=5*100, width=5*100)
DimPlot(epi, label = FALSE, group.by="pheno")#, cols=c("blue", "red"))
dev.off()
png("dimplot.gleason.png", height=5*100, width=5*100)
DimPlot(epi, label = FALSE, group.by="gleason")#, cols=c("blue", "red"))
dev.off()
epi@meta.data$temp=NA
epi@meta.data$temp[epi@meta.data$epitype=="tL0"]="tL0"
png("dimplot.png", height=5*100, width=5*100)
DimPlot(epi, label = FALSE, group.by="temp")#, cols=c("blue", "red"))
dev.off()
str(epi@meta.data$pheno)
table(epi@meta.data$pheno)
  # ADT  Aggr  CRPC  Indo mHSPC  NEPC  Norm  Prim
  # 154  1753  2124  1649  1629    78  2613  2419
table(is.na(epi@meta.data$pheno))


### barplot the percentage of each samples in cell types
{
  # give sample ID as phenotype + number
  pheno.sample=epi@meta.data[, c("fig.pheno", "fig.sample")]
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
  epi@meta.data$pheno.sample=pheno.sample[rownames(epi@meta.data), ]$pheno.sample

  celltype.sta=table(epi@meta.data$pheno.sample, epi@meta.data$fig.celltype)
  celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
  celltype.sta=celltype.sta[ c("Clike", "TypeC", "Basal", "NE"), !grepl("incidental|Normal", colnames(celltype.sta), perl=TRUE) ]
  barN=ncol(celltype.sta)
  n_top = nrow(celltype.sta) 
  pdf(paste("temp.barplot.celltypeInSample.pdf",sep=""))
  # the width between paper edge and plot ( down,left,up,right : four directions)
  par(mar = c(7, 4, 1, 8), xpd=T)
  # color=DiscretePalette(n_top, palette = "polychrome")
  require("scales")
  color=hue_pal()(n_top)
  barplot( celltype.sta, xlab="", ylab="Percentage", col=color, names.arg=colnames(celltype.sta), las=2, cex.names=0.5) # , xaxt="n"
  # facts=barN+1
  end_point = barN #0.5 + n_top *facts-1
  # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
  legend( barN+5,60, rownames(celltype.sta), fill=color );
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
  dev.off()
}


### boxplot Clike percentage across samples (mHSPC vs Primary, CRPC/NEPC vs ADT)
{
  celltype.sta=table(epi@meta.data$fig.sample, epi@meta.data$fig.celltype)
  celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
  # give phenotype annotation to the percentage table
  epi@meta.data$fig.pheno.simp=as.character(epi@meta.data$fig.pheno)
  epi@meta.data$fig.pheno.simp[epi@meta.data$fig.pheno.simp%in%c("CRPC", "NEPC")]="CRPC/NEPC"
  epi@meta.data$fig.pheno.simp=factor(epi@meta.data$fig.pheno.simp, levels=c("Normal prostate", "Indolent incidental PCa", "Aggressive incidental PCa", "Primary PCa", "mHSPC", "post-ADT", "CRPC/NEPC"))
  temp.df=data.frame(sample=epi@meta.data$fig.sample, pheno=epi@meta.data$fig.pheno.simp)
  temp.df=temp.df[!duplicated(temp.df), ]
  current.cluster.ids <- as.character(temp.df$sample)
  new.cluster.ids <- as.character(temp.df$pheno)
  temp <- plyr::mapvalues(x = colnames(celltype.sta), from = current.cluster.ids, to = new.cluster.ids)
  temp=factor(temp, levels=c("Normal prostate", "Indolent incidental PCa", "Aggressive incidental PCa", "Primary PCa", "mHSPC", "post-ADT", "CRPC/NEPC"))
  celltype.sta=data.frame(Percentage=celltype.sta["Clike", ], Phenotype=temp)


  # Clike percentage in aggressive vs indolent incidental samples
  # only consider primary, mHSPC, CRPC/NEPC
  df=celltype.sta[celltype.sta$Phenotype%in%c("Indolent incidental PCa", "Aggressive incidental PCa"), ]
  df$Phenotype=factor(as.character(df$Phenotype), levels=c("Indolent incidental PCa", "Aggressive incidental PCa"))
  require(ggplot2)
  pdf("f.3d.boxplot.ClikePCTacrossPheno.pdf", width=5, height=5)
  p <- ggplot(df, aes(x=Phenotype, y=Percentage)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  theme_classic() 
  # geom_jitter(shape=16, position=position_jitter(0.2))
  print(p)
  # boxplot( pct~pheno, data=df, ylab="Percentage of Clike cells" )
  #  main=paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test"
  res <- wilcox.test(df$Percentage[df$Phenotype=="Aggressive incidental PCa"], df$Percentage[df$Phenotype=="Indolent incidental PCa"])
  # legend( "topright", paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test", sep="") ,cex=0.8);
  dev.off()


  # Clike percentage in mHSPC vs primary PCa samples
  # only consider mHSPC and primary PCa samples
  df=celltype.sta[celltype.sta$Phenotype%in%c("Primary PCa", "mHSPC","post-ADT", "CRPC/NEPC"), ]
  df$Phenotype=factor(as.character(df$Phenotype), levels=c("Primary PCa", "mHSPC","post-ADT", "CRPC/NEPC"))
  require(ggplot2)
  pdf("f.3d.boxplot.ClikePCTacrossPheno.pdf", width=5, height=5)
  p <- ggplot(df, aes(x=Phenotype, y=Percentage)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  theme_classic() 
  # geom_jitter(shape=16, position=position_jitter(0.2))
  print(p)
  # boxplot( pct~pheno, data=df, ylab="Percentage of Clike cells" )
  #  main=paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test"
  res <- wilcox.test(df$Percentage[df$Phenotype=="Primary PCa"], df$Percentage[df$Phenotype=="mHSPC"])
  res <- wilcox.test(df$Percentage[df$Phenotype=="Primary PCa"], df$Percentage[df$Phenotype=="CRPC/NEPC"])
  res <- wilcox.test(df$Percentage[df$Phenotype=="mHSPC"], df$Percentage[df$Phenotype=="CRPC/NEPC"])
  res <- wilcox.test(df$Percentage[df$Phenotype=="Primary PCa"], df$Percentage[df$Phenotype=="post-ADT"])

  # legend( "topright", paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test", sep="") ,cex=0.8);
  dev.off()


 
}




### inferCNV by "part11__inferCNV.onlyEpithelial.V3.R"
{
  ### read inferCNV result (based on $SCT@counts), and classify cancer and normal cells

  ### downsampled epithelial cells's result with "nL0" and "nB2" cells as reference is under "../inferCNV/allcelltypes.sct.counts.out/", this is the best downsampled result
  ### downsample rules:   epi = SetIdent(epi, value = "epitype"); table(epi@meta.data$epitype); epi=subset(epi, downsample=300)

  # epi=readRDS("epi.fullsamples.rds")
  cnv=list()
  cnv[[1]]=read.table("../inferCNV/epitype.V4.sct.counts.out/infercnv.observations.txt")
  cnv[[2]]=read.table("../inferCNV/epitype.V4.sct.counts.out/infercnv.references.txt")
  for(i in 1:2)
  {
    cnv[[i]]=round((cnv[[i]]-1)^2, 3)
    cnv[[i]]=apply(cnv[[i]], 2, mean)
  }
  str(cnv)
  cnv=c(cnv[[1]], cnv[[2]])
  # saveRDS(cnv, "cnv.intensity.epitype.V4.sct.counts.rds")
  # cnv=readRDS("cnv.intensity.epitype.V4.sct.counts.rds")

  table(names(cnv) %in% colnames(epi@assays$SCT@counts))
  table(colnames(epi@assays$SCT@counts) %in% names(cnv))
  ### only consider the downsampled cells
  # epi=subset(epi, cells=names(cnv))

  # str(names(cnv)[!names(cnv)%in%colnames(epi@assays$SCT@counts)])
  # str(colnames(epi@assays$SCT@counts)[grepl("^G", colnames(epi@assays$SCT@counts), perl=TRUE)])

  # names.temp=gsub("\\.", "-", names(cnv))
  # str(names.temp[!names.temp%in%colnames(epi@assays$SCT@counts)])
  # # chr(0)
  # table(names.temp%in%colnames(epi@assays$SCT@counts))
  # #  TRUE
  # # 14439
  # table(colnames(epi@assays$SCT@counts)%in%names.temp)
  # #  TRUE
  # # 14439
  # names(cnv)=names.temp
  epi@meta.data$cnv.counts=rep(0, nrow(epi@meta.data))
  epi@meta.data[names(cnv), "cnv.counts"]=cnv
  # get median CNV intensity of primary tumor epithelial cells
  table(epi@meta.data$fig.pheno)
    # Normal prostate   Indolent incidental PCa Aggressive incidental PCa
    #         11945                      3589                      2869
    #   Primary PCa                     mHSPC                  post-ADT
    #          6970                      7344                       327
    #          CRPC                      NEPC
    #          6250                       169
  temp.median=median(epi@meta.data$cnv.counts[epi@meta.data$fig.pheno%in%c("mHSPC", "Primary PCa")])
  temp.median
  # 0.001974255
  quantile(epi@meta.data$cnv.counts[epi@meta.data$epitype%in%c("nL5")], 0.90)
  # 0.002349702
  quantile(epi@meta.data$cnv.counts[epi@meta.data$epitype%in%c("tC24")], 0.90)
  # 0.002051158

  # png("temp.VlnPlot.CNV.png", width=12*100, height=8*100)
  pdf("epi.VlnPlot.CNV_COUNTS.pdf", width=12, height=5)
  check.genes="cnv.counts"
  temp.fun <- function(x)
  {return(quantile(x, 0.90))}
  temp=VlnPlot(epi, features = check.genes, pt.size = 0, ncol = 1, group.by = "epitype")
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 6, colour = "black", shape = 95)
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = temp.fun, geom='point', size = 6, colour = "blue", shape = 95)
  }
  temp=temp+geom_hline(yintercept = temp.median, col = 'red', size = 1)# dashed line: , linetype = 'dotted'
  print(temp)
  dev.off()

}


### only consider high-CNV clusters ( 90th percentile > temp.median) as tumor cells and re-analyze them 
{ 
  temp.median=median(epi@meta.data$cnv.counts[epi@meta.data$fig.pheno%in%c("mHSPC", "Primary PCa")])
  temp.median
  # [1] 0.001974255
  print(dim(epi))
  # [1] 21539 39463
  ### identify whether a cluster contains tumor cells
  unicluster=unique(epi@meta.data$epitype)
  cluster.cnv=rep(1, length(unicluster))
  names(cluster.cnv)=unicluster
  for(i in unicluster)
  {
    cluster.iden = epi@meta.data$epitype==i
    cluster.cnv[i]=quantile(epi@meta.data$cnv.counts[cluster.iden], 0.90)
  }
  print("Normal clusters:")
  print(unicluster[cluster.cnv<temp.median])
  ### 95th percentile
  ### [1] "nL4" "nL1" "nL0" "nB2" "nC3"
  ### 90th percentile
  ### [1] "nL4"  "nL1"  "nL0"  "nB2"  "nC3"  "nL5"  "iL0"  "tC20" "tB21" "tC24"
  unicluster=unicluster[cluster.cnv>temp.median]
  table(epi@meta.data$epitype)
  epi.tumor=epi
  epi.tumor=subset(epi.tumor, subset=epitype%in%unicluster)
  print(dim(epi.tumor))
  ### 95th percentile
  ### [1] 21539 27631
  ### 90th percentile
  ### [1] 21539 23706
  # table(epi.tumor@meta.data$epitype)
  # table(epi.tumor@meta.data$fig.pheno)
}



# ### only consider high-CNV cells (> temp.median) as tumor cells and re-analyze them 
# {
#   temp.median=median(epi@meta.data$cnv.counts[epi@meta.data$fig.pheno%in%c("mHSPC", "Primary PCa")])
#   dim(epi)
#   # [1] 21539 39463
#   table(epi@meta.data$epitype)
#   epi.tumor=epi
#   epi.tumor@meta.data$temp.iden=epi.tumor@meta.data$cnv.counts>temp.median
#   epi.tumor=subset(epi.tumor, subset=temp.iden==TRUE)
#   dim(epi.tumor)
#   # [1] 21539 15437
#   table(epi.tumor@meta.data$epitype)
#   table(epi.tumor@meta.data$fig.pheno)
# }


### Harmony-based clustering (parameters adjusted), pair-wise DEGs as var genes 
{
  epi.tumor@meta.data$celltype[epi.tumor@meta.data$celltype=="Clike2"]="Intermediate"
  epi.tumor = SetIdent(epi.tumor, value = "epitype")
  table(epi.tumor@meta.data$epitype)
  # epi.tumor=subset(epi.tumor, downsample=150)
  # epi.tumor=subset(epi.tumor, downsample=300)
  markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
  str(epi.tumor@assays$SCT@var.features)
  # chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "LCN2" "SCGB1A1" "MMP7" ...
  geneset=c()
  print("Including clusters:")
  unicluster=unique(epi.tumor@meta.data$epitype)
  print(unicluster)
  for(i in names(markerlist.epi))
  {
    j=markerlist.epi[[i]]
    mn=strsplit(i, split="__")[[1]]
    m=mn[1]
    n=mn[2]
    if( (m%in%unicluster)&(n%in%unicluster) )
    {
      geneset=union( geneset, rownames(j)[ (j$p_val_adj<0.05) & (abs(j$avg_log2FC)>log2(1.5)) ] )
      #  )
    }
  }
  str(geneset)
  # chr [1:519] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "UPK1A" ... log2(4)
  # chr [1:907] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(3)
  # chr [1:2127] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(2)
  # chr [1:4782] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" ... log2(1.5)
  # chr [1:4367] "GATA3" "PIGR" "KRT7" ... including only part of clusters
  DefaultAssay(epi.tumor)="SCT"
  # var.features.bak=epi.tumor@assays$SCT@var.features
  # epi.tumor@assays$SCT@var.features=var.features.bak
  epi.tumor@assays$SCT@var.features=geneset

  epi.tumor <- RunPCA(epi.tumor, assay="SCT", verbose = FALSE, npcs = 100)
  pdf("temp.PCA.pdf")
  temp=ElbowPlot(epi.tumor,  ndims = 100)
  print(temp)
  dev.off()

  require(harmony)
  DefaultAssay(epi.tumor)
  # [1] "SCT"
  # Sys.time()
  table(is.na(epi.tumor@meta.data$fig.sample))
  pdf("temp.runHarmony.pdf")
  # epi.tumor <- RunHarmony(epi.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
  # epi.tumor <- RunHarmony(epi.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25, theta=0, lambda=1, sigma=0.05, nclust=10, block.size=0.002) # for fullsample cells
  # epi.tumor <- RunHarmony(epi.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25, theta=2, lambda=1, sigma=0.1, nclust=10, block.size=0.05) ### for cnv-high cells
  # epi.tumor <- RunHarmony(epi.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:35, theta=2, lambda=1, sigma=0.02, nclust=10, block.size=0.05) # for cnv-high clusters
  # epi.tumor <- RunHarmony(epi.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:35, theta=2, lambda=1, sigma=0.01, nclust=10, block.size=0.05) # for cnv-high clusters
  epi.tumor <- RunHarmony(epi.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE,   dims.use=1:50, theta=4, lambda=1, sigma=0.02, block.size=0.05, nclust=10     ) # for cnv-high clusters
  dev.off()
  # Sys.time()

  epi.tumor <- RunUMAP(epi.tumor, reduction = "harmony", dims = 1:50, verbose = FALSE)
  epi.tumor <- FindNeighbors(epi.tumor, reduction = "harmony", dims = 1:50, verbose = FALSE)
  epi.tumor <- FindClusters(epi.tumor, verbose = FALSE, resolution = 0.2)
  # resolution=0.2
  # table(epi.tumor@meta.data$seurat_clusters)

  epi.tumor@meta.data$fig.cluster=epi.tumor@meta.data$seurat_clusters
  png("dimplot.epitype.png", height=5*140, width=5*140)
  temp=DimPlot(epi.tumor, label = TRUE, group.by="epitype", label.size=4) # + NoLegend()
  print(temp)
  dev.off()
  png("dimplot.phenotype.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = FALSE, group.by="pheno")#, cols=c("blue", "red"))
  print(temp)
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = FALSE, group.by="gleason")#, cols=c("blue", "red"))
  print(temp)
  dev.off()
  png("dimplot.cluster.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  print(temp)
  dev.off()
  png("dimplot.celltype.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = TRUE, group.by="celltype", label.size=8) # + NoLegend()
  print(temp)
  dev.off()

  epi.tumor@meta.data$for_temp_compare=NA
  epi.tumor@meta.data$for_temp_compare[epi.tumor@meta.data$fig.pheno=="Indolent incidental PCa"]="Indolent incidental PCa"
  png("dimplot.indo.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = TRUE, group.by="for_temp_compare", label.size=8) # + NoLegend()
  print(temp)
  dev.off()

  epi.tumor@meta.data$for_temp_compare=NA # Aggressive 
  epi.tumor@meta.data$for_temp_compare[epi.tumor@meta.data$fig.pheno=="Aggressive incidental PCa"]="Aggressive incidental PCa"
  png("dimplot.aggr.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = TRUE, group.by="for_temp_compare", label.size=8) # + NoLegend()
  print(temp)
  dev.off()

  epi.tumor@meta.data$for_temp_compare=NA # NEPC
  epi.tumor@meta.data$for_temp_compare[epi.tumor@meta.data$fig.pheno=="NEPC"]="NEPC"
  png("dimplot.NEPC.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = TRUE, group.by="for_temp_compare", label.size=8) # + NoLegend()
  print(temp)
  dev.off()

  epi.tumor@meta.data$for_temp_compare=NA # post-ADT
  epi.tumor@meta.data$for_temp_compare[epi.tumor@meta.data$fig.pheno=="post-ADT"]="post-ADT"
  png("dimplot.ADT.png", height=5*100, width=5*100)
  temp=DimPlot(epi.tumor, label = TRUE, group.by="for_temp_compare", label.size=8) # + NoLegend()
  print(temp)
  dev.off()

  epi.tumor <- SetIdent(epi.tumor, value = "fig.cluster")
  # epi.tumor.markers <- FindAllMarkers(epi.tumor, assay="SCT")
  # saveRDS(epi.tumor.highCNV.markers, "epi.tumor.markers.as.fig.cluster.byHarmony.rds")
  # epi.tumor.markers =readRDS("epi.tumor.highCNV.markers.as.fig.cluster.byHarmony.rds")

  # saveRDS(epi.tumor, "epi.tumor.highCNV.harmony.adjed.rds")
  # epi.tumor=readRDS("epi.tumor.highCNV.harmony.adjed.rds")
}


### finer-cluster basal & typeC cells
{
  bas.tumor=subset(epi.tumor, subset=(fig.cluster%in%c("4", "9")))
  dim(bas.tumor)
  # [1] 20066  2325
  DefaultAssay(bas.tumor)="SCT"
  str(bas.tumor@assays$SCT@var.features)
  markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
  str(bas.tumor@assays$SCT@var.features)
  # chr [1:3000] "MSMB" "NPY" "OLFM4" "PLA2G2A" "LTF" "LCN2" "SCGB1A1" "MMP7" ...
  geneset=c()
  unicluster=unique(bas.tumor@meta.data$epitype)
  print("Including clusters:")
  print(unicluster)
  # only include previously identified Basal/C/Clike clusters from all phenotypes
  # unicluster= c("iaC11", "iaC10", "iaC12", "iB2", "tC23", "crC6")
  for(i in names(markerlist.epi))
  {
    j=markerlist.epi[[i]]
    mn=strsplit(i, split="__")[[1]]
    m=mn[1]
    n=mn[2]
    if( (m%in%unicluster)&(n%in%unicluster) )
    {
      geneset=union( geneset, rownames(j)[ (j$p_val_adj<0.05)& (abs(j$avg_log2FC)>log2(1.5)) ] )
    }
  }
  str(geneset)
  # chr [1:1690] "REG3A" "PADI3" "DHRS2" "LINC01764" "MXRA5" "DUOX2" "EMX2" ...
  DefaultAssay(bas.tumor)="SCT"
  # var.features.bak=bas.tumor@assays$SCT@var.features
  # bas.tumor@assays$SCT@var.features=var.features.bak

  bas.tumor=SCTransform(bas.tumor, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  bas.tumor@assays$SCT@var.features=geneset
  bas.tumor <- RunPCA(bas.tumor, assay="SCT", verbose = FALSE, npcs = 100)
  pdf("temp.PCA.pdf")
  temp=ElbowPlot(bas.tumor,  ndims = 100)
  print(temp)

  ### clustering basd on PCA
  {
    # DefaultAssay(bas.tumor)="SCT"
    # bas.tumor <- RunUMAP(bas.tumor, dims = 1:20, verbose = FALSE)
    # bas.tumor <- FindNeighbors(bas.tumor, dims = 1:20, verbose = FALSE)

    # DefaultAssay(bas.tumor)="SCT" 
    # bas.tumor <- FindClusters(bas.tumor, verbose = FALSE, resolution = 0.2)
    # table(bas.tumor@meta.data$seurat_clusters)
    # bas.tumor.markers <- FindAllMarkers(bas.tumor, assay="SCT")
    # saveRDS(bas.tumor.markers, "bas.tumor.markers.as.fig.cluster.byPCA.rds")
    # bas.tumor.markers =readRDS("bas.tumor.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(bas.tumor)
    # [1] "SCT"
    Sys.time()
    # bas.tumor <- RunHarmony(bas.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    bas.tumor <- RunHarmony(bas.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE,     dims.use=1:50, theta=4, lambda=1, sigma=0.02, block.size=0.05, nclust=10     ) # for cnv-high clusters
    dev.off()
    Sys.time()

    # bas.tumor <- RunUMAP(bas.tumor, reduction = "harmony", dims = 1:50, verbose = FALSE)
    bas.tumor <- FindNeighbors(bas.tumor, reduction = "harmony", dims = 1:50, verbose = FALSE)
    bas.tumor <- FindClusters(bas.tumor, verbose = FALSE, resolution = 0.2)
    # resolution=0.2
    table(bas.tumor@meta.data$seurat_clusters)
    #    0    1    2    3    4
    # 1294  834  558  225  147
    bas.tumor.markers <- FindAllMarkers(bas.tumor, assay="SCT")
    # saveRDS(bas.tumor.markers, "bas.tumor.markers.as.fig.cluster.byHarmony.rds")
    # bas.tumor.markers =readRDS("bas.tumor.markers.as.fig.cluster.byHarmony.rds")
    markers=bas.tumor.markers
    markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
    markers$foldchange=2^(markers$avg_log2FC)
    write.csv(markers[,c("foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.bas.tumor.markers.as.fig.cluster","csv",sep="."), row.name=T)
  }


  bas.tumor@meta.data$fig.cluster=bas.tumor@meta.data$seurat_clusters
  png("dimplot.cluster.png", height=5*100, width=5*100)
  temp=DimPlot(bas.tumor, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  print(temp)
  dev.off()
  png("dimplot.patient.png", height=5*100, width=5*100)
  temp=DimPlot(bas.tumor, label = TRUE, group.by="fig.patient")
  print(temp)
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  temp=DimPlot(bas.tumor, label = TRUE, group.by="gleason")
  print(temp)
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  temp=DimPlot(bas.tumor, label = TRUE, group.by="pheno")
  print(temp)
  dev.off()
  png("dimplot.oriCelltype.png", height=5*100, width=5*100)
  temp=DimPlot(bas.tumor, label = TRUE, group.by="celltype")
  print(temp)
  dev.off()

  # saveRDS(bas.tumor, "bas.tumor.rds")
  # bas.tumor=saveRDS("bas.tumor.rds")

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

  bas.tumor <- SetIdent(bas.tumor, value = "fig.cluster")
  DefaultAssay(bas.tumor)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(bas.tumor, features = check.genes)
  dev.off()
  DefaultAssay(bas.tumor)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=15*100)
  temp=VlnPlot(bas.tumor, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  ### name clusters
  current.cluster.ids <- c(0:4)
  new.cluster.ids <- c("Basal", "Intermediate", "TypeC", "Clike", "Basal")
  # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
  bas.tumor@meta.data$fig.celltype <- plyr::mapvalues(x = as.character(bas.tumor@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  png("dimplot.fig.celltype.png", height=5*100, width=5*100)
  DimPlot(bas.tumor, label = TRUE, group.by="fig.celltype", label.size=8)
  dev.off()

  # indi.plot(bas.tumor, tocheck="fig.patient")
  # indi.plot(bas.tumor, tocheck="fig.cluster")

  # saveRDS(bas.tumor, "bas.tumor.rds")
  # bas.tumor=readRDS("bas.tumor.rds")

  btype=bas.tumor@meta.data$fig.celltype
  names(btype)=rownames(bas.tumor@meta.data)
  print(table(btype))
  table(names(btype)%in%rownames(epi.tumor@meta.data))
  # TRUE
  # 3058

}


### finer-cluster cells in  cluster with both typeC and Clike cells
{
  typec.tumor=subset(bas.tumor, subset=(fig.cluster%in%c("2")))
  dim(typec.tumor)
  # [1] 16793   558
  DefaultAssay(typec.tumor)="SCT"
  str(typec.tumor@assays$SCT@var.features)
  # chr [1:3689] "VGLL1" "CLIC6" "GATA3" "UPK1B" "WFDC2" "GSTP1" "RARRES1" ...
  markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
  geneset=c()
  unicluster=unique(typec.tumor@meta.data$epitype)
  print("Including clusters:")
  print(unicluster)
  # only include previously identified Basal/C/Clike clusters from all phenotypes
  # unicluster= c("iaC11", "iaC10", "iaC12", "iB2", "tC23", "crC6")
  for(i in names(markerlist.epi))
  {
    j=markerlist.epi[[i]]
    mn=strsplit(i, split="__")[[1]]
    m=mn[1]
    n=mn[2]
    if( (m%in%unicluster)&(n%in%unicluster) )
    {
      geneset=union( geneset, rownames(j)[ (j$p_val_adj<0.05)& (abs(j$avg_log2FC)>log2(1.5)) ] )
    }
  }
  str(geneset)
  # chr [1:2956] "VGLL1" "CLIC6" "GATA3" "UPK1B" "WFDC2" "GSTP1" "RARRES1" ...
  DefaultAssay(typec.tumor)="SCT"
  typec.tumor=SCTransform(typec.tumor, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used to speed up model fitting
  typec.tumor@assays$SCT@var.features=geneset
  typec.tumor <- RunPCA(typec.tumor, assay="SCT", verbose = FALSE, npcs = 100)
  pdf("temp.PCA.pdf")
  temp=ElbowPlot(typec.tumor,  ndims = 100)
  print(temp)

  ### clustering basd on PCA
  {
    # DefaultAssay(typec.tumor)="SCT"
    # typec.tumor <- RunUMAP(typec.tumor, dims = 1:20, verbose = FALSE)
    # typec.tumor <- FindNeighbors(typec.tumor, dims = 1:20, verbose = FALSE)

    # DefaultAssay(typec.tumor)="SCT" 
    # typec.tumor <- FindClusters(typec.tumor, verbose = FALSE, resolution = 0.2)
    # table(typec.tumor@meta.data$seurat_clusters)
    # typec.tumor.markers <- FindAllMarkers(typec.tumor, assay="SCT")
    # saveRDS(typec.tumor.markers, "typec.tumor.markers.as.fig.cluster.byPCA.rds")
    # typec.tumor.markers =readRDS("typec.tumor.markers.as.fig.cluster.byPCA.rds")

  }

  ### clustering basd on Harmony
  {
    require(harmony)
    DefaultAssay(typec.tumor)
    # [1] "SCT"
    Sys.time()
    # typec.tumor <- RunHarmony(typec.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25) 
    typec.tumor <- RunHarmony(typec.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE,     dims.use=1:25, theta=4, lambda=1, sigma=0.02, block.size=0.05, nclust=10     ) # for cnv-high clusters
    dev.off()
    Sys.time()

    # typec.tumor <- RunUMAP(typec.tumor, reduction = "harmony", dims = 1:20, verbose = FALSE)
    typec.tumor <- FindNeighbors(typec.tumor, reduction = "harmony", dims = 1:25, verbose = FALSE)
    typec.tumor <- FindClusters(typec.tumor, verbose = FALSE, resolution = 0.2)
    # resolution=0.2
    table(typec.tumor@meta.data$seurat_clusters)
    #    0    1    2    3    4
    # 1294  834  558  225  147
    typec.tumor.markers <- FindAllMarkers(typec.tumor, assay="SCT")
    # saveRDS(typec.tumor.markers, "typec.tumor.markers.as.fig.cluster.byHarmony.rds")
    # typec.tumor.markers =readRDS("typec.tumor.markers.as.fig.cluster.byHarmony.rds")
    markers=typec.tumor.markers
    markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
    markers$foldchange=2^(markers$avg_log2FC)
    write.csv(markers[,c("foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.typec.tumor.markers.as.fig.cluster","csv",sep="."), row.name=T)
  }


  typec.tumor@meta.data$fig.cluster=typec.tumor@meta.data$seurat_clusters
  png("dimplot.cluster.png", height=5*100, width=5*100)
  temp=DimPlot(typec.tumor, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
  print(temp)
  dev.off()
  png("dimplot.patient.png", height=5*100, width=5*100)
  temp=DimPlot(typec.tumor, label = TRUE, group.by="fig.patient")
  print(temp)
  dev.off()
  png("dimplot.gleason.png", height=5*100, width=5*100)
  temp=DimPlot(typec.tumor, label = TRUE, group.by="gleason")
  print(temp)
  dev.off()
  png("dimplot.pheno.png", height=5*100, width=5*100)
  temp=DimPlot(typec.tumor, label = TRUE, group.by="pheno")
  print(temp)
  dev.off()
  png("dimplot.oriCelltype.png", height=5*100, width=5*100)
  temp=DimPlot(typec.tumor, label = TRUE, group.by="celltype")
  print(temp)
  dev.off()

  # saveRDS(typec.tumor, "typec.tumor.rds")
  # typec.tumor=readRDS("typec.tumor.rds")

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

  table(typec.tumor@meta.data$fig.cluster)
  #   0   1   2   3
  # 215 195  76  72
  typec.tumor <- SetIdent(typec.tumor, value = "fig.cluster")
  DefaultAssay(typec.tumor)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(typec.tumor, features = check.genes)
  dev.off()
  DefaultAssay(typec.tumor)="SCT"
  png("temp.VlnPlot.png", width=10*100, height=15*100)
  temp=VlnPlot(typec.tumor, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()

  ### name clusters
  current.cluster.ids <- c(0:3)
  new.cluster.ids <- c("TypeC", "Clike", "Basal", "Clike")
  # tB: tumor basal, tI: tumor intermediate, tH: tumor hillock, tC: tumor typeC
  typec.tumor@meta.data$fig.celltype <- plyr::mapvalues(x = as.character(typec.tumor@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)
  png("dimplot.fig.celltype.png", height=5*100, width=5*100)
  DimPlot(typec.tumor, label = TRUE, group.by="fig.celltype", label.size=8)
  dev.off()

  # indi.plot(typec.tumor, tocheck="fig.patient")
  # indi.plot(typec.tumor, tocheck="fig.cluster")

  # saveRDS(typec.tumor, "typec.tumor.rds")
  # typec.tumor=readRDS("typec.tumor.rds")

  ctype=typec.tumor@meta.data$fig.celltype
  names(ctype)=rownames(typec.tumor@meta.data)
  print(table(ctype))
  table(names(ctype)%in%rownames(epi.tumor@meta.data))

}


### to show basal typeC Clike intermediate clear in UMAP
### only ARL14 is notablly high, the others are not, so consider only show the primary PCa data
{
  sub.tumor=subset(epi.tumor, celltype%in%c("Clike", "TypeC", "Basal"))
  sub.tumor=subset(sub.tumor,  fig.pheno)
  # sub.tumor=SCTransform(sub.tumor, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi")# glmGamPoi is used
  str(sub.tumor@assays$SCT@var.features)
  # chr [1:3000] "BPIFB1" "SCGB1A1" "LCN2" "S100A2" "LTF" "FCGBP" "SPINK1" ...

  ### used pairwise DEGs as variable genes
  # {
  #   markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
  #   geneset=c()
  #   unicluster=unique(sub.tumor@meta.data$epitype)
  #   print("Including clusters:")
  #   print(unicluster)
  #   # only include previously identified Basal/C/Clike clusters from all phenotypes
  #   # unicluster= c("iaC11", "iaC10", "iaC12", "iB2", "tC23", "crC6")
  #   for(i in names(markerlist.epi))
  #   {
  #     j=markerlist.epi[[i]]
  #     mn=strsplit(i, split="__")[[1]]
  #     m=mn[1]
  #     n=mn[2]
  #     if( (m%in%unicluster)&(n%in%unicluster) )
  #     {
  #       geneset=union( geneset, rownames(j)[ (j$p_val_adj<0.05)& (abs(j$avg_log2FC)>log2(1.5)) ] )
  #     }
  #   }
  #   str(geneset)
  #   # chr [1:1556] "REG3A" "PADI3" "DHRS2" "LINC01764" "MXRA5" "DUOX2" "EMX2" ...
  #   sub.tumor@assays$SCT@var.features=geneset
  #   DefaultAssay(sub.tumor)="SCT"
  # }
 
  str(sub.tumor@assays$SCT@var.features)
  # chr [1:3000] "BPIFB1" "SCGB1A1" "LCN2" "S100A2" "LTF" "FCGBP" "SPINK1" ...

  ### Harmony-based clustering
  {
    sub.tumor <- RunPCA(sub.tumor, assay="SCT", verbose = FALSE, npcs = 100)
    pdf("temp.PCA.pdf")
    ElbowPlot(sub.tumor,  ndims = 100)
    dev.off()

    require(harmony)
    DefaultAssay(sub.tumor)
    # [1] "SCT"
    Sys.time()
    table(is.na(sub.tumor@meta.data$fig.sample))
    pdf("temp.runHarmony.pdf")
    sub.tumor <- RunHarmony(sub.tumor, group.by.vars = "fig.sample", assay.use="SCT", plot_convergence = TRUE, dims.use=1:25)
    dev.off()
    Sys.time()

    sub.tumor <- RunUMAP(sub.tumor, reduction = "harmony", dims = 1:25, verbose = FALSE)
    sub.tumor <- FindNeighbors(sub.tumor, reduction = "harmony", dims = 1:25, verbose = FALSE)
    sub.tumor <- FindClusters(sub.tumor, verbose = FALSE, resolution = 0.2)
    table(sub.tumor@meta.data$seurat_clusters)
    sub.tumor@meta.data$fig.cluster=sub.tumor@meta.data$seurat_clusters
    #   0   1   2   3   4   5   6
    # 730 429 429 309 233 101  12

    # sub.tumor.markers <- FindAllMarkers(sub.tumor, assay="SCT")
    # saveRDS(sub.tumor.markers, "sub.tumor.markers.as.fig.cluster.byHarmony.rds")
    # sub.tumor.markers =readRDS("sub.tumor.markers.as.fig.cluster.byHarmony.rds")

    # saveRDS(sub.tumor, "sub.tumor.harmony.rds")
    # sub.tumor=readRDS("sub.harmony.rds")
  }

  # keep colors consistent with all-epithelial-celltypes image
  oorder=c("Clike", "TypeC",  "Basal" , "Intermediate", "Luminal (low-grade PCa)", "Luminal (high-grade PCa)", "Luminal (CRPC)", "Luminal (Cell Cycling)", "Neuroendocrine")
  # oorder=oorder[length(oorder):1]
  # get colour
  require("scales")
  color=hue_pal()(length(oorder))
  color=color[length(oorder):1]# 
  # specify Clike's colour
  color[length(oorder)]="#FF6600"

  # only keep Clike, TypeC, and Basal
  color=color[length(oorder):(length(oorder)-2)]
  oorder=oorder[3:1]

  pdf("sub.tumor.celltype_dimplot.pdf", width=5.5, height=5)
  DimPlot(sub.tumor, label = FALSE, group.by="celltype", raster = FALSE, order=oorder, cols=color, pt.size=0.1, )# + NoLegend()
  dev.off()

  pdf("sub.tumor.cluster_dimplot.pdf", width=5.5, height=5)
  DimPlot(sub.tumor, label = FALSE, group.by="fig.cluster", raster = TRUE)
  # + NoLegend()  
  dev.off()

  pdf("sub.tumor.dimplot.phenotype.pdf", height=5, width=7)
  temp=DimPlot(sub.tumor, label = FALSE, group.by="fig.pheno", raster = TRUE)#, cols=c("blue", "red"))
  print(temp)
  dev.off()

  # check known epithelial cell types
  check.genes=c("KRT4", "PSCA", "KRT5", "IL1A", "ARL14", "UPK1A", "CD55",  "FGFR3"  )
  # # prostate/urinary-cancer-specific marker
  # check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
  # # tumor-initiating markers
  # check.genes=c("NKX3-1", "LY6D")
  # # zhangbo's markers
  # check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")

  DefaultAssay(sub.tumor)="SCT"
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(sub.tumor, features = check.genes)
  dev.off()
  DefaultAssay(sub.tumor)="SCT"
  png("temp.VlnPlot.png", width=5*100, height=20*100)
  sub.tumor=SetIdent(sub.tumor, value="fig.celltype")
  temp=VlnPlot(sub.tumor, features = check.genes, pt.size = 0, group.by="celltype", ncol = 2)
  for(i in 1:length(check.genes)) 
  {
    temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
  }
  print(temp)
  dev.off()


}


current.cluster.ids <- c(0:11)
new.cluster.ids <- c("Luminal (high-grade PCa)", "Luminal (high-grade PCa)", "Luminal (high-grade PCa)", "Luminal (low-grade PCa)", "Basal/intermediate", "Luminal (CRPC)", "Luminal (CRPC)", "Luminal (Cell Cycling)", "Luminal (CRPC)", "Basal/intermediate", "Neuroendocrine", "Luminal (high-grade PCa)")
epi.tumor@meta.data$fig.celltype <- plyr::mapvalues(x = as.character(epi.tumor@meta.data$fig.cluster), from = current.cluster.ids, to = new.cluster.ids)


### add basal/intermediate's finer annotation
epi.tumor@meta.data[names(btype), "fig.celltype"]=btype
table(epi.tumor@meta.data$fig.celltype)
#          Basal                  TypeC                  Clike
#           1441                    558                    225
#   Intermediate    Luminal (low-grade PCa)   Luminal (high-grade PCa)
#            834                   2939                  12761
# Luminal (CRPC) Luminal (Cell Cycling)         Neuroendocrine
#           3869                    976                    103

### add basal/intermediate's finer annotation
epi.tumor@meta.data[names(ctype), "fig.celltype"]=ctype
table(epi.tumor@meta.data$fig.celltype)
#          Basal                    TypeC                    Clike
#           1517                      215                      492
#   Intermediate  Luminal (low-grade PCa) Luminal (high-grade PCa)
#            834                     2939                    12761
# Luminal (CRPC)   Luminal (Cell Cycling)           Neuroendocrine
#           3869                      976                      103

epi.tumor@meta.data$fig.celltype=factor(epi.tumor@meta.data$fig.celltype, levels=c("Clike", "TypeC", "Basal", "Intermediate", "Luminal (low-grade PCa)", "Luminal (high-grade PCa)", "Luminal (CRPC)", "Luminal (Cell Cycling)", "Neuroendocrine"))
png("dimplot.celltype.png", height=5*100, width=5*100)
DimPlot(epi.tumor, label = TRUE, group.by="fig.celltype", label.size=5)
dev.off()


### f2e
table(epi.tumor@meta.data$celltype)
    # Basal CellCycle     Clike   Luminal        NE     TypeC
    #   652       795       241     18315       132       925
oorder=c("Clike", "TypeC",  "Basal" , "Luminal", "CellCycle", "NE")
# oorder=oorder[length(oorder):1]
# get colour
require("scales")
color=hue_pal()(length(oorder))
color=color[length(oorder):1]# 
names(color)=oorder[length(oorder):1]
# specify Clike's colour
# color[length(oorder)]="#FF6600"
pdf("f.2e_celltype_dimplot.pdf", width=5, height=4.5)
DimPlot(epi.tumor, label = FALSE, group.by="celltype", raster = TRUE, order=oorder, cols=color)
# + NoLegend()  
dev.off()
png("f.1e_celltype_dimplot.png", width=6.5*100, height=5*100)
DimPlot(epi.tumor, label = FALSE, group.by="fig.celltype", raster = FALSE, order=oorder, cols=color, pt.size=0.5)# + NoLegend()
dev.off()
pdf("f.2es_cluster_dimplot.pdf", width=5.5, height=5)
DimPlot(epi.tumor, label = FALSE, group.by="fig.cluster", raster = TRUE)
# + NoLegend()  
dev.off()


### f2f
pdf("f.2f_dimplot.phenotype.pdf", height=5, width=7)
temp=DimPlot(epi.tumor, label = FALSE, group.by="fig.pheno", raster = TRUE)#, cols=c("blue", "red"))
print(temp)
dev.off()

png("f.2f_dimplot.phenotype.png", height=5*100, width=7*100)
temp=DimPlot(epi.tumor, label = FALSE, group.by="fig.pheno")#, cols=c("blue", "red"))
print(temp)
dev.off()

pdf("f.2f_dimplot.gleason.pdf", height=5, width=5)
temp=DimPlot(epi.tumor, label = FALSE, group.by="gleason", raster = TRUE)#, cols=c("blue", "red"))
print(temp)
dev.off()

png("f.2f_dimplot.gleason.png", height=5*100, width=5*100)
temp=DimPlot(epi.tumor, label = FALSE, group.by="gleason")#, cols=c("blue", "red"))
print(temp)
dev.off()

pdf("f.2f_dimplot.cluster.pdf", height=5, width=5)
temp=DimPlot(epi.tumor, label = TRUE, group.by="fig.cluster", label.size=8, raster = TRUE) # + NoLegend()
print(temp)
dev.off()

png("f.2f_dimplot.cluster.png", height=5*100, width=5*100)
temp=DimPlot(epi.tumor, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
print(temp)
dev.off()

### only keep NE
temp.epi.tumor=subset(epi.tumor, fig.pheno=="NEPC") 
pdf("f.2h_NEsample_dimplot.pdf", width=7, height=5)
DimPlot(temp.epi.tumor, label = FALSE, group.by="fig.celltype", pt.size=0.01, raster = FALSE)# + NoLegend()
dev.off()


epi.tumor=SetIdent(epi.tumor, value="fig.celltype")
epi.tumor.markers <- FindAllMarkers(epi.tumor, assay="SCT")
# saveRDS(epi.tumor.markers, "epi.tumor.markers.as.celltype.byHarmony.rds")
# epi.tumor.markers =readRDS("epi.tumor.markers.as.celltype.byHarmony.rds")

# epi.tumor.markers =readRDS("epi.tumor.markers.as.celltype.byHarmony.rds")
markers=epi.tumor.markers
markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5) & markers$pct.1>0.1, ]
markers$foldchange=2^(markers$avg_log2FC)
write.csv(markers[,c("foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("markers.epi.tumor.markers.as.fig.cluster","csv",sep="."), row.name=T)

# & markers$pct.1>0.3 & markers$pct.2<0.1
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = (avg_log2FC) )
top10$cluster=factor(top10$cluster, levels=c("Basal/intermediate", "Luminal (low-grade PCa)", "Luminal (high-grade PCa)", "Luminal (CRPC)", "Luminal (Cell Cycling)", "Neuroendocrine"))
top10=top10[order(top10$cluster), ]
temp=subset(epi.tumor, downsample=300)
temp <- SetIdent(temp, value = "fig.celltype")
levels(temp)=c("Basal/intermediate", "Luminal (low-grade PCa)", "Luminal (high-grade PCa)", "Luminal (CRPC)", "Luminal (Cell Cycling)", "Neuroendocrine")
tocheck=c(top10$gene, "KLK3", "AR", "FOXA1", "PSCA", "KRT4", "AMACR")

pdf(paste("heatmap.epi.tumor.celltype.top10.pdf", sep=""),  width=(length(unique(markers$cluster))/4+5), height=length(tocheck)/8)
temp.heatmap<-DoHeatmap( temp, features = tocheck, slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()



### plot heatmap for special categories of DEGs
{
  # epi.tumor.bak=readRDS("epi.tumor.highCNV.harmony.adjed.rds")
  # epi.tumor=epi.tumor.bak
  setwd("..")
  subtype_markers=subtypeMarker(epi.tumor)
  setwd("./ ... ")
  epi.tumor.markers =readRDS("epi.tumor.markers.as.celltype.byHarmony.rds")
  markers=epi.tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
  markers$foldchange=exp(markers$avg_log2FC)
  markers=markers[ markers$avg_log2FC>log2(1.5) & markers$pct.1>0.3, ]
  # for TFs, change log2(1.5) to log(1.2)
  # epi.tumor <- ScaleData(epi.tumor, features = rownames(epi.tumor), assay="RNA" )
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
      DefaultAssay(epi.tumor)="SCT"
      pdf(paste("temp.heatmap.", j ,"_fig.cluster.pdf", sep=""), width=(length(unique(temp.markers$cluster))/4+5), height=nrow(top10)/10+3)
      temp=subset(epi.tumor, downsample=300)
      temp <- SetIdent(temp, value = "fig.celltype")
      temp.heatmap<-DoHeatmap(temp, features = top10$gene, slot="scale.data", assay="SCT", angle = 0, size=3) 
      # , subset=(fig.cluster%in%temp.clusters)
      print(temp.heatmap)
      dev.off()
      DefaultAssay(epi.tumor)="SCT"
    }
  } 
}


### f2f
### pick up 2 genes per group
DefaultAssay(epi.tumor)="SCT"
tocheck=c( "IL1RN", "IL1A", "PMAIP1",  "PSCA", "KRT4" ,   "KRT5",   "KLK3", "AR",  "TOP2A", "STMN1",    "ENO2"  )
# tocheck=c("TACSTD2", "GPRC5A", "UPK1A", "CD55", "FGFR3", "TGFBR3", "IL1A", "ARL14", "FABP4")
# tocheck=c( "KRT5", "KRT19", "PSCA", "KRT4" ,   "KLK3", "ACPP",    "AMACR", "AR",     "TOP2A", "STMN1",      "ENO2", "FOXA1",  "cnv.counts", "nFeature_RNA") # MMP7, AMACR
# pdf("f1_epithelial.dotmap_celltype_markers.pdf", width = 12, height = 2);
# DotPlot(temp, features = tocheck, cols = c('white', 'blue'), dot.scale = 6, group.by="fig.celltype");
# dev.off();
# tocheck=c("KRT4", "PSCA", "KRT5", "IL1A", "ARL14", "UPK1A", "CD55",  "FGFR3"  )

pdf("epi.tumor_VlnPlot_celltype.pdf", width = 5, height =20)
# oorder=c("Clike", "TypeC",  "Basal" , "Intermediate", "Luminal (low-grade PCa)", "Luminal (high-grade PCa)", "Luminal (CRPC)", "Luminal (Cell Cycling)", "Neuroendocrine")
temp.oorder=oorder
temp.color=color
# specify Clike's colour
epi.tumor@meta.data$factor.fig.celltype=factor(epi.tumor@meta.data$celltype, levels=oorder)
temp=VlnPlot(epi.tumor, features = tocheck, pt.size = 0, group.by="factor.fig.celltype", ncol = 1, cols = temp.color)
for(i in 1:length(tocheck)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
# adjust the border line's width, not succesfully
# temp$layers[[1]]$aes_params$size = 0.1
# adjust = 0.1 # refer to ggplot2.tidyverse.org/reference/geom_violin.html
print(temp)
dev.off()
DefaultAssay(epi.tumor)="SCT"


# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")

DefaultAssay(epi.tumor)="SCT"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(epi.tumor, features = check.genes)
dev.off()
DefaultAssay(epi.tumor)="SCT"
png("temp.VlnPlot.png", width=5*100, height=20*100)
epi.tumor=SetIdent(epi.tumor, value="fig.celltype")
temp=VlnPlot(epi.tumor, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


### barplot the percentage of each samples in cell types
{
  # give sample ID as phenotype + number
  pheno.sample=epi.tumor@meta.data[, c("fig.pheno", "fig.sample")]
  pheno.sample$fig.sample=as.character(pheno.sample$fig.sample)
  # str( pheno.sample)
  # 'data.frame':   23706 obs. of  2 variables:
  #  $ fig.pheno : Factor w/ 8 levels "Normal prostate",..: 2 2 2 2 2 2 2 2 2 2 ...
  #  $ fig.sample: chr  "waizhouyouGJC" "waizhouyouGJC" "waizhouyouGJC" 
  pheno.sample=pheno.sample[order(pheno.sample$fig.pheno, pheno.sample$fig.sample, decreasing=c(FALSE, FALSE)), ]
  unique.sample=unique(pheno.sample$fig.sample)
  pheno.sample$fig.sample=factor(pheno.sample$fig.sample, levels=unique.sample)
  pheno.sample$fig.sample.id=as.numeric(pheno.sample$fig.sample)
  pheno.sample$pheno.sample=paste(as.character(pheno.sample$fig.pheno), pheno.sample$fig.sample.id)
  unique.pheno.sample=unique(pheno.sample$pheno.sample)
  pheno.sample$pheno.sample=factor(pheno.sample$pheno.sample, levels=unique.pheno.sample)
  epi.tumor@meta.data$pheno.sample=pheno.sample[rownames(epi.tumor@meta.data), ]$pheno.sample

  celltype.sta=table(epi.tumor@meta.data$pheno.sample, epi.tumor@meta.data$fig.celltype)
  celltype.sta=apply(celltype.sta, 1, function(x) 100*x/sum(x))
  barN=ncol(celltype.sta)
  n_top = nrow(celltype.sta) 
  pdf(paste("temp.barplot.celltypeInSample.pdf",sep=""))
  # the width between paper edge and plot ( down,left,up,right : four directions)
  par(mar = c(7, 4, 1, 8), xpd=T)
  # color=DiscretePalette(n_top, palette = "polychrome")
  require("scales")
  color=hue_pal()(n_top)
  barplot( celltype.sta, xlab="", ylab="Percentage", col=color, names.arg=colnames(celltype.sta), las=2, cex.names=0.5) # , xaxt="n"
  # facts=barN+1
  end_point = barN #0.5 + n_top *facts-1
  # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
  legend( barN+10,60, rownames(celltype.sta), fill=color );
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
  dev.off()

  # only consider Clike cells
  str(celltype.sta)
  pdf(paste("temp.barplot.celltypeInSample.pdf",sep=""))
  # the width between paper edge and plot ( down,left,up,right : four directions)
  par(mar = c(7, 4, 1, 8), xpd=T)
  # color=DiscretePalette(n_top, palette = "polychrome")
  require("scales")
  color=hue_pal()(n_top)
  barplot( celltype.sta, xlab="", ylab="Percentage", col=color, names.arg=colnames(celltype.sta), las=2, cex.names=0.5) # , xaxt="n"
  # facts=barN+1
  end_point = barN #0.5 + n_top *facts-1
  # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
  legend( barN+10,60, rownames(celltype.sta), fill=color );
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
  # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
  dev.off()

}


### specially mark mHSPC samples in featureplot
table(epi.tumor@meta.data$pheno)
  # ADT  Aggr  CRPC  Indo mHSPC  NEPC  Norm  Prim
  # 327  2869  6250  3589  6316   169 11945  7998
unique(epi.tumor@meta.data$fig.patient)
#  [1] "CBH"        "GHL"        "GJC"        "SQB"        "ZYM"
#  [6] "CZK-T"      "DHB-T"      "HYQ"        "LPX"        "PXL"
# [11] "QLX-M"      "R-M"        "RSC1023"    "S6-TUMOR-M" "SHSCR"
# [16] "SJ-M"       "T"          "ZQF-M"      "DGY-PG"     "HYL"
# [21] "QJZ"        "XYM"
mHSPC=c("DHB-T", "HYQ", "PXL", "QLX-M", "RSC1023")# DHB-T wasn't annotated as HSPC before, but we later knew this tumor had stage M1 
epi@meta.data$temp <- epi@meta.data$pheno
epi@meta.data$temp[epi@meta.data$fig.patient%in%mHSPC]="mHSPC"
epi@meta.data$temp[epi@meta.data$temp=="Indo"]="Indolent incidental PCa"
epi@meta.data$temp[epi@meta.data$temp=="Aggr"]="Aggressive incidental PCa"
epi@meta.data$temp[epi@meta.data$temp=="Prim"]="Primary PCa"
epi@meta.data$temp[epi@meta.data$temp=="Norm"]="Normal prostate"
epi@meta.data$temp[epi@meta.data$temp=="ADT"]="post-ADT"
epi@meta.data$temp=factor(epi@meta.data$temp, levels=c("Normal prostate", "Indolent incidental PCa", "Aggressive incidental PCa", "Primary PCa", "mHSPC", "post-ADT", "CRPC", "NEPC"))
epi@meta.data$fig.pheno=epi@meta.data$temp
table(epi@meta.data$fig.pheno)
# Aggressive incidental PCa                      CRPC   Indolent incidental PCa
#                      2869                      6250                      3589
#                     mHSPC                      NEPC           Normal prostate
#                      7344                       169                     11945
#                  post-ADT               Primary PCa
#                       327                      6970
pdf("f.1d_phenotype_dimplot.pdf", width=7.5, height=5)
DimPlot(epi, label = FALSE, group.by="fig.pheno", pt.size=0.001, raster = FALSE)# + NoLegend()
dev.off()
epi@meta.data$temp=NA
epi@meta.data$temp[epi@meta.data$fig.pheno=="Aggressive incidental PCa"]="Aggressive incidental PCa"
png("dimplot.png", height=5*100, width=5*100)
DimPlot(epi, label = FALSE, group.by="temp")#, cols=c("blue", "red"))
dev.off()
### fe.1c fe.1d
pdf("fe.1c_cluster_dimplot.pdf", width=5.5, height=5)
DimPlot(epi, label = FALSE, group.by="fig.cluster", pt.size=0.01, raster = FALSE)# + NoLegend()
dev.off()
epi@meta.data$fig.sample.numb=as.numeric(as.factor(epi@meta.data$fig.sample))
pdf("fe.1d_sample_dimplot.pdf", width=7, height=5)
DimPlot(epi, label = FALSE, group.by="fig.sample.numb", pt.size=0.01, raster = FALSE)# + NoLegend()
dev.off()


### f1e
# dotmap, epithelial markers
imp.genes=c()
imp.genes=c( "STMN1", "HMMR", "PTTG1",   "CHGA", "SYP", "ENO2",   "AR", "FOXA1", "AMACR", "KLK3", "KRT8", "KRT18",    "KRT19", "TP63", "KRT5",   "KRT4", "TACSTD2", "PSCA")
# DES: myofibroblast
imp.genes=unique(imp.genes)
imp.genes=imp.genes[length(imp.genes):1]# reverse
# delete celltype with less than 100 cells
# epi=subset(epi, subset=celltype!="Myocyte")
epi@meta.data$fig.celltype=factor(epi@meta.data$fig.celltype, levels=c("Basal/typeC", "Luminal", "CRPC luminal", "NE", "CellCycle"))
pdf("f1_epithelial.dotmap_celltype_markers.pdf", width = 18, height = 3);
DotPlot(epi, features = imp.genes, cols = c('white', 'blue'), dot.scale = 6, group.by="fig.celltype");
dev.off();


### qusage to find activated pathways
require(qusage)
packageVersion("qusage")
# [1] '2.24.0'
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
epi.tumor@meta.data$for_temp_color=as.numeric(epi.tumor@meta.data$fig.celltype)
run_qusage_heatmap.seurat3(epi.tumor, nm = 'fig.celltype', kegg.list, my.seed=100)
run_qusage_heatmap.seurat3(epi.tumor, nm = 'fig.celltype', hall.list, my.seed=100)


### monocle 3
{
  ### cell state transition trajectory analysis
  unloadNamespace("monocle")
  require(monocle3)
  require(ggplot2)
  require(dplyr)
  # load data
  # expression_matrix, a numeric matrix of expression values, where rows are genes, and columns are cells
  # cell_metadata, a data frame, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
  gene_metadata=data.frame(gene_short_name=rownames(epi.tumor@assays$SCT@counts))
  rownames(gene_metadata)=rownames(epi.tumor@assays$SCT@counts)
  cds <- new_cell_data_set( epi.tumor@assays$SCT@counts, 
                            cell_metadata = epi.tumor@meta.data,
                            gene_metadata = gene_metadata)
  # pca 
  cds <- preprocess_cds(cds, num_dim = 100)
  plot_pc_variance_explained(cds)
  dev.off()
  # umap
  cds <- reduce_dimension(cds, cores=5)
  plot_cells(cds)
  dev.off()
  pdf("temp.pdf", width=5, height=5)
  plot_cells(cds, color_cells_by="spec.inte.comm", group_label_size = 7, cell_size=1)
  dev.off()
  pdf("temp.pdf", width=5, height=5)
  plot_cells(cds, color_cells_by="phenotype", group_label_size = 7, cell_size=1)
  dev.off()
  # clustering
  cds <- cluster_cells(cds)
  # plot monocle's clustering result
  pdf("temp.pdf", width=5, height=5)
  plot_cells(cds, color_cells_by = "partition", group_label_size = 7, cell_size=0.35)
  dev.off()
  # plot seurat's clustering result
  pdf("temp.pdf", width=5, height=5)
  plot_cells(cds, color_cells_by = "fig.cluster", group_label_size = 2, cell_size=0.35)
  dev.off()

  # as monocle's cluster result is not well, and as I like to keep consistent with seurat's result, replace PCA and UMAP in monocle with seurat's results.
  str( cds @ int_colData @ listData $ reducedDims @ listData $ PCA )
  str( epi.tumor @ reductions $ pca @ cell.embeddings )
  table( rownames( cds @ int_colData @ listData $ reducedDims @ listData $ PCA )==rownames( epi.tumor @ reductions $ pca @ cell.embeddings ) )
  cds @ int_colData @ listData $ reducedDims @ listData $ PCA = epi.tumor@reductions$harmony@cell.embeddings[, 1:25]
  str( cds @ int_colData @ listData $ reducedDims @ listData $ UMAP )
  str( epi.tumor @ reductions $ umap @ cell.embeddings )
  cds @ int_colData @ listData $ reducedDims @ listData $ UMAP = epi.tumor @ reductions $ umap @ cell.embeddings
  pdf("temp.pdf", width=5, height=5)
  plot_cells(cds, color_cells_by = "fig.cluster", group_label_size = 2, cell_size=0.35)
  dev.off()
  # learn trajectory
  # replace monocle result's "partition" with seurat result's "fig.cluster"
  str(cds @ clusters @ listData $ UMAP $ partitions)
   # Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
   # - attr(*, "names")= chr [1:3994] "AAACCCACAAGAGCTG-1_1" "AAACCCACAGTCTGGC-1_1" "AAACCCACATGGCTGC-1_1" "AAACCCAGTCTTGCTC-1_1" ...
  str(cds @ colData @ listData $ fig.cluster)
   # Factor w/ 14 levels "0","1","2","3",..: 3 2 12 12 1 2 2 13 2 1 ..
  table(names(cds @ clusters @ listData $ UMAP $ partitions)==rownames(cds @ colData) )
  # TRUE
  # 3994
  cds @ clusters @ listData $ UMAP $ partitions=cds @ colData @ listData $ fig.cluster
  names(cds @ clusters @ listData $ UMAP $ partitions) = rownames(cds @ colData)
  # idenfiy trajectory
  cds <- cluster_cells(cds, reduction_method="UMAP")
  cds <- learn_graph(cds, use_partition = FALSE, close_loop=FALSE, learn_graph_control=list(minimal_branch_len=6))
  plot_cells(cds,
             color_cells_by = "fig.celltype",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  dev.off()

  cds <- order_cells(cds, root_cells=colnames(cds)[cds@colData$fig.celltype=="Clike"])
  plot_cells(cds,
             color_cells_by = "pseudotime",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1.5,
             label_groups_by_cluster = FALSE,
             label_principal_points = FALSE)
  dev.off()

  # identify genes correlated with pseudotime

  # correlated with pseudotime?
  ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
  pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.01))
  plot_cells(cds, genes=c("TIE1", "CDH2", "KRT8"),
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
  # chr [1:13002] "AL669831.5" "AL645608.3" "AL645608.5" "AL645608.1" "SAMD11" ...
  # too many genes
}

### slingshot
{
  require(slingshot, quietly = FALSE)
  require(SingleCellExperiment)
  sce <- SingleCellExperiment(assays = List(counts = epi@assays$RNA@counts))
  rd1=epi@reductions$pca@cell.embeddings
  rd2=epi@reductions$umap@cell.embeddings
  rd3=epi@reductions$harmony@cell.embeddings[, 1:25]
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2, harmony=rd3)
  epitype=epi@meta.data$epitype
  names(epitype)=rownames(epi@meta.data)
  colData(sce)$epitype <- epitype
  sce <- slingshot(sce, clusterLabels = 'epitype', reducedDim = 'harmony')
  # not completed
}


### monocle 2
{
  # unloadNamespace("monocle3")
  require(monocle)

  # epi.tumor=readRDS("epi.tumor.highCNV.harmony.adjed.rds")
  temp=epi.tumor
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

  ### use all significant pair-wise epitype's DEGs as ordering genes
  {
    unicluster=unique(epi.tumor@meta.data$epitype)
    markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
    str(epi.tumor@assays$SCT@var.features)
    geneset=c()
    print("Including clusters:")
    print(unicluster)
    for(i in names(markerlist.epi))
    {
      j=markerlist.epi[[i]]
      mn=strsplit(i, split="__")[[1]]
      m=mn[1]
      n=mn[2]
      if( (m%in%unicluster)&(n%in%unicluster) )
      {
        j=j[ (j$p_val_adj<0.05) & ( (j$avg_log2FC)>log2(1.5) ), ]
        j=j[ order(j$avg_log2FC, decreasing=TRUE), ]
        # if(nrow(j)>5)
        # {j=j[1:5, ]}
        newset=rownames(j)
        geneset=union( geneset, newset )
        #  )
      }
    }
    print(str(geneset))
    # chr [1:4367] "GATA3" "PIGR" "KRT7" "CX3CL1" "LCN2" "CLIC6" "UPK1A" "ALOX5" ... p_val_adj<0.05 & abs(log2(FC))>log2(1.5)
    # chr [1:155] "LCN2" "UPK1A" "PSCA" "CD74" "SPINK1" "DHRS2" "S100A9" "HPGD" ... p_val_adj<0.05 & FC>1.5 & top5
    # chr [1:269] "LCN2" "UPK1A" "PSCA" "CD74" "SPINK1" "S100A9" "DHRS2" "HLA-A" ... p_val_adj<0.05 & FC>1.5 & top10
    # chr [1:4088] "LCN2" "UPK1A" "PSCA" "CD74" "SPINK1" "S100A9" "DHRS2" ...p_val_adj<0.05 & abs(log2(FC))>log2(1.5) -----BEST
    ordering.genes=unique(geneset)
  }

  # ### use epi.tumor's fig.cluster's markers as ordering genes
  # {
  #   tumor.markers =readRDS("epi.tumor.markers.as.celltype.byHarmony.rds")
  #   markers=tumor.markers
  #   markers=markers[markers$p_val_adj<0.05, ]
  #   markers=markers[abs(markers$avg_log2FC)>log2(1.5), ]
  #   markers$foldChange=2^(markers$avg_log2FC)
  #   # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
  #   markers=markers[order(markers$foldChange, decreasing=TRUE), ]
  #   ordering.genes=unique(markers$gene)
  #   # ordering.genes=unique(rownames(markers))
  #   str(ordering.genes)
  #   # chr [1:1781] "SCGB1A1" "LCN2" "LTF" "BPIFB1" "S100P" "PSCA" "OLFM4" ...
  # }

 
  temp.m <- setOrderingFilter(temp.m,  ordering.genes) # 
  pdf("plot_ordering_genes.pdf")
  plot_ordering_genes(temp.m)
  dev.off()

  # temp.m.ica <- reduceDimension(temp.m, max_components = 2, method = 'ICA'
      # ) # method = 'DDRTree' / 'ICA'
  temp.m.tree <- reduceDimension(temp.m, max_components = 2, method = 'DDRTree'
      ) # method = 'DDRTree' / 'ICA'
  temp.m=temp.m.tree
  # saveRDS( temp.m.ica, "temp.m.ica" )
  # saveRDS( temp.m.tree, "temp.m.tree" )
  # temp.m=readRDS( "temp.m.ica" )
  # temp.m=readRDS( "temp.m.tree" )

  temp.m <- orderCells(temp.m)
  pdf("temp_trajectory_epitype.pdf")
  temp=plot_cell_trajectory(temp.m, color_by = "State")
  print(temp)
  dev.off()
  pdf("temp_trajectory_celltype.pdf")
  temp=plot_cell_trajectory(temp.m, color_by = "celltype")
  print(temp)
  dev.off()

  GM_state <- function(temp.m){
    if (length(unique(pData(temp.m)$State)) > 1){
      T0_counts <- table(pData(temp.m)$State, pData(temp.m)$fig.celltype)[,"Clike"]
      return(as.numeric(names(T0_counts)[which
            (T0_counts == max(T0_counts))]))
    } else {
      return (1)
    }
  }
  table(GM_state(temp.m))

  # set a indicator for time, and sort again
  str(pData(temp.m)$State)
  table(pData(temp.m)$State)
  # pData(temp.m)$State=pData(temp.m)$fig.celltype
  temp.m <- orderCells(temp.m, root_state = GM_state(temp.m))
  saveRDS(temp.m, "epi.tumor.monocle.rds")
  # temp.m=readRDS("epi.tumor.monocle.rds")

  ### f2g
  pdf("f.2gs_plot_cell_trajectory_byState_temp.m.pdf") # , width=4, height=5
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "State", cell_size=0.5)
  print(temp)
  dev.off()

  png("plot_cell_trajectory_byPseudotime_temp.m.png") # , width=4, height=5
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "Pseudotime", cell_size=0.5)
  print(temp)
  dev.off()

  pdf("f.2gs_plot_cell_trajectory_byPseudotime_temp.m.pdf") # , width=4, height=5
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "Pseudotime", cell_size=0.5)
  print(temp)
  dev.off()

  png("f.2g_plot_cell_trajectory_fig.celltype_temp.m.png") # , width=4, height=5
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "fig.celltype", cell_size=0.5)
  print(temp)
  dev.off()

  pdf("f.2g_plot_cell_trajectory_fig.celltype_temp.m.pdf") # , width=4, height=5
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "factor.fig.celltype", cell_size=0.5)
  print(temp)
  dev.off()

  png("plot_cell_trajectory_fig.pheno_temp.m.png") # , width=5, height=5
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "fig.pheno", cell_size=0.2, cell_name_size = 8)
  print(temp)
  dev.off()

  # str(pData(temp.m)$fig.pheno)
  # temp.levels=levels(pData(temp.m)$fig.pheno)
  # # delete "Normal prostate"
  # temp.levels=temp.levels[2:length(temp.levels)]
  # pData(temp.m)$fig.pheno=factor(pData(temp.m)$fig.pheno, levels=temp.levels)
  pdf("plot_cell_trajectory_fig.pheno_temp.m.pdf", width=5, height=6)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "fig.pheno", cell_size=0.2, cell_name_size = 8)
  print(temp)
  dev.off()

  pdf("plot_cell_trajectory_fig.zone_temp.m.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_temp.m.png")
  temp=plot_cell_trajectory(temp.m, color_by = "fig.zone", cell_size=0.05, cell_name_size = 8)
  print(temp)
  dev.off()

  temp.m.time=temp.m@phenoData@data$Pseudotime
  names(temp.m.time)=rownames(temp.m@phenoData@data)
  table(rownames(temp.m@phenoData@data)==rownames(epi.tumor@meta.data))
  # TRUE
  # 5255
  temp.m@phenoData@data$fig.celltype=epi.tumor@meta.data$celltype
  table(temp.m@phenoData@data$fig.celltype)
  oorder=c("Clike", "TypeC",  "Basal" , "Intermediate", "Luminal", "CellCycle", "NE")
  temp.m@phenoData@data$factor.fig.celltype=factor(temp.m@phenoData@data$fig.celltype, levels=oorder)
  require("scales")
  color=hue_pal()(length(oorder))
  color[length(oorder)]="#FF6600"


  png("plot_cell_trajectory_details_fig.celltype_temp.m.png", width=30*100, height=6*100)
  temp=plot_cell_trajectory(temp.m, color_by = "fig.celltype") +
      facet_wrap(~fig.celltype, nrow = 1)
  print(temp)
  dev.off()

  pdf("plot_cell_trajectory_details_fig.celltype_temp.m.pdf", width=40, height=6)
  temp=plot_cell_trajectory(temp.m, color_by = "factor.fig.celltype", cell_size=0.2) +
      facet_wrap(~fig.celltype, nrow = 1)
  print(temp)
  dev.off()

  png("f.2h_plot_cell_trajectory_details_fig.pheno_temp.m.png", width=22*100, height=6*100)
  temp=plot_cell_trajectory(temp.m, color_by = "fig.pheno") +
      facet_wrap(~fig.pheno, nrow = 1)
  print(temp)
  dev.off()

  pdf("f.2h_plot_cell_trajectory_details_fig.pheno_temp.m.pdf", width=22, height=6)
  temp=plot_cell_trajectory(temp.m, color_by = "fig.pheno", cell_size=0.2) +
      facet_wrap(~fig.pheno, nrow = 1)
  print(temp)
  dev.off()


  ### barplot the percentage of each phenotype in states
  {
    # give temp.m's State to epi.tumor
    temp.m.state=temp.m@phenoData@data$State
    names(temp.m.state)=rownames(temp.m@phenoData@data)
    epi.tumor@meta.data[names(temp.m.state), "State"]=temp.m.state
    epi.tumor@meta.data$State=as.character(epi.tumor@meta.data$State)

    # group monocle's state
    table(epi.tumor@meta.data$State)
    current.cluster.ids <- sort(unique(epi.tumor@meta.data$State))
    print(current.cluster.ids)
    # [1] "1" "2" "3" "4" "5" "6" "7"
    new.cluster.ids <- c("Starting state (state 1)", "State associated with reccurence of CRPC/NEPC (state 2, 3 and 4)", "State associated with reccurence of CRPC/NEPC (state 2, 3 and 4)", "State associated with reccurence of CRPC/NEPC (state 2, 3 and 4)", "State responsive to ADT (state 5, 6 and 7", "State responsive to ADT (state 5, 6 and 7", "State responsive to ADT (state 5, 6 and 7")
    epi.tumor@meta.data$State3 <- plyr::mapvalues(x = epi.tumor@meta.data$State, from = current.cluster.ids, to = new.cluster.ids)
    table(epi.tumor@meta.data$State3)
    str(epi.tumor@meta.data$State3)

    # barplot phenotype vs monocle's predicted states 
    state.sta=table(epi.tumor@meta.data$fig.pheno, epi.tumor@meta.data$State3)
    state.sta=apply(state.sta, 1, function(x) 100*x/sum(x))
    state.sta=t(state.sta)
    # ### delete "Normal prostate"
    # str(epi.tumor@meta.data$fig.pheno)
    # temp.levels=levels(epi.tumor@meta.data$fig.pheno)
    # temp.levels=temp.levels[2:length(temp.levels)]
    # epi.tumor@meta.data$fig.pheno=factor(epi.tumor@meta.data$fig.pheno, levels=temp.levels)
    barN=ncol(state.sta)
    n_top = nrow(state.sta)
    pdf(paste("temp.barplot.statesInSample.pdf",sep=""), width=10, height=3)
    # the width between paper edge and plot ( down,left,up,right : four directions)
    par(mar = c(7, 4, 1, 12), xpd=T)
    # color=DiscretePalette(n_top, palette = "polychrome")
    require("scales")
    color=hue_pal()(n_top)
    barplot( state.sta, xlab="", ylab="Percentage", col=color, names.arg=colnames(state.sta), las=1, cex.names=0.5, beside=TRUE) # , xaxt="n"
    # facts=barN+1
    end_point = barN*n_top #0.5 + n_top *facts-1
    # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(state.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
    legend( end_point+1,100, rownames(state.sta), fill=color, bty='n', cex = 0.4 );
    # , bty='n': no-frame
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = state.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(state.sta[1,],1),"%") , cex=0.7 )
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-state.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(state.sta[2,],1),"%") , cex=0.7 )
    dev.off()
  }

  # density plot of cell types' in pseudotime
  {
    # refer to www.r-graph-gallery.com/135-stacked-density-graph.html
    # refer to www.sthda.com/english/wiki/ggplot2-density-plot-quick-start-guide-r-software-and-data-visualization
    library(ggplot2)
    # library(hrbrthemes)
    library(dplyr)
    library(tidyr)
    library(viridis)

    # give monocle's pseudotime to seurat object
    temp.m.time=temp.m@phenoData@data$Pseudotime
    names(temp.m.time)=rownames(temp.m@phenoData@data)
    epi.tumor@meta.data[names(temp.m.time), "pseudotime"]=temp.m.time

    temp.df=epi.tumor@meta.data[, c("factor.fig.celltype", "pseudotime")]
    # get colors
    require("scales")
    uni.color=levels(epi.tumor@meta.data$factor.fig.celltype)
    color=hue_pal()(length(uni.color))
    # color=color[length(oorder):1]# 
    # specify Clike's colour
    # color[length(oorder)]="#FF6600"
    pdf("f.2h.CellTypeDensityAlongPseudotime.pdf", width=10, height=5)
    # transparancy desnity plot. alpha: transparent degree
    # {
    #   p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=factor.fig.celltype, fill=factor.fig.celltype)) +
    #     geom_density(adjust=1.5, alpha=.4)
    # }
    # stacked density plot
    # {
      # p1 <- ggplot(data=temp.df, aes(x=pseudotime, group=factor.fig.celltype, fill=factor.fig.celltype)) +
      # geom_density(adjust=1.5, position="fill")
    # }
    {
      p1<- ggplot(temp.df, aes(x=pseudotime, color=factor.fig.celltype)) +
      geom_density() +
      labs(title="Cell type density curve",x="Pseudotime", y = "Density") +
      scale_color_manual(values=color)+
      theme_classic()
    }
    print(p1)
    dev.off()
  }

  ### boxplot Clike percentage for CRPC+NEPC vs mHSPC vs Primary and incidental aggressive vs indolent
  {
    # give sample ID as phenotype + number
    {
      pheno.sample=epi.tumor@meta.data[, c("fig.pheno", "fig.sample")]
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
      epi.tumor@meta.data$pheno.sample=pheno.sample[rownames(epi.tumor@meta.data), ]$pheno.sample
    }

    # get percentage table of phenotype+sampleID vs celltype
    state.sta=table(epi.tumor@meta.data$pheno.sample, epi.tumor@meta.data$fig.celltype)
    state.sta=apply(state.sta, 1, function(x) 100*x/sum(x))
    state.sta=t(state.sta)
    # map the colnames of state.sta from pheno.sample into phenotype 
    {
      # map table state.sta's pheno.sample to phenotype 
      temp.map=epi.tumor@meta.data[, c("pheno.sample", "fig.pheno")]
      table(duplicated(temp.map))
      temp.map=temp.map[!duplicated(temp.map), ]
      current.cluster.ids <-  as.character(temp.map$pheno.sample)
      print(current.cluster.ids)
      new.cluster.ids <- as.character(temp.map$fig.pheno)
      rownames(state.sta) <- plyr::mapvalues(x = rownames(state.sta), from = current.cluster.ids, to = new.cluster.ids)
      temp=levels(temp.map$fig.pheno)
      print(temp)
      # temp=temp[2:length(temp)]# delete "Normal prostate"
      rownames(state.sta)[rownames(state.sta)%in%c("CRPC", "NEPC")]="CRPC/NEPC"
    }
    ### only consider Clike's percentage 
    temp.df=data.frame(pheno=rownames(state.sta), pct=state.sta[, "Clike"])
    # only consider primary, mHSPC, CRPC/NEPC
    df=temp.df[temp.df$pheno%in%c("Primary PCa", "mHSPC", "CRPC/NEPC"), ]
    df$pheno=factor(df$pheno, levels=c("Primary PCa", "mHSPC", "CRPC/NEPC"))
    require(ggplot2)
    pdf("f.2j.boxplot.ClikePCTacrossPheno.pdf", width=5, height=5)
    p <- ggplot(df, aes(x=pheno, y=pct)) + geom_boxplot() + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    theme_classic() 
    # geom_jitter(shape=16, position=position_jitter(0.2))
    print(p)
    # boxplot( pct~pheno, data=df, ylab="Percentage of Clike cells" )
    #  main=paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test"
    res <- wilcox.test(df$pct[df$pheno=="CRPC/NEPC"], df$pct[df$pheno=="mHSPC"])
    res <- wilcox.test(df$pct[df$pheno=="mHSPC"], df$pct[df$pheno=="Primary PCa"])
    # legend( "topright", paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test", sep="") ,cex=0.8);
    dev.off()

    # only consider primary, mHSPC, CRPC/NEPC
    df=temp.df[temp.df$pheno%in%c("Indolent incidental PCa", "Aggressive incidental PCa"), ]
    df$pheno=factor(df$pheno, levels=c("Indolent incidental PCa", "Aggressive incidental PCa"))
    require(ggplot2)
    pdf("f.2k.boxplot.ClikePCTacrossPheno.pdf", width=4, height=5)
    p <- ggplot(df, aes(x=pheno, y=pct)) + geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    theme_classic()
    # geom_jitter(shape=16, position=position_jitter(0.2))
    print(p)
    # boxplot( pct~pheno, data=df, ylab="Percentage of Clike cells" )
    #  main=paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test"
    res <- t.test(df$pct[df$pheno=="Aggressive incidental PCa"], df$pct[df$pheno=="Indolent incidental PCa"])
    # legend( "topright", paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test", sep="") ,cex=0.8);
    dev.off()

    # barplot
    df=temp.df
    df$pheno=factor(df$pheno, levels=c("Indolent incidental PCa", "Aggressive incidental PCa", "Primary PCa", "mHSPC", "CRPC/NEPC"))
    require(ggplot2)
    pdf("f.2j.boxplot.ClikePCTacrossPheno.pdf", width=5, height=5)
    p<-ggplot(data=df, aes(x=pheno, y=pct)) +
            geom_bar(stat="identity")

    print(p)
    # boxplot( pct~pheno, data=df, ylab="Percentage of Clike cells" )
    #  main=paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test"
    res <- t.test(df$pct[df$pheno=="Aggressive incidental PCa"], df$pct[df$pheno=="Indolent incidental PCa"])
    # legend( "topright", paste("p-value: ", formatC(res$p.value, digits=3), " Wilcoxon rank sum test", sep="") ,cex=0.8);
    dev.off()
  }

  ### barplot the percentage of each phenotype+sampleID in states
  {
    # give temp.m's State to epi.tumor
    temp.m.state=temp.m@phenoData@data$State
    names(temp.m.state)=rownames(temp.m@phenoData@data)
    epi.tumor@meta.data[names(temp.m.state), "State"]=temp.m.state
    epi.tumor@meta.data$State=as.character(epi.tumor@meta.data$State)

    # give sample ID as phenotype + number
    {
      pheno.sample=epi.tumor@meta.data[, c("fig.pheno", "fig.sample")]
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
      epi.tumor@meta.data$pheno.sample=pheno.sample[rownames(epi.tumor@meta.data), ]$pheno.sample
    }

    # group monocle's state
    table(epi.tumor@meta.data$State)
    current.cluster.ids <- sort(unique(epi.tumor@meta.data$State))
    print(current.cluster.ids)
    # [1] "1" "2" "3" "4" "5" "6" "7"
    new.cluster.ids <- c("Starting state (state 1)", "State associated with reccurence of CRPC/NEPC (state 2, 3 and 4)", "State associated with reccurence of CRPC/NEPC (state 2, 3 and 4)", "State associated with reccurence of CRPC/NEPC (state 2, 3 and 4)", "State responsive to ADT (state 5, 6 and 7", "State responsive to ADT (state 5, 6 and 7", "State responsive to ADT (state 5, 6 and 7")
    epi.tumor@meta.data$State3 <- plyr::mapvalues(x = epi.tumor@meta.data$State, from = current.cluster.ids, to = new.cluster.ids)
    table(epi.tumor@meta.data$State3)
    str(epi.tumor@meta.data$State3)

    # barplot phenotype vs monocle's predicted states 
    # state.sta=table(epi.tumor@meta.data$fig.pheno, epi.tumor@meta.data$State3)
    state.sta=table(epi.tumor@meta.data$pheno.sample, epi.tumor@meta.data$State3)
    state.sta=apply(state.sta, 1, function(x) 100*x/sum(x))
    state.sta=t(state.sta)
    # map the colnames of state.sta from pheno.sample into phenotype 
    {
      # map table state.sta's pheno.sample to phenotype 
      temp.map=epi.tumor@meta.data[, c("pheno.sample", "fig.pheno")]
      table(duplicated(temp.map))
      temp.map=temp.map[!duplicated(temp.map), ]
      current.cluster.ids <-  as.character(temp.map$pheno.sample)
      print(current.cluster.ids)
      new.cluster.ids <- as.character(temp.map$fig.pheno)
      rownames(state.sta) <- plyr::mapvalues(x = rownames(state.sta), from = current.cluster.ids, to = new.cluster.ids)
      temp=levels(temp.map$fig.pheno)
      print(temp)
      temp=temp[2:length(temp)]# delete "Normal prostate"
    }
    barN=ncol(state.sta)
    n_top = nrow(state.sta)
    n_pheno=length(unique(rownames(state.sta)))
    pheno.count=table(rownames(state.sta))
    # order the rows as levels of fig.pheno
    pheno.count=pheno.count[temp]
    pdf(paste("temp.barplot.statesInSample.pdf",sep=""), width=10, height=3)
    # the width between paper edge and plot ( down,left,up,right : four directions)
    par(mar = c(7, 4, 1, 12), xpd=T)
    # color=DiscretePalette(n_top, palette = "polychrome")
    require("scales")
    temp.col=hue_pal()(n_pheno)
    color=c()
    for(i in 1:n_pheno)
    {
      temp.name=names(pheno.count)[i]
      color=c(color, rep(temp.col[i], pheno.count[temp.name]))
    }
    barplot( state.sta, xlab="", ylab="Percentage", col=color, names.arg=colnames(state.sta), las=1, cex.names=0.5, beside=TRUE) # , xaxt="n"
    # facts=barN+1
    end_point = barN*n_top #0.5 + n_top *facts-1
    # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(state.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
    legend( end_point+1,100, rownames(state.sta), fill=color, bty='n', cex = 0.4 );
    # , bty='n': no-frame
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = state.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(state.sta[1,],1),"%") , cex=0.7 )
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-state.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(state.sta[2,],1),"%") , cex=0.7 )
    dev.off()
  }


  ### barplot the percentage of each phenotype in cell types
  {
    # give sample ID as phenotype + number
    {
      pheno.sample=epi.tumor@meta.data[, c("fig.pheno", "fig.sample")]
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
      epi.tumor@meta.data$pheno.sample=pheno.sample[rownames(epi.tumor@meta.data), ]$pheno.sample
    }

    # get the percentage table
    state.sta=table(epi.tumor@meta.data$pheno.sample, epi.tumor@meta.data$fig.celltype)
    state.sta=apply(state.sta, 1, function(x) 100*x/sum(x))
    # state.sta=t(state.sta)
    # state.sta=state.sta[, c("Clike", "Neuroendocrine")]
    # map the colnames of state.sta from pheno.sample into phenotype 
    {
      # map table state.sta's pheno.sample to phenotype 
      temp.map=epi.tumor@meta.data[, c("pheno.sample", "fig.pheno")]
      table(duplicated(temp.map))
      temp.map=temp.map[!duplicated(temp.map), ]
      current.cluster.ids <-  as.character(temp.map$pheno.sample)
      print(current.cluster.ids)
      new.cluster.ids <- as.character(temp.map$fig.pheno)
      # rownames(state.sta) <- plyr::mapvalues(x = rownames(state.sta), from = current.cluster.ids, to = new.cluster.ids)
      colnames(state.sta) <- plyr::mapvalues(x = colnames(state.sta), from = current.cluster.ids, to = new.cluster.ids)
      temp=levels(temp.map$fig.pheno)
      print(temp)
      temp=temp[2:length(temp)]# delete "Normal prostate"
    }
    barN=ncol(state.sta)
    n_top = nrow(state.sta)
    # n_pheno=length(unique(rownames(state.sta)))
    # pheno.count=table(rownames(state.sta))
    n_pheno=length(unique(colnames(state.sta)))
    pheno.count=table(colnames(state.sta))
    # order the rows as levels of fig.pheno
    pheno.count=pheno.count[temp]
    pdf(paste("temp.barplot.CelltypeInSample.pdf",sep=""), width=10, height=3)
    # the width between paper edge and plot ( down,left,up,right : four directions)
    par(mar = c(6, 4, 1, 12), xpd=T)
    # color=DiscretePalette(n_top, palette = "polychrome")
    require("scales")
    color=hue_pal()(n_top)
    barplot( state.sta, xlab="", ylab="Percentage", col=color, names.arg=colnames(state.sta), las=2, cex.names=0.5, beside=FALSE) # , xaxt="n"
    # las=2 : vertical x-names. las=1 : horizontal x-names 
    # facts=barN+1
    end_point = barN+8 #0.5 + n_top *facts-1
    # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(state.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
    legend( end_point,100, rownames(state.sta), fill=color, bty='n', cex = 0.8 );
    # , bty='n': no-frame
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = state.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(state.sta[1,],1),"%") , cex=0.7 )
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-state.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(state.sta[2,],1),"%") , cex=0.7 )
    dev.off()
  }


  png("plot_cell_trajectory_details_epitype_temp.m.png", width=60*100, height=6*100)
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


### dimplot as Gleason
table(epi@meta.data$fig.patient)
current.cluster.ids <- sort(unique(epi@meta.data$fig.patient))
print(current.cluster.ids)
#  [1] "53"  "54"  "55"  "56"  "59"  "62"  "71"  "73"  "74"  "75"  "76"  "77"
# [13] "CBH" "GHL" "GJC" "HYL" "QJZ" "SQB" "XYM" "ZYM"
new.cluster.ids <- c("9", "10", "9", "9", "7", "7", "9", "9", "3+4", "9", "7", "9", "NORM", "NORM", "IND", "CRPC", "NEPC", "AGG", "CRPC", "IND")
epi@meta.data$gleason <- plyr::mapvalues(x = as.character(epi@meta.data$fig.patient), from = current.cluster.ids, to = new.cluster.ids)
table(epi@meta.data$gleason)
#    0   10  3+4    6    7    9 CRPC NEPC
# 3618   64  243 1994 3545 2219 2196  560
png("dimplot.gleason.png", height=5*140, width=5*140)
DimPlot(epi, label = FALSE, group.by="gleason", pt.size=1, cols=c("#CCFFCC", "#00FF00", "#66FFFF", "#0000FF", "#FF33FF", "#FF0000", "BLACK", "#CC9933", "#FFFF33"), order=c("NEPC", "CRPC", "10", "9", "7", "3+4", "AGG", "IND", "0")
)
dev.off()

# check known epithelial cell types
check.genes=c("KLK3", "KRT14", "KRT13", "SCGB1A1", "KRT19", "KRT5", "CHGA", "SYP", "ENO2", "TP63", "AR", "nFeature_RNA", "nFeature_SCT", "KRT4", "TACSTD2", "PSCA")
# # prostate/urinary-cancer-specific marker
# check.genes=c("GATA3", "MMP9", "SLC45A3", "NKX3-1" )#"UPK2"
# # tumor-initiating markers
# check.genes=c("NKX3-1", "LY6D")
# # zhangbo's markers
# check.genes=c("COL6A1", "PCA3", "SPON2", "TDRD1")

DefaultAssay(epi)="RNA"
png("temp.featureplot.png", height=10*100, width=14*100)
FeaturePlot(epi, features = check.genes)
dev.off()
DefaultAssay(epi)="SCT"
png("temp.VlnPlot.png.png", width=10*100, height=8*100)
epi=SetIdent(epi, value="epitype")
temp=VlnPlot(epi, features = check.genes, pt.size = 0, group.by=NULL, ncol = 2)
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()

# epi.bak=epi
epi@assays$RNA@scale.data=matrix(0, 0, 0)
epi@assays$SCT.bak=NULL
epi@assays$SCT.bak.2=NULL
epi@commands=list()
table(epi@meta.data$epitype)
# saveRDS(epi, "epi.downsampled.rds")
# epi=readRDS("epi.downsampled.rds")


### get closest distance between inter-phenotype and intra-phenotype clusters
markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
# epi.indagg=readRDS("epi.indagg.seurat3.indivSCTed.celltyped.rds")
# epi.crpc=readRDS("epi.crpc.seurat3.indivSCTed.celltyped.rds")
# normale=readRDS("normale.seurat3.indivSCTed.celltyped.rds")
# epi13=readRDS("epi13.seurat3.indivSCTed.celltyped.rds")


# greedy search shortest evolution paths via combing:

### check distance between specified two layers
dis(markerlist=markerlist.epi, "tL1", "iL0")

### consider all 4 levels: normal p, early agg/indo PCa, primary PCa, and CRPC/NEPC
lv.1=sort(unique(normale@meta.data$epitype))
lv.2=sort(unique(incie@meta.data$epitype))
prime=subset(prime, subset=epitype!="tC24")
lv.3=sort(unique(prime@meta.data$epitype))
lv.4=sort(unique(crpce@meta.data$epitype))
lv=list()
for(i in 1:4)
{
  lv[[i]]=get(paste("lv", i, sep="."))
}

shortest.4=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=4)
# print(shortest.4)
print(easy.shortest(shortest.4))

shortest.3=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=3)
# print(shortest.3)
print(easy.shortest(shortest.3))

shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
# print(shortest.2)
print(easy.shortest(shortest.2))


### excluding early tumor samples, only 3 levels: normal p, primary PCa, CRPC/NEPC
lv.1=sort(unique(normale@meta.data$epitype))
# lv.1=sort(unique(epi.indagg@meta.data$epitype))
lv.2=sort(unique(epi13@meta.data$epitype))
lv.3=sort(unique(epi.crpc@meta.data$epitype))
lv=list()
for(i in 1:3)
{
  lv[[i]]=get(paste("lv", i, sep="."))
}

shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
print(shortest)
print(easy.shortest(shortest.2))

shortest.3=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=3)
# print(shortest)
print(easy.shortest(shortest.3))



### only plot UMAP for specified layers to show evolution paths and phenotype more clear
### show evolution of normal p --> primary PCa
n2t=subset(epi, subset=(pheno%in%c("Norm", "Prim")))
png("n2t_celltype_dimplot.png", width=5*150, height=5*150)
DimPlot(n2t, label = FALSE, group.by="epitype", pt.size=1)# + NoLegend()
dev.off()
png("n2t_gleason_dimplot.png", height=5*150, width=5*150)
DimPlot(n2t, label = FALSE, group.by="gleason", pt.size=1, cols=c("#CCFFCC", "#0000FF", "#FF33FF", "#FF0000", "BLACK"), order=c("10", "9", "7", "3+4", "0")
)
dev.off()


### show evolution of normal --> CRPC/NEPC
n2c=subset(epi, subset=(pheno%in%c("Norm", "CRPC", "NEPC")))
png("n2c_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(n2c, label = TRUE, group.by="epitype", pt.size=0.01)# + NoLegend()
dev.off()
png("n2c_gleason_dimplot.png", height=5*140, width=5*140)
DimPlot(n2c, label = FALSE, group.by="gleason", pt.size=1, cols=c("#CCFFCC", "#CC9933", "#FFFF33"), order=c("NEPC", "CRPC", "NORM")
)
dev.off()


### show evolution of within early indolent/aggressive PCa
n2et=subset(epi, subset=(pheno%in%c("Indo", "Aggr")))
png("n2et_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(n2et, label = TRUE, group.by="epitype", pt.size=0.01)# + NoLegend()
dev.off()
png("n2et_gleason_dimplot.png", height=5*140, width=5*140)
DimPlot(n2et, label = FALSE, group.by="gleason", pt.size=1, cols=c("#00FF00", "#66FFFF"), order=c("Aggr", "Indo"))
dev.off()


### consider all 4 levels: normal p, early agg/indo PCa, primary PCa, and CRPC/NEPC
lv.1=c("etC1")
lv.2=sort(unique(epi.indagg@meta.data$epitype))
lv=list()
for(i in 1:2)
{
  lv[[i]]=get(paste("lv", i, sep="."))
}

shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
# print(shortest.2)
print(easy.shortest(shortest.2, show.second.shortest=F))


lv.1=sort(unique(epi.indagg@meta.data$epitype))
lv.2=sort(unique(epi.indagg@meta.data$epitype))
lv=list()
for(i in 1:2)
{
  lv[[i]]=get(paste("lv", i, sep="."))
}

shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
# print(shortest.2)
print(easy.shortest(shortest.2, show.second.shortest=TRUE))


lv.1=sort(unique(prime@meta.data$epitype))
lv.2=sort(unique(crpce@meta.data$epitype))
lv=list()
for(i in 1:2)
{
  lv[[i]]=get(paste("lv", i, sep="."))
}

shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
# print(shortest.2)
print(easy.shortest(shortest.2, show.second.shortest=FALSE))

### show evolution from normal p to early indolent/aggressive PCa
n2et=subset(epi, subset=(pheno%in%c("Norm", "Indo", "Aggr")))
png("n2et_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(n2et, label = TRUE, group.by="epitype", pt.size=0.01)# + NoLegend()
dev.off()
png("n2et_gleason_dimplot.png", height=5*140, width=5*140)
DimPlot(n2et, label = FALSE, group.by="gleason", pt.size=1, cols=c("#CCFFCC", "#00FF00", "#66FFFF"), order=c("AGG", "IND", "NORM"))
dev.off()


lv.1=sort(unique(epi.indagg@meta.data$epitype))
lv.2=sort(unique(epi13@meta.data$epitype))
lv=list()
for(i in 1:2)
{
  lv[[i]]=get(paste("lv", i, sep="."))
}

shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
# print(shortest.2)
print(easy.shortest(shortest.2, show.second.shortest=FALSE))

### show evolution from normal p to early indolent/aggressive PCa
et2t=subset(epi, subset=(pheno%in%c("Prim", "Indo", "Aggr")))
png("n2et_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(et2t, label = TRUE, group.by="epitype", pt.size=0.01)# + NoLegend()
dev.off()
png("n2et_gleason_dimplot.png", height=5*140, width=5*140)
DimPlot(et2t, label = FALSE, group.by="gleason", pt.size=1, cols=c("#00FF00", "#66FFFF", "#0000FF", "#FF33FF", "#FF0000", "BLACK"), order=c("10", "9", "7", "3+4", "AGG", "IND"))
dev.off()


# split CNV-high clusters and CNV-low clusters in early tumor and predict evolution paths
lv.1=sort(unique(epi.indagg@meta.data$epitype))
lv.2=sort(unique(epi13@meta.data$epitype))
lv=list()
for(i in 1:2)
{
  lv[[i]]=get(paste("lv", i, sep="."))
}

shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
# print(shortest.2)
print(easy.shortest(shortest.2, show.second.shortest=FALSE))

### show evolution from normal p to early indolent/aggressive PCa
et2t=subset(epi, subset=(pheno%in%c("Prim", "Indo", "Aggr")))
png("n2et_celltype_dimplot.png", width=5*100, height=5*100)
DimPlot(et2t, label = TRUE, group.by="epitype", pt.size=0.01)# + NoLegend()
dev.off()
png("n2et_gleason_dimplot.png", height=5*140, width=5*140)
DimPlot(et2t, label = FALSE, group.by="gleason", pt.size=1, cols=c("#00FF00", "#66FFFF", "#0000FF", "#FF33FF", "#FF0000", "BLACK"), order=c("10", "9", "7", "3+4", "AGG", "IND"))
dev.off()


### deconvlution to validate association between percentage of specific celltype and clinical features, e.g. Gleason, T stage, N stage, M stage, and BCR/OS.
{
  ### deconvlution's signature preparation
  # epi13=readRDS("epi13.seurat3.indivSCTed.celltyped.rds")
  ### mark epi13's cells as tumor-initiating celltypes
  table(epi13@meta.data$epitype)
    # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
    # 213   165   109   115   219    87   632 11931  6763  4929   446
  current.cluster.ids <- c("tB1", "tB2", "tB4", "tBC3", "tBI0", "tBI5", "tCyc", "tL0", "tL1", "tL2", "tLI3")
  new.cluster.ids <- c("TypeC", "TypeC", "Basal", "TypeC", "Basal", "Basal", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal") # c("Peripheral", "Transition", "Central")
  epi13@meta.data$ori <- plyr::mapvalues(x = as.character(epi13@meta.data$epitype), from = current.cluster.ids, to = new.cluster.ids)
  epi13@meta.data$ori=as.factor(epi13@meta.data$ori)
  epi13@meta.data$ori=as.character(epi13@meta.data$ori)
  table(epi13@meta.data$ori)
  # saveRDS(epi13, "epi13.seurat3.indivSCTed.celltyped.rds")
    # Basal Luminal   TypeC
    #   415   24701     493
  # consider celltypes in all cell types
  prim13=readRDS("prim13.seurat3.indivSCTed.celltyped.rds")
  table(prim13@meta.data$celltype)
      #   Basal cell Endothelial cell  Epithelial cell       Fibroblast
      #          964             3212              651              408
      # Luminal cell       Macrophage        Mast cell    Myofibroblast
      #        23994             2142             1443             1798
      #       T cell
      #         3135
  prim13@meta.data$ori=prim13@meta.data$celltype
  epi13.ori=epi13@meta.data$ori
  names(epi13.ori)=rownames(epi13@meta.data)
  table(epi13.ori)
    # Basal Luminal   TypeC
    #   415   24701     493
  prim13@meta.data[names(epi13.ori), "ori"]=epi13.ori
  table(prim13@meta.data$ori)
        #      Basal Endothelial cell       Fibroblast          Luminal
        #        415             3212              408            24701
        # Macrophage        Mast cell    Myofibroblast           T cell
        #       2142             1443             1798             3135
        #      TypeC
        #        493
  # prim13.sct.data.bak=prim13@assays$SCT@data
  # prim13@assays$SCT@data=prim13@assays$SCT@scale.data
  prim13=SetIdent(prim13, value = "ori")
  DefaultAssay(prim13)="SCT"
  # saveRDS(prim13, "prim13.seurat3.indivSCTed.celltyped.rds")
  # pairwise comparison
  celltypes=sort(unique(prim13@meta.data$ori))
  markerlist=list()
  for(i in celltypes)
  {
    for(j in celltypes)
    {
      if(i!=j)
      {
        tempname=paste(i, j, sep="_")
        print("tempname:")
        print(tempname)
        markerlist[[tempname]]=FindMarkers(prim13, min.pct = 0, ident.1 = i, ident.2=j, assay="SCT")
      }
    }
  }
  # prim13@assays$SCT@data=prim13.sct.data.bak
  ######
  # @data
  # saveRDS(markerlist, "markerlist.prim13.pairwise.initiating.types.rds")
  # @scale.data
  # saveRDS(markerlist, "markerlist.prim13.pairwise.initiating.types.scaled.data.rds")


  # @data
  markerlist=readRDS("markerlist.prim13.pairwise.initiating.types.rds")
  # @scale.data
  # markerlist=readRDS("markerlist.prim13.pairwise.initiating.types.scaled.data.rds")
  ######
  celltypes=sort(unique(prim13@meta.data$ori))
  markerlist.filter=list()
  for(j in celltypes)
  {
    print("celltype:")
    print(j)
    
    print("------------>>>>>>>>>>>")
    sp.markers=rownames(prim13@assays$SCT@scale.data)
    for(i in names(markerlist))
    {
      temp=strsplit(i, "_")[[1]][1]
      if(temp==j)
      {
        markers=markerlist[[i]]
        markers=markers[markers$p_val_adj<0.05, ]
        markers=markers[order(markers$avg_logFC, decreasing = TRUE), ]
        markers=markers[markers$avg_logFC>log(1.5) & markers$pct.1>0.2, ] # markers$avg_logFC>log(1.5) & markers$pct.1>0.2 # 
        sp.markers=intersect(sp.markers, rownames(markers))
        print("pairwise:")
        print(i)
        print("length(sp.markers):")
        print(length(sp.markers))
      }
    }
    markerlist.filter[[j]]=sp.markers
  }
  for(i in names(markerlist.filter))
  {
    markerlist.filter[[i]]=markerlist.filter[[i]][1:10]
    markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
  }
  markers=Reduce(union, markerlist.filter)
  str(markers)
  # markers$avg_logFC>log(1.5) & markers$pct.1>0.2 # This s the BEST group
  # chr [1:552] "RBP7" "IFI6" "TIE1" "PIK3R3" "PLPP3" "TGFBR3" "S1PR1" "KCNN3" ...
  # @scale.data markers$avg_logFC>log(1.6)
  # chr [1:468] "MBOAT2" "TRIB2" "HES1" "APC" "PNISR" "RARRES3" "N4BP2L2" ...
  # @scale.data markers$avg_logFC>log(1.4)
  # chr [1:516] "MBOAT2" "TRIB2" "HES1" "APC" "PNISR" "AKAP9" "RARRES3" ...
  # @scale.data markers$avg_logFC>log(1.5)
  # chr [1:493] "MBOAT2" "TRIB2" "HES1" "APC" "PNISR" "AKAP9" "RARRES3" ...
  # @scale.data markers$avg_logFC>log(2)
  # chr [1:929] "MOB3C" "MBOAT2" "KLF11" "TRIB2" "HES1" "PLRG1" "APC" ...
  # markers$avg_logFC>log(1.75)
  # chr [1:322] "RBP7" "PLPP3" "CD34" "VAMP5" "CALCRL" "CAVIN2" "TGFBR2" ...
  # markers$avg_logFC>0.25 & markers$pct.1>0.2
  # chr [1:1038] "TMEM63A" "RBP7" "CDA" "ID3" "IFI6" "TIE1" "NASP" "PIK3R3" ...
  pdf(paste("temp.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/10 )
  temp.heatmap<-DoHeatmap( subset(prim13, downsample=300), features = markers, group.by = "ori", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()

  # @data-based markers
  # write.table(markers, file="marker.union.prim13.pairwise.initiating.types.txt", col.name=F, row.name=F, quote=F, sep='\t')
  # @scale.data-based markers:
  write.table(markers, file="markerlist.prim13.pairwise.initiating.types.scaled.data.txt", col.name=F, row.name=F, quote=F, sep='\t')

  # @data-based markers
  # markers=read.table("marker.union.prim13.pairwise.initiating.types.txt", sep='\t')$V1
  # @scale.data-based markers:
  markers=read.table("markerlist.prim13.pairwise.initiating.types.scaled.data.txt", sep='\t')$V1
  # make Cibersort's Signature matrix
  celltypes=sort(unique(prim13@meta.data$ori))
  # [1] "Basal"            "Endothelial cell" "Fibroblast"       "Luminal"
  # [5] "Macrophage"       "Mast cell"        "Myofibroblast"    "T cell"
  # [9] "TypeC"
  # create signature matrix "sigm"
  sigm=prim13@assays$SCT@scale.data[markers, ]
  colnames(sigm)=prim13@meta.data$ori
  sigm=t(sigm)

  sigm=aggregate(sigm, by = list(rownames(sigm)), FUN = mean, na.rm = TRUE)
  sigm[1:3, 1:3]
  rownames(sigm)=sigm$"Group.1"
  sigm=sigm[, c(-1)]
  sigm=as.matrix(sigm)
  sigm=t(sigm)
  sigm[1:3, 1:3]
  # saveRDS(sigm, "prim13.ori.average.sig.matrix.rds")
  # saveRDS(sigm, "prim13.ori.average.sig.matrix.scale.data.rds")

  ### use FARDEEP to deconvolution
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522071/
  # Rscript stage.R ../prim13.ori.average.sig.matrix.rds fardeep no no TypeC TypeC
  # Rscript survival2.R ../prim13.ori.average.sig.matrix.rds fardeep no no TypeC TypeC
}


#### split indolent and aggresive samples CNV intensity on violin plot
epi@meta.data$celltype.split=epi@meta.data$epitype
epi@meta.data$celltype.split[epi@meta.data$pheno=="Indo"]=paste("Indo", epi@meta.data$celltype.split[epi@meta.data$pheno=="Indo"], sep=".")
epi@meta.data$celltype.split[epi@meta.data$pheno=="Aggr"]=paste("Aggr", epi@meta.data$celltype.split[epi@meta.data$pheno=="Aggr"], sep=".")
table(epi@meta.data$epitype)
table(epi@meta.data$celltype.split)

pdf("temp.VlnPlot.CNV_COUNTS.IndoAggrSplited.pdf", width=20, height=8)
check.genes="cnv.counts"
temp=VlnPlot(epi, features = check.genes, pt.size = 0, ncol = 1, group.by = "celltype.split")
for(i in 1:length(check.genes)) 
{
  temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
}
print(temp)
dev.off()


### cell communication with NicheNet
{
  ### start from comparing GS-high luminal cells and GS-low luminal cells
  # markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
  # epi=readRDS("epi.seurat3.indivSCTed.celltyped.rds")
  # epi13=readRDS("epi13.seurat3.indivSCTed.celltyped.rds")

  table( epi13@meta.data$epitype )
    # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
    # 213   165   109   115   219    87   632 11931  6763  4929   446
  table( epi13@meta.data$gleason )
    #  10   3+4     7     9
    # 267  2617  6855 15870
  lum13=subset(epi13, subset=(epitype%in%c("tL0", "tL1", "tL2", "tLI3")))
  # compare Gleason 9/10 VS. Gleason 7's luminal cells
  {
    lum13@meta.data$gleason.easy=rep( "low", nrow(lum13@meta.data) )
    lum13@meta.data$gleason.easy[lum13@meta.data$gleason%in%c("9", "10")]="high"
    lum13=SetIdent(lum13, value = "gleason.easy")
    markers.lum13.gs=FindMarkers(lum13, ident.1 = "high", ident.2 = "low", assay="SCT")
    table(lum13@meta.data$gleason)
    #  10   3+4     7     9
    # 218  2538  6656 14657
    table(lum13@meta.data$gleason.easy)
    #  high   low
    # 14875  9194
    # saveRDS(markers.lum13.gs, "markers.lum13.gs.rds")
    # markers.lum13.gs=readRDS("markers.lum13.gs.rds")

    markers=markers.lum13.gs
    markers=markers[order(markers$avg_logFC, decreasing = TRUE), ]
    markers$foldchange=exp(markers$avg_logFC)
    markers.1=markers[markers$p_val_adj<0.05  & markers$avg_logFC>log(1.5) & markers$pct.1>0.1,]#
    markers.2=markers[markers$p_val_adj<0.05  & markers$avg_logFC<(-log(1.5)) & markers$pct.2>0.1,]
    markers=rbind(markers.1, markers.2)
    str(markers)
    # 'data.frame':    562 obs. of  6 variables:
    write.csv(markers[,c("foldchange","pct.1", "pct.2", "p_val_adj")], file=paste("markers.lum13.gs","csv",sep="."), row.name=T)

    # lum13 <- ScaleData(lum13, features = rownames(lum13@assays$RNA@data), assay="RNA" )
    # lum13 <- ScaleData(lum13, features = rownames(lum13@assays$SCT@data), assay="SCT" )
    temp.markers=markers[c(1:20, (nrow(markers)-19):nrow(markers)), ]
    lum13=SetIdent(lum13, value = "gleason.easy")
    pdf(paste("lum13.gs.top20.pdf", sep=""),  width=(2)/4+5, height=nrow(temp.markers)/8)
    temp.heatmap<-DoHeatmap( subset(lum13, downsample=500), features = rownames(temp.markers), group.by = "gleason.easy", slot="scale.data", assay="RNA", angle = 0 , size=2)
    # , subset=(fig.cluster%in%temp.clusters)
    print(temp.heatmap)
    dev.off()

    # enrichment
    str(rownames(markers)[markers$avg_logFC>0])
    enrichment.gene.set.plot(rownames(markers)[markers$avg_logFC>0], topNterm=20)
    str(rownames(markers)[markers$avg_logFC<0])
    enrichment.gene.set.plot(rownames(markers)[markers$avg_logFC<0], topNterm=20)
  }

  ### Gleason 9/10 VS. 7 NicheNet to do cell communication analysis
  {
    require("nichenetr")
    require("tidyverse")

    ### Step 1: Define expressed genes in sender and receiver cell populations
    # get expressed genes in high-GS luminal, low-GS luminal, all-typeC-origining cells

    # genes expressed in >10% high-GS luminal cells
    table(lum13@meta.data$gleason)
      #  10   3+4     7      9
      # 218  2538  6656  14657
    png("dimplot.gleason.png", height=5*100, width=5*100)
    DimPlot(lum13, label = TRUE, group.by="gleason", label.size=8) # + NoLegend()
    dev.off()

    png("dimplot.cluster.png", height=5*100, width=5*100)
    DimPlot(lum13, label = TRUE, group.by="fig.cluster", label.size=8) # + NoLegend()
    dev.off()

    high.gs.lum13=subset(lum13, subset=(gleason.easy%in%c("high")))
    temp=high.gs.lum13
    temp.ncol=ncol(temp@assays$SCT@counts)
    temp=as.matrix(temp@assays$SCT@counts)
    temp=temp[rowMaxs(temp)>0, ]
    temp.exped.genes.pct=apply(temp, 1, function(x) (table(x>0)[["TRUE"]]/temp.ncol))
    temp.exped.genes=rownames(temp)[temp.exped.genes.pct>0.1]
    high.gs.lum13.exped.genes=temp.exped.genes
    str(high.gs.lum13.exped.genes)
     # chr [1:7163] "NOC2L" "HES4" "ISG15" "AGRN" "SDF4" "B3GALT6" "INTS11" ...
    # saveRDS(high.gs.lum13.exped.genes, "high.gs.lum13.exped.genes.rds")

    # genes expressed in >10% low-GS luminal cells
    low.gs.lum13=subset(lum13, subset=(gleason.easy%in%c("low")))
    temp=low.gs.lum13
    temp.ncol=ncol(temp@assays$SCT@counts)
    temp=as.matrix(temp@assays$SCT@counts)
    temp=temp[rowMaxs(temp)>0, ]
    temp.exped.genes.pct=apply(temp, 1, function(x) (table(x>0)[["TRUE"]]/temp.ncol))
    temp.exped.genes=rownames(temp)[temp.exped.genes.pct>0.1]
    low.gs.lum13.exped.genes=temp.exped.genes
    str(low.gs.lum13.exped.genes)
     # chr [1:4703] "NOC2L" "HES4" "SDF4" "INTS11" "AURKAIP1" "CCNL2" "MRPL20" ...
    # saveRDS(low.gs.lum13.exped.genes, "low.gs.lum13.exped.genes.rds")
    expressed_genes_receiver=intersect(low.gs.lum13.exped.genes, high.gs.lum13.exped.genes)
    str(expressed_genes_receiver)
     # chr [1:4626] "NOC2L" "HES4" "SDF4" "INTS11" "AURKAIP1" "CCNL2" "MRPL20" ...

    ### (Not employed) Define expressed_genes_sender
    ### Version 1, genes expressed in >10% cells as expressed_genes_sender
    {
      # genes expressed in >10% typeC-origining cells
      typec13=subset(epi13, subset=(epitype%in%c("tB1", "tB2", "tBC3")))
      temp=typec13
      temp.ncol=ncol(temp@assays$SCT@counts)
      temp=as.matrix(temp@assays$SCT@counts)
      temp=temp[rowMaxs(temp)>0, ]
      temp.exped.genes.pct=apply(temp, 1, function(x) (table(x>0)[["TRUE"]]/temp.ncol))
      temp.exped.genes=rownames(temp)[temp.exped.genes.pct>0.1]
      typec.exped.genes=temp.exped.genes
      str(typec.exped.genes)
       # chr [1:4644] "NOC2L" "HES4" "ISG15" "AGRN" "SDF4" "INTS11" "AURKAIP1" ...
      # saveRDS(typec.exped.genes, "typec.exped.genes.rds")
      expressed_genes_sender=typec.exped.genes
    }

    ### (Employed) Version 2: genes up-regulated in typeC-origining cells compared to high/low-Gleason luminal cells as expressed_genes_sender
    {
      # epi.indagg=readRDS("epi.indagg.seurat3.indivSCTed.celltyped.rds")
      # epi.crpc=readRDS("epi.crpc.seurat3.indivSCTed.celltyped.rds")
      # normale=readRDS("normale.seurat3.indivSCTed.celltyped.rds")
      # epi13=readRDS("epi13.seurat3.indivSCTed.celltyped.rds")
      epi=merge(x = normale, y = list(et=epi.indagg, t=epi13, crpc=epi.crpc))
      str(epi)
      str(epi@assays$SCT)
      epi@assays$SCT@scale.data=cbind(normale@assays$SCT@scale.data, epi.indagg@assays$SCT@scale.data, epi13@assays$SCT@scale.data, epi.crpc@assays$SCT@scale.data)
      table(epi@meta.data$epitype)
       #       B2        B3        B6        C4      crI5      crL0      crL1      crL2
       #     4056      4001       435      1992       256      3977       910       861
       #     crL4     crNE3      etB1      etC1  etC4.AGG  etC5.AGG  etI1.AGG  etI3.IND
       #      662       685      2746      4308       747       413      1282      1526
       # etL0.AGG  etL0.IND  etL2.AGG etNE6.AGG        I5        I7        L0        L1
       #      517      8237      2890       298       603       123      6268      6201
       #       L8       tB1       tB2       tB4      tBC3      tBI0      tBI5      tCyc
       #       60       213       165       109       115       219        87       632
       #      tL0       tL1       tL2      tLI3
       #    11931      6763      4929       446
      # genes up-regulated in tumor type C compared with normal type C + Basal
      epi@meta.data$for_temp_compare=rep("other", ncol(epi))
      epi@meta.data$for_temp_compare[epi@meta.data$epitype%in%c("tB1", "tB2", "tBC3")]="tC"
      epi@meta.data$for_temp_compare[epi@meta.data$epitype%in%c("B2", "B3", "C4", "B6")]="nBC"
      epi@meta.data$for_temp_compare[epi@meta.data$epitype%in%c("tL0", "tL1", "tL2", "tLI3")]="tL"
      table(epi@meta.data$for_temp_compare)
      #   nBC other    tC    tL
      # 10484 44617   493 24069
      epi=SetIdent(epi, value="for_temp_compare")
      tC.up.1=FindMarkers(epi, ident.1 = "tC", ident.2 = "nBC", assay="SCT")
      # saveRDS(tC.up.1, "tC.up.1.rds")
      tC.up.2=FindMarkers(epi, ident.1 = "tC", ident.2 = "tL", assay="SCT")
      # saveRDS(tC.up.2, "tC.up.2.rds")
      # tC.up.2=readRDS("tC.up.2.rds")
      tC.up.genes.1=rownames(tC.up.1)[tC.up.1$pct.1>0.1 & tC.up.1$pct.2<0.1 & tC.up.1$avg_logFC>log(1.5)]
      tC.up.genes.2=rownames(tC.up.2)[tC.up.2$pct.1>0.1 & tC.up.2$pct.2<0.05 & tC.up.2$avg_logFC>log(1.5)]
      # tC.up.genes=intersect(tC.up.genes.1, tC.up.genes.2)
      # expressed_genes_sender=tC.up.genes
      expressed_genes_sender=tC.up.genes.2
      str(expressed_genes_sender)
      #  chr [1:87] "KRT17" "MMP7" "WFDC2" "CCL2" "OLFM4" "KRT15" "GSTP1" "S100A2" ...
    }

    ### (Not employed) Version 3: up-regulated genes in Type C tumor cells than any cell type (including immune cell types and stromal cell types) as expressed_genes_sender
    {
      epi13=readRDS("epi13.seurat3.indivSCTed.celltyped.rds")
      ### mark epi13's cells as tumor-initiating celltypes
      table(epi13@meta.data$epitype)
        # tB1   tB2   tB4  tBC3  tBI0  tBI5  tCyc   tL0   tL1   tL2  tLI3
        # 213   165   109   115   219    87   632 11931  6763  4929   446
      current.cluster.ids <- c("tB1", "tB2", "tB4", "tBC3", "tBI0", "tBI5", "tCyc", "tL0", "tL1", "tL2", "tLI3")
      new.cluster.ids <- c("TypeC", "TypeC", "Basal", "TypeC", "Basal", "Basal", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal") # c("Peripheral", "Transition", "Central")
      epi13@meta.data$ori <- plyr::mapvalues(x = as.character(epi13@meta.data$epitype), from = current.cluster.ids, to = new.cluster.ids)
      epi13@meta.data$ori=as.factor(epi13@meta.data$ori)
      epi13@meta.data$ori=as.character(epi13@meta.data$ori)
      table(epi13@meta.data$ori)
        # Basal Luminal   TypeC
        #   415   24701     493
      # consider celltypes in all cell types
      prim13=readRDS("prim13.seurat3.indivSCTed.celltyped.rds")
      table(prim13@meta.data$celltype)
          #   Basal cell Endothelial cell  Epithelial cell       Fibroblast
          #          964             3212              651              408
          # Luminal cell       Macrophage        Mast cell    Myofibroblast
          #        23994             2142             1443             1798
          #       T cell
          #         3135
      prim13@meta.data$ori=prim13@meta.data$celltype
      epi13.ori=epi13@meta.data$ori
      names(epi13.ori)=rownames(epi13@meta.data)
      table(epi13.ori)
        # Basal Luminal   TypeC
        #   415   24701     493
      prim13@meta.data[names(epi13.ori), "ori"]=epi13.ori
      table(prim13@meta.data$ori)
            #      Basal Endothelial cell       Fibroblast          Luminal
            #        415             3212              408            24701
            # Macrophage        Mast cell    Myofibroblast           T cell
            #       2142             1443             1798             3135
            #      TypeC
            #        493
      prim13=SetIdent(prim13, value = "ori")
      markerlist=readRDS("markerlist.prim13.pairwise.initiating.types.rds")
      ######
      celltypes=sort(unique(prim13@meta.data$ori))
      print(celltypes)
      # [1] "Basal"            "Endothelial cell" "Fibroblast"       "Luminal"
      # [5] "Macrophage"       "Mast cell"        "Myofibroblast"    "T cell"
      # [9] "TypeC"      
      j="TypeC"
      {
        print("celltype:")
        print(j)
        print("------------>>>>>>>>>>>")
        sp.markers=rownames(prim13@assays$SCT@scale.data)
        for(i in names(markerlist))
        {
          temp=strsplit(i, "_")[[1]][1]
          temp.2=strsplit(i, "_")[[1]][2]
          if((temp==j)&(temp.2!="Basal"))
          {
            markers=markerlist[[i]]
            markers=markers[markers$p_val_adj<0.05, ]
            markers=markers[order(markers$avg_logFC, decreasing = TRUE), ]
            markers=markers[markers$avg_logFC>log(1.5) & markers$pct.1>0.1, ] # markers$avg_logFC>log(1.5) & markers$pct.1>0.2 # 
            sp.markers=intersect(sp.markers, rownames(markers))
            print("pairwise:")
            print(i)
            print("length(sp.markers):")
            print(length(sp.markers))
          }
        }
      }
      str(sp.markers)
      # chr [1:26] "EPHA2" "GADD45A" "S100A6" "CRABP2" "LAMB3" "FHL2" "RARRES1" ...
      expressed_genes_sender=sp.markers
    }

    ### Step 2: Define the gene set of interest and a background of genes
    # gene set of interests as DEGs between high- and low-GS luminal cells
    {
      markers.lum13.gs=readRDS("markers.lum13.gs.rds")
      markers=markers.lum13.gs
      markers=markers[order(markers$avg_logFC, decreasing = TRUE), ]
      markers$foldchange=exp(markers$avg_logFC)
      markers.1=markers[markers$p_val_adj<0.05  & markers$avg_logFC>log(1.5) & markers$pct.1>0.1,]#
      markers.2=markers[markers$p_val_adj<0.05  & markers$avg_logFC<(-log(1.5)) & markers$pct.2>0.1,]
      markers=rbind(markers.1, markers.2)
      str(rownames(markers))
       # chr [1:562] "OGN" "IGF1" "KLK11" "GAPDH" "RPS3" "RPL7A" "RPS18" "RPS2" ...
      geneset_oi=rownames(markers.1)
      # geneset_oi=rownames(markers.2)
      str(geneset_oi)
      # chr [1:553] "OGN" "IGF1" "KLK11" "GAPDH" "RPS3" "RPL7A" "RPS18" "RPS2" ...
      # chr [1:9] "CLU" "FAM3B" "GDF15" "AZGP1" "MT1G" "LINC01088" "VGLL3" ...

      # read ligand-target regulation matirx
      ligand_target_matrix = readRDS("./NicheNet/ligand_target_matrix.rds")
      background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
      str(background_expressed_genes)
      # chr [1:4575] "NOC2L" "HES4" "SDF4" "INTS11" "AURKAIP1" "CCNL2" "MRPL20" ...
    }

    # Step 3: Define a set of potential ligands
    {
      lr_network = readRDS("./NicheNet/lr_network.rds")
      # If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions and only keep ligand-receptor interactions that are described in curated databases. To do this: uncomment following line of code:
      # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
      ligands = lr_network %>% pull(from) %>% unique()
      expressed_ligands = intersect(ligands, expressed_genes_sender)
      receptors = lr_network %>% pull(to) %>% unique()
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
      lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
      str(lr_network_expressed)
      # tibble [20 x 4] (S3: tbl_df/tbl/data.frame)
      #  $ from    : chr [1:20] "AREG" "AREG" "LAMC2" "LTF" ...
      #  $ to      : chr [1:20] "EGFR" "EGFR" "ITGB1" "LRP11" ...
      potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
      str(potential_ligands)
      # chr [1:10] "AREG" "LAMC2" "LTF" "CTGF" "PLAU" "SAA1" "GSTP1" "NTN4" ...
    }

    # Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest
    {
      ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
      ligand_activities %>% arrange(-pearson) 
      # # A tibble: 3 x 4
      #   test_ligand auroc  aupr pearson
      #   <chr>       <dbl> <dbl>   <dbl>
      # 1 SAA1        0.575 0.167  0.0988
      # 2 CX3CL1      0.561 0.162  0.0866
      # 3 PLAU        0.562 0.157  0.0812
      best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
      str(best_upstream_ligands)
       # chr [1:3] "SAA1" "CX3CL1" "PLAU"
      write.table(best_upstream_ligands, file="best_upstream_ligands.txt", col.name=F, row.name=F, quote=F, sep='\t')

      # normale=readRDS("normale.seurat3.indivSCTed.celltyped.rds")
      # epi13=readRDS("epi13.seurat3.indivSCTed.celltyped.rds")
      # prim13=readRDS("prim13.seurat3.indivSCTed.celltyped.rds")
      tocheck=best_upstream_ligands
      pdf("VlnPlot.ligand_celltype.normale.pdf", width = 10, height =ceiling(length(tocheck)/2)*2)
      temp=VlnPlot(prim13, features = tocheck, pt.size = 0, group.by="ori", ncol = 2)
      for(i in 1:length(tocheck)) 
      {
        temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
      }
      print(temp)
      dev.off()

      pdf("VlnPlot.ligand_celltype.epi13.pdf", width = 10, height =ceiling(length(tocheck)/2)*2)
      temp=VlnPlot(epi13, features = tocheck, pt.size = 0, group.by="epitype", ncol = 2)
      for(i in 1:length(tocheck)) 
      {
        temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
      }
      print(temp)
      dev.off()

      # show histogram of ligand activity scores
      p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
        geom_histogram(color="black", fill="darkorange")  + 
        # geom_density(alpha=.1, fill="orange") +
        geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
        labs(x="ligand activity (PCC)", y = "# ligands") +
        theme_classic()
      print(p_hist_lig_activity)
      dev.off()


      # Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap
      active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
      active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
      order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
      order_targets = active_ligand_target_links_df$target %>% unique()
      vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

      p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands in TypeC-initiating tumor cells","High-Gleason related genes in luminal cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
      pdf("ligand.2.targets.regulation.pdf", width=10, height=8)
      print(p_ligand_target_network)
      dev.off()

      write.table(order_targets, file="order_targets.txt", col.name=F, row.name=F, quote=F, sep='\t')
      enrichment.gene.set.plot(order_targets, topNterm=20)

      temp.markers=order_targets
      lum13=SetIdent(lum13, value = "gleason.easy")
      pdf(paste("lum13.gs.top20.pdf", sep=""),  width=(2)/4+5, height=length(temp.markers)/8)
      temp.heatmap<-DoHeatmap( subset(lum13, downsample=500), features = temp.markers, group.by = "gleason.easy", slot="scale.data", assay="SCT", angle = 0 , size=2)
      # , subset=(fig.cluster%in%temp.clusters)
      print(temp.heatmap)
      dev.off()

      # check expression of strong regulated targets in lum13 across high- and low-Gleason samples
      tocheck=c("AR")
      pdf("VlnPlot.target_celltype.epi13.pdf", width = 10, height =30)
      temp=VlnPlot(lum13, features = tocheck, pt.size = 0, group.by="gleason.easy", ncol = 2)
      for(i in 1:length(tocheck)) 
      {
        temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
      }
      print(temp)
      dev.off()
    }

    # Follow-up analysis 1: Ligand-receptor network inference for top-ranked ligands
    {
      # get the ligand-receptor network of the top-ranked ligands
      lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
      best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

      # get the weights of the ligand-receptor interactions as used in the NicheNet model
      weighted_networks = readRDS("./NicheNet/weighted_networks.rds")
      lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

      # convert to a matrix
      lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
      lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

      # perform hierarchical clustering to order the ligands and receptors
      dist_receptors = dist(lr_network_top_matrix, method = "binary")
      hclust_receptors = hclust(dist_receptors, method = "ward.D2")
      order_receptors = hclust_receptors$labels[hclust_receptors$order]

      dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
      hclust_ligands = hclust(dist_ligands, method = "ward.D2")
      order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

      vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
      p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized ligands in TypeC-initiating tumor cells", "Receptors expressed by tumor luminal cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

      pdf("ligand.2.receptor.regulation.pdf", width=10, height=8)
      print(p_ligand_receptor_network)
      dev.off()

      temp.markers=best_upstream_receptors
      lum13=SetIdent(lum13, value = "gleason.easy")
      pdf(paste("lum13.gs.top20.pdf", sep=""),  width=(2)/4+5, height=length(temp.markers)/8)
      temp.heatmap<-DoHeatmap( subset(lum13, downsample=500), features = temp.markers, group.by = "gleason.easy", slot="scale.data", assay="SCT", angle = 0 , size=2)
      # , subset=(fig.cluster%in%temp.clusters)
      print(temp.heatmap)
      dev.off()

      write.table(best_upstream_receptors, file="best_upstream_receptors.txt", col.name=F, row.name=F, quote=F, sep='\t')
    }

  }
}


# get type C tumor cell surface marker
markerlist.prim13=readRDS("markerlist.prim13.pairwise.initiating.types.rds")
markerlist=markerlist.prim13
celltypes=sort(unique(prim13@meta.data$ori))
print(celltypes)
# [1] "Basal"            "Endothelial cell" "Fibroblast"       "Luminal"
# [5] "Macrophage"       "Mast cell"        "Myofibroblast"    "T cell"
# [9] "TypeC"      
j="TypeC"
{
  print("celltype:")
  print(j)
  print("------------>>>>>>>>>>>")
  sp.markers=rownames(prim13@assays$SCT@scale.data)
  for(i in names(markerlist))
  {
    temp=strsplit(i, "_")[[1]][1]
    temp.2=strsplit(i, "_")[[1]][2]
    if( (temp==j) ) # & (temp.2!="Basal") )
    {
      markers=markerlist[[i]]
      markers=markers[markers$p_val_adj<0.05, ]
      markers=markers[order(markers$avg_logFC, decreasing = TRUE), ]
      markers=markers[markers$avg_logFC>log(1.5) & markers$pct.1>0.5, ] # & markers$pct.2<0.01
      sp.markers=intersect(sp.markers, rownames(markers))
      print("pairwise:")
      print(i)
      print("length(sp.markers):")
      print(length(sp.markers))
    }
  }
}
str(sp.markers)
#  chr [1:11] "GADD45A" "LAMB3" "FHL2" "MAST4" "SQSTM1" "PHLDA2" "GSTP1" ...
markerlist.epi=readRDS("epi.markerlist.pairwise.rds")
markerlist=markerlist.epi$"C4__tBC3"
markerlist=markerlist[order(markerlist$avg_logFC),]
markerlist=markerlist[markerlist$avg_logFC<(-log(1.5))&markerlist$pct.1<0.01,]#&markerlist$pct.2>0.1
deg.tBC3.C4=rownames(markerlist)
str(deg.tBC3.C4)
# chr [1:8] "MIF" "PCBP1" "TUBA1B" "HES4" "ARL6IP4" "GABARAP" "PRELID1" ...
intersect(sp.markers, deg.tBC3.C4) # surface: MIF
# character(0)

# keep cell-surface markers
require(msigdbr)
cc = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
cc = cc %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
cdgenes=sort(unique(cc$gene_symbol[cc$gs_name=="GO_CELL_SURFACE"]))
sp.markers=sp.markers[sp.markers%in%cdgenes]

pdf(paste("lum13.gs.top20.pdf", sep=""),  width=6, height=6)
temp.heatmap<-DoHeatmap( subset(epi13, downsample=500), features = sp.markers, group.by = "ori.type", slot="scale.data", assay="SCT", angle = 0 , size=2)
# , subset=(fig.cluster%in%temp.clusters)
print(temp.heatmap)
dev.off()


### additionally consider normal typeC with all primary PCa cell types together to make Deconvolution Signature Matrix
{
  ### deconvlution to validate association between percentage of specific celltype and clinical features, e.g. Gleason, T stage, N stage, M stage, and BCR/OS.
  ### deconvlution's signature preparation
  # normale=readRDS("normale.seurat3.indivSCTed.celltyped.rds")
  # prim13=readRDS("prim13.seurat3.indivSCTed.celltyped.rds")
  ### combine normal TypeC and to 13-PCa data
  temp=subset(normale, subset=(celltype=="C4"))
  prim13plus=merge(x = prim13, y = temp)
  str(prim13plus)
  str(prim13plus@assays$SCT)
  prim13plus@assays$SCT@scale.data=cbind(prim13@assays$SCT@scale.data, temp@assays$SCT@scale.data)
  table(is.na(prim13plus@meta.data$ori))
  prim13plus@meta.data$ori[is.na(prim13plus@meta.data$ori)]="NormC"
  # pairwise comparison
  celltypes=sort(unique(prim13plus@meta.data$ori))
  markerlist.plus=list()
  prim13plus=SetIdent(prim13plus, value = "ori")
  # for(i in celltypes)
  i="NormC"
  {
    for(j in celltypes)
    {
      if(i!=j)
      {
        tempname=paste(i, j, sep="_")
        print("tempname:")
        print(tempname)
        markerlist.plus[[tempname]]=FindMarkers(prim13plus, min.pct = 0, ident.1 = i, ident.2=j, assay="SCT")
      }
    }
  }
  markerlist.mirror=list()
   # for(i in celltypes)
  i="NormC"
  {
    for(j in celltypes)
    {
      if(i!=j)
      {
        tempname=paste(j, i, sep="_")
        print("tempname:")
        print(tempname)
        markerlist.mirror[[tempname]]=FindMarkers(prim13plus, min.pct = 0, ident.1 = j, ident.2=i, assay="SCT")
      }
    }
  }
  markerlist=readRDS("markerlist.prim13.pairwise.initiating.types.rds")
  markerlist=c(markerlist, markerlist.plus, markerlist.mirror)
  # saveRDS(markerlist, "markerlist.plus.prim13.pairwise.initiating.types.rds")

  # get celltype-specific markers
  markerlist=readRDS("markerlist.plus.prim13.pairwise.initiating.types.rds")
  {
    markerlist=readRDS("markerlist.plus.prim13.pairwise.initiating.types.rds")
    celltypes=sort(unique(prim13plus@meta.data$ori))
    markerlist.filter=list()
    for(j in celltypes)
    {
      print("celltype:")
      print(j)
      print("------------>>>>>>>>>>>")
      sp.markers=rownames(prim13plus@assays$SCT@scale.data)
      for(i in names(markerlist))
      {
        temp=strsplit(i, "_")[[1]][1]
        if(temp==j)
        {
          markers=markerlist[[i]]
          markers=markers[markers$p_val_adj<0.05, ]
          markers=markers[order(markers$avg_logFC, decreasing = TRUE), ]
          markers=markers[markers$avg_logFC>log(1.5), ] # markers$avg_logFC>log(1.5) & markers$pct.1>0.2 # 
          sp.markers=intersect(sp.markers, rownames(markers))
          print("pairwise:")
          print(i)
          print("length(sp.markers):")
          print(length(sp.markers))
        }
      }
      markerlist.filter[[j]]=sp.markers
    }
    # for(i in names(markerlist.filter))
    # {
    #   markerlist.filter[[i]]=markerlist.filter[[i]][1:10]
    #   markerlist.filter[[i]]=markerlist.filter[[i]][!is.na(markerlist.filter[[i]])]
    # }
    markers=Reduce(union, markerlist.filter)
    str(markers)
  }

  pdf(paste("temp.heatmap.pdf", sep=""), width=length(celltypes)/3+5, height=length(markers)/10 )
  temp.heatmap<-DoHeatmap( subset(prim13plus, downsample=300), features = markers, group.by = "ori", slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  print(temp.heatmap)
  dev.off()
  
  # write.table(markers, file="markerlist.prim13plus.pairwise.initiating.types.scaled.data.txt", col.name=F, row.name=F, quote=F, sep='\t')
  # markers=read.table("markerlist.prim13plus.pairwise.initiating.types.scaled.data.txt", sep='\t')$V1

  # make Cibersort's Signature matrix
  celltypes=sort(unique(prim13plus@meta.data$ori))
  # [1] "Basal"            "Endothelial cell" "Fibroblast"       "Luminal"
  # [5] "Macrophage"       "Mast cell"        "Myofibroblast"    "T cell"
  # [9] "TypeC"
  # create signature matrix "sigm"
  sigm=prim13plus@assays$SCT@scale.data[markers, ]
  colnames(sigm)=prim13plus@meta.data$ori
  sigm=t(sigm)

  sigm=aggregate(sigm, by = list(rownames(sigm)), FUN = mean, na.rm = TRUE)
  sigm[1:3, 1:3]
  rownames(sigm)=sigm$"Group.1"
  sigm=sigm[, c(-1)]
  sigm=as.matrix(sigm)
  sigm=t(sigm)
  sigm[1:3, 1:3]
  # saveRDS(sigm, "prim13plus.ori.average.sig.matrix.rds")
  # saveRDS(sigm, "prim13plus.ori.average.sig.matrix.scale.data.rds")

  ### use FARDEEP to deconvolution
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522071/
  # Rscript stage.R ../prim13plus.ori.average.sig.matrix.rds fardeep no no TypeC TypeC
  # Rscript stage.R ../prim13plus.ori.average.sig.matrix.rds fardeep no no NormC NormC
  # Rscript survival2.R ../prim13plus.ori.average.sig.matrix.rds fardeep no no TypeC TypeC
}


### cell state transition trajectory analysis
### use Monocle2
# refer to: http://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle
# BiocManager::install("monocle")
require(monocle)
packageVersion("monocle")
# [1] ‘2.16.0’
# epi.indagg=readRDS("epi.indagg.seurat3.indivSCTed.celltyped.rds")
# epi.crpc=readRDS("epi.crpc.seurat3.indivSCTed.celltyped.rds")
# normale=readRDS("normale.seurat3.indivSCTed.celltyped.rds")
# epi13=readRDS("epi13.seurat3.indivSCTed.celltyped.rds")
epi=merge(x = normale, y = list(et=epi.indagg, t=epi13, crpc=epi.crpc))
epi@assays$SCT@scale.data=cbind(normale@assays$SCT@scale.data, epi.indagg@assays$SCT@scale.data, epi13@assays$SCT@scale.data, epi.crpc@assays$SCT@scale.data)
# get cross-dataset pariwise DEGs
markerlist.epi=readRDS("epi.markerlist.pairwise.rds")

epi = SetIdent(epi, value = "epitype")
{

  typec=subset( epi, subset=(epitype%in%c("tBC3", "C4")) )
  table(typec@meta.data$epitype)
  typec@meta.data$epitype[typec@meta.data$epitype=="C4"]="NormC"
  typec@meta.data$epitype[typec@meta.data$epitype=="tBC3"]="TumorC"
  gene_metadata=data.frame(gene_short_name=rownames(typec@assays$SCT@data))
  rownames(gene_metadata)=rownames(typec@assays$SCT@data)
  typec.evo <- newCellDataSet(  typec@assays$SCT@data,
                            phenoData = new("AnnotatedDataFrame", typec@meta.data),
                            featureData = new("AnnotatedDataFrame", gene_metadata), 
                            expressionFamily=negbinomial.size() )
  table(typec.evo@phenoData@data$epitype)
  table(typec.evo@phenoData@data$pheno)
  typec.evo <- estimateSizeFactors(typec.evo)
  typec.evo <- estimateDispersions(typec.evo)

  typec.evo <- detectGenes(typec.evo, min_expr = 0.1)
  ### only keep expressed genes
  expressed_genes <- row.names(typec.evo)[typec.evo@featureData@data$num_cells_expressed>= 10]
  typec.evo <- typec.evo[expressed_genes,]

  ### use all significant markers of clusters as ordering genes
  tumor.markers.ori=markerlist.epi$"C4__tBC3"
  tumor.markers=tumor.markers.ori
  tumor.markers$pct.1=tumor.markers.ori$pct.2
  tumor.markers$pct.2=tumor.markers.ori$pct.1
  tumor.markers$avg_logFC=-1*tumor.markers$avg_logFC

  markers=tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[abs(markers$avg_logFC)>log(1.5), ]
  markers$foldChange=exp(markers$avg_logFC)
  # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
  markers=markers[order(markers$foldChange, decreasing=TRUE), ]
  # ordering.genes=unique(markers$gene)
  ordering.genes=unique(rownames(markers))

  typec.evo <- setOrderingFilter(typec.evo,  ordering.genes) # 
  pdf("plot_ordering_genes.pdf")
  plot_ordering_genes(typec.evo)
  dev.off()

  typec.evo <- reduceDimension(typec.evo, max_components = 2,
      method = 'DDRTree')

  typec.evo <- orderCells(typec.evo)

  plot_cell_trajectory(typec.evo, color_by = "epitype")
  dev.off()

  ### Gap exists between normal type C cells and tumor type C cells
  ### try to fix gap with indolent/aggressive samples

  # lv.1=sort(unique(normale@meta.data$epitype))
  lv.1=sort(unique(epi.indagg@meta.data$epitype))
  lv.2=sort(unique(epi13@meta.data$epitype))
  # lv.3=sort(unique(epi.crpc@meta.data$epitype))
  lv=list()
  for(i in 1:2)
  {
    lv[[i]]=get(paste("lv", i, sep="."))
  }
  shortest.2=dis.greedy(markerlist=markerlist.epi, lv=lv, layer=2)
  print(shortest)
  print(easy.shortest(shortest.2))
  # $inter.path
  #      etB1                    etC1             etC4.AGG etC5.AGG etI1.AGG
  # tB1  ""                      "tB1__etC1"      ""       ""       ""
  # tB2  ""                      "tB2__tB1__etC1" ""       ""       ""
  # tB4  "tB4__etB1"             ""               ""       ""       ""
  # tBC3 "tBC3__etB1"            ""               ""       ""       ""
  ### the cluster etC1 is the origin

  typec=subset( epi, subset=(epitype%in%c("tB1", "tB2", "tBC3", "C4", "etC1")) )
  table(typec@meta.data$epitype)
  typec@meta.data$epitype[typec@meta.data$epitype=="C4"]="NormC"
  typec@meta.data$epitype[typec@meta.data$epitype%in%c("tB1", "tB2", "tBC3")]="TumorC"
  typec@meta.data$epitype[typec@meta.data$epitype=="etC1"]="IntC"
  gene_metadata=data.frame(gene_short_name=rownames(typec@assays$SCT@data))
  rownames(gene_metadata)=rownames(typec@assays$SCT@data)
  typec.evo <- newCellDataSet(  typec@assays$SCT@data,
                            phenoData = new("AnnotatedDataFrame", typec@meta.data),
                            featureData = new("AnnotatedDataFrame", gene_metadata), 
                            expressionFamily=negbinomial.size() )
  table(typec.evo@phenoData@data$epitype)
  table(typec.evo@phenoData@data$pheno)
  typec.evo <- estimateSizeFactors(typec.evo)
  typec.evo <- estimateDispersions(typec.evo)

  typec.evo <- detectGenes(typec.evo, min_expr = 0.1)
  ### only keep expressed genes
  expressed_genes <- row.names(typec.evo)[typec.evo@featureData@data$num_cells_expressed>= 10]
  typec.evo <- typec.evo[expressed_genes,]

  ### use all significant markers of clusters as ordering genes
  tumor.markers.ori=markerlist.epi$"C4__tBC3"
  tumor.markers=tumor.markers.ori
  tumor.markers$pct.1=tumor.markers.ori$pct.2
  tumor.markers$pct.2=tumor.markers.ori$pct.1
  tumor.markers$avg_logFC=-1*tumor.markers$avg_logFC

  markers=tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[abs(markers$avg_logFC)>log(1.5), ]
  markers$foldChange=exp(markers$avg_logFC)
  # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
  markers=markers[order(markers$foldChange, decreasing=TRUE), ]
  # ordering.genes=unique(markers$gene)
  ordering.genes=unique(rownames(markers))

  typec.evo <- setOrderingFilter(typec.evo,  ordering.genes) # 
  pdf("plot_ordering_genes.pdf")
  plot_ordering_genes(typec.evo)
  dev.off()

  typec.evo <- reduceDimension(typec.evo, max_components = 2,
      method = 'DDRTree')
  typec.evo <- orderCells(typec.evo)
  plot_cell_trajectory(typec.evo, color_by = "epitype")
  dev.off()

  # set a indicator for time, and sort again
  typec.evo@phenoData@data$timeIset=rep(1, ncol(typec.evo))
  typec.evo@phenoData@data$timeIset[typec.evo@phenoData@data$phenocluster=="typec.evo.7"]=0
  typec.evo <- orderCells(typec.evo, root_state = typec.evo@phenoData@data$timeIset)

  # eFig.D.1
  pdf("plot_cell_trajectory_byPseudotime_typec.evo.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_typec.evo.png")
  plot_cell_trajectory(typec.evo, color_by = "Pseudotime", cell_size=0.4)
  dev.off()

  # eFig.D.2
  pdf("plot_cell_trajectory_fig.cluster_typec.evo.pdf", width=5, height=5)
  # png("plot_cell_trajectory_byPseudotime_typec.evo.png")
  plot_cell_trajectory(typec.evo, color_by = "fig.cluster", cell_size=0.4)
  dev.off()

  typec.evo.time=typec.evo@phenoData@data$Pseudotime
  names(typec.evo.time)=rownames(typec.evo@phenoData@data)

  png("plot_cell_trajectory_allIn1_typec.evo.png")
  plot_cell_trajectory(typec.evo, color_by = "phenocluster")
  dev.off()

  png("plot_cell_trajectory_details_phenocluster_typec.evo.png", width=30*100, height=6*100)
  typec=plot_cell_trajectory(typec.evo, color_by = "phenocluster") +
      facet_wrap(~fig.cluster, nrow = 1)
  print(typec)
  dev.off()

  # eFig.D.3
  png("plot_cell_trajectory_details_cluster_typec.evo.png", width=30*100, height=5*100)
  typec=plot_cell_trajectory(typec.evo, color_by = "fig.cluster", cell_size = 0.8) +
      facet_wrap(~fig.cluster, nrow = 1)
  print(typec)
  dev.off()
  pdf("plot_cell_trajectory_details_cluster_typec.evo.pdf", width=30, height=5)
  typec=plot_cell_trajectory(typec.evo, color_by = "fig.cluster", cell_size = 0.4) + facet_wrap(~fig.cluster, nrow = 1)
  print(typec)
  dev.off()

  saveRDS(typec.evo, "monocle.typec.evo.rds")
  typec.evo=readRDS("monocle.typec.evo.rds")

  # find DEGs along Pseudotime 
  diff_test_res.typec.evo <- differentialGeneTest(typec.evo, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  # saveRDS(diff_test_res.typec.evo, "diff_test_res.typec.evo.rds")
  diff_test_res.typec.evo=readRDS("diff_test_res.typec.evo.rds")
  # diff_test_res.bak=diff_test_res
  diff_test_res =  diff_test_res.typec.evo
  str(diff_test_res[,c("gene_short_name", "pval", "qval")])
  diff_test_res=subset(diff_test_res, qval < 0.1 & pval < 0.01)
  sig_gene_names=intersect(diff_test_res[,"gene_short_name"], subtype_markers$"TF_fromLieBing")
  imp.genes=c()
  imp.genes=c(imp.genes, c("BMPR2", "CCND2", "SOX4"))
  imp.genes=c(imp.genes, c("CDH2", "TWIST2", "HHIP"))
  imp.genes=c(imp.genes, c("COL5A1", "ITGA2", "MATN2", "SEMA3A"))
  imp.genes=c(imp.genes, sig_gene_names)
  imp.gen es=intersect(rownames(typec.evo), imp.genes)
  pdf("heatmap.trajectory.typec.evo.pdf")
  plot_pseudotime_heatmap(typec.evo[imp.genes,],
                  num_clusters = 3,
                  cores = 1,
                  show_rownames = T)
  dev.off()

  pdf("exp.along.time.typec.evo.pdf")
  imp.genes=c()
  imp.genes=c(imp.genes, c("BMPR2", "CCND2", "SOX4"))
  imp.genes=c(imp.genes, c("CDH2",  "HHIP"))# "TWIST2",
  imp.genes=c(imp.genes, c("COL5A1", "ITGA2", "MATN2", "SEMA3A"))
  genes <- imp.genes
  plot_genes_branched_pseudotime(typec.evo[genes,],
                       branch_point = 1,
                       color_by = "fig.cluster",
                       ncol = 1)
  dev.off()
}


### NMF new (library "NNMF")
{
  require(Seurat)
  # epi.tumor.bak=readRDS("epi.tumor.highCNV.harmony.adjed.rds")
  epi.tumor=epi.tumor.bak

  ### NNMF
  # refer to https://mran.microsoft.com/snapshot/2017-01-23/web/packages/NNLM/vignettes/Fast-And-Versatile-NMF.html
  {
    require(NNLM)
    epi.tumor=SetIdent(epi.tumor, value = "fig.pheno")
    table(epi.tumor@meta.data$fig.pheno)
    # Normal prostate   Indolent incidental PCa Aggressive incidental PCa
    #               0                      1786                      2437
    #     Primary PCa                     mHSPC                  post-ADT
    #            6449                      6439                       176
    #            CRPC                      NEPC
    #            6250                       169
    temp.levels=levels(epi.tumor@meta.data$fig.pheno)
    # delete "Normal prostate"
    temp.levels=temp.levels[-1]
    epi.tumor@meta.data$fig.pheno=as.character(epi.tumor@meta.data$fig.pheno)
    epi.tumor@meta.data$fig.pheno=factor(epi.tumor@meta.data$fig.pheno, levels=temp.levels)
    table(epi.tumor@meta.data$fig.pheno)
    ### get variable genes at first, consider var.genes across phenotypes and var.genes in each phenotype 
    ## select second small phenotype's sample number to balance classes
    var.genes=c()
    ## class-balanced samples 
    epi.tumor.sub=subset(epi.tumor, downsample=1786)
    ## consider var.genes across phenotypes 
    epi.tumor.sub=SCTransform(epi.tumor.sub, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = TRUE, do.scale=TRUE, method="glmGamPoi", do.center = TRUE, variable.features.n = 1000)
    var.genes=union(var.genes, epi.tumor.sub@assays$SCT@var.features)
    str(var.genes)
    ## consider var.genes in each phenotype 
    for(i in as.character(unique(epi.tumor@meta.data$fig.sample)))
    {
      epi.tumor.sub=subset(epi.tumor, fig.sample==i)
      epi.tumor.sub=SCTransform(epi.tumor.sub, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = TRUE, do.scale=TRUE, method="glmGamPoi", do.center = TRUE, variable.features.n = 500)
      var.genes=union(var.genes, epi.tumor.sub@assays$SCT@var.features)
    }
    str(var.genes)
    # chr [1:5245] "MSMB" "NPY" "PLA2G2A" "OLFM4" "LCN2" "LTF" "ACPP" "MMP7" ...

    epi.tumor.sub=subset(epi.tumor, downsample=1786)
    epi.tumor.sub=SCTransform(epi.tumor.sub, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi", do.center = TRUE)
    table(var.genes%in%rownames(epi.tumor.sub@assays$SCT@scale.data))
    # FALSE  TRUE
    #     7  5238
    # var.gene.bak=var.genes
    # var.genes=var.gene.bak
    temp.var.genes=var.genes[var.genes%in%rownames(epi.tumor.sub@assays$SCT@scale.data)]
    epi.tumor.sub@assays$SCT@scale.data=epi.tumor.sub@assays$SCT@scale.data[temp.var.genes, ]   
    epi.tumor.sub@assays$SCT@scale.data[epi.tumor.sub@assays$SCT@scale.data<0]=0
    decomp=list()
    set.seed(100)
    decomp[["all"]] <- nnmf(epi.tumor.sub@assays$SCT@scale.data, 50, rel.tol = 1e-5)
    ### NMF for each sample
    for(i in as.character(unique(epi.tumor@meta.data$fig.sample)))
    {
      epi.tumor.sub=subset(epi.tumor, fig.sample==i)
      epi.tumor.sub=SCTransform(epi.tumor.sub, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE, do.scale=TRUE, method="glmGamPoi", do.center = TRUE, variable.features.n = 1000)
      temp.var.genes=var.genes[var.genes%in%rownames(epi.tumor.sub@assays$SCT@scale.data)]
      epi.tumor.sub@assays$SCT@scale.data=epi.tumor.sub@assays$SCT@scale.data[temp.var.genes, ]
      epi.tumor.sub@assays$SCT@scale.data[epi.tumor.sub@assays$SCT@scale.data<0]=0
      table(epi.tumor.sub@meta.data$fig.sample)
      set.seed(100)
      decomp[[i]]=nnmf(epi.tumor.sub@assays$SCT@scale.data, 20, rel.tol = 1e-5)
    }
    # decomp[[1]]
    str(decomp)
    saveRDS(decomp, "decomp.epi.tumor.v2.rds")


    ### combine all W matrices, and do matrix decomposition again to get H matrix 
    ## combine 
    epi.tumor=epi.tumor.bak
    epi.tumor@assays$SCT@scale.data[epi.tumor@assays$SCT@scale.data<0]=0
    temp.var.genes=var.genes[var.genes%in%rownames(epi.tumor@assays$SCT@scale.data)]
    epi.tumor@assays$SCT@scale.data=epi.tumor@assays$SCT@scale.data[temp.var.genes, ]
    # initialize 0 matrix for combined W
    wcombine=matrix(0, nrow=length(temp.var.genes), ncol=0)
    rownames(wcombine)=temp.var.genes
    for(i in names(decomp))
    {
      toadd=decomp[[i]]$W
      temp=matrix(0, nrow=length(temp.var.genes), ncol=ncol(toadd))
      rownames(temp)=temp.var.genes
      temp[rownames(toadd), ]=toadd
      wcombine=cbind(wcombine, temp)
    }
    ## matrix decomposition
    decomp.final = nnlm(wcombine, epi.tumor@assays$SCT@scale.data);
    # check similarity between origin expression matrix and W*H
    temp.exp.m=wcombine%*%decomp.final$coefficients
    temp.exp.m.2=epi.tumor@assays$SCT@scale.data
    table(temp.exp.m.2>=0)
    table(colnames(temp.exp.m.2)==colnames(temp.exp.m))
    temp.exp.m[1:4, 1:4]
    temp.exp.m.2[1:4, 1:4]
    cor(temp.exp.m[1, ], temp.exp.m.2[1, ])
    cor(temp.exp.m[, 1], temp.exp.m.2[, 1])

    ### Create a DimReduc object
    ### refer to https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/CreateDimReducObject
    str(epi.tumor@reductions$PCA)
    m.embeddings=t(decomp.final$coefficients)
    m.loadings=wcombine
    colnames.m=paste("nmf", 1:ncol(m.embeddings), sep="_")
    colnames(m.embeddings)=colnames.m
    colnames(m.loadings)=colnames.m
    temp=CreateDimReducObject(
      embeddings = m.embeddings,
      loadings = m.loadings,
      projected = matrix(0, 0, 0),
      assay = "SCT",
      stdev = numeric(),
      key = "nmf_",
      global = FALSE,
      jackstraw = NULL,
      misc = list()
    )
    epi.tumor=epi.tumor.bak
    epi.tumor@reductions$nmf=temp

    # saveRDS(epi.tumor, "epi.tumor.nmf.v2.rds")
    # epi.tumor=readRDS("epi.tumor.nmf.v2.rds")
  }

  ### clustering programs
  {
    str(epi.tumor@reductions$nmf)
    table(epi.tumor@reductions$nmf@cell.embeddings>0)
    table(epi.tumor@reductions$nmf@cell.embeddings>=0)
    str(epi.tumor@reductions$nmf@cell.embeddings)

    ### clustering of nmf components
    {
      # get distance between programs
      # refer to https://www.rdocumentation.org/packages/amap/versions/0.8-18/topics/Dist
      require(amap)
      # scale every iNMF program's weights, and get distance between pair of programs
      {
        ### clustring with mst.knn
        # refer to https://cran.r-project.org/web/packages/mstknnclust/vignettes/guide.html
        # pdist=Dist(t(scale(epi.tumor@reductions$nmf@feature.loadings)), method = "correlation", nbproc = 2, diag = TRUE, upper = TRUE)
        # pdist.m = as.matrix(pdist)
        # pcluster = mst.knn(pdist.m)
      }
      # pdist=dist(t(scale(epi.tumor@reductions$nmf@feature.loadings)), method = "correlation")
      pdist=Dist(t((epi.tumor@reductions$nmf@feature.loadings)), method = "euclidean", nbproc = 2, diag = TRUE, upper = TRUE)
      # do hierarchical clustering
      dend = hclust(pdist) # hierarchical clustering
      # refer to: https://towardsdatascience.com/10-tips-for-choosing-the-optimal-number-of-clusters-277e93d72d92
      # refer to: https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
      # refer to: https://uc-r.github.io/hc_clustering
      require("factoextra")
      toprint1=fviz_nbclust(t((epi.tumor@reductions$nmf@feature.loadings)), FUN = hcut, method = "silhouette", k.max=100)
      toprint2=toprint1+theme(text = element_text(size=5), axis.text.x = element_text(angle = 45))
      print(toprint2)
      dev.off()

      require(ComplexHeatmap)
      # require(dendextend)
      # https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
      pdf("heatmap.distance.pdf")
      # pdist=Dist(t(scale(epi.tumor@reductions$nmf@feature.loadings)), method = "euclidean", nbproc = 2, diag = TRUE, upper = TRUE)
      # dend = hclust(pdist) # hierarchical clustering
      gdend.cut <- cutree(dend , k = 32)
      # row_dend = hclust((pdist))
      temp = Heatmap(as.matrix(pdist), row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5), col=c("red", "white"),  column_split = gdend.cut,  row_split = gdend.cut, row_dend_gp = gpar(fontsize = 5), column_dend_gp = gpar(fontsize = 5))  # split =gdend.cut, km=7,  column_split = gdend.cut, cluster_rows = color_branches(row_dend, k = 7), cluster_columns = color_branches(row_dend, k = 7)
      print(temp)
      dev.off()
    }


    ### average top200 genes' weight for program cluster ("cweight")
    {
      pweight=(t(epi.tumor@reductions$nmf@feature.loadings))
      str(pweight)
      # num [1:630, 1:5245] 0 0 3.57 0 0 ...
      current.cluster.ids <- names(gdend.cut)
      new.cluster.ids <- sprintf("P_%02d", gdend.cut)
      pweight.cluster=plyr::mapvalues(x = rownames(pweight), from = current.cluster.ids, to = new.cluster.ids)
      table(pweight.cluster)
      table(gdend.cut)

      pweight[1:5, 1:5]
      str(pweight)
      str(pweight.cluster)
      ### make non-top200 genes' weight = 0
      !!!
      cweight=aggregate(pweight, by = list(pweight.cluster), FUN = mean, na.rm = TRUE)
      rownames(cweight)=cweight$Group.1
      cweight=cweight[, colnames(cweight)!="Group.1"]
      cweight=as.matrix(cweight)
      table(colnames(cweight)%in%rownames(epi.tumor@assays$SCT@scale.data))
    }


    ### get program's top 100 genes
    {
      ### get program cluster's gene loading matrix (cweight)
      ## get each gene's priority score in each program
      ## refer to: https://doi.org/10.1016/j.cell.2021.08.003
      cid=sort(unique(pweight.cluster))
      cweight=matrix(0, nrow=length(cid), ncol=ncol(pweight))
      rownames(cweight)=cid
      colnames(cweight)=colnames(pweight)
      for(i in cid)
      {
        wi=pweight[pweight.cluster==i, ]
        if(is.matrix(wi))
        {wimean=colMeans(wi)}else
        {wimean=wi}
        uniratio=rep(Inf, ncol(pweight))
        names(uniratio)=colnames(pweight)
        for(j in cid)
        {
          if(j!=i)
          {
            wj=pweight[pweight.cluster==j, ]
            if(is.matrix(wj))
            {wjmean=colMeans(wj)}else
            {wjmean=wj}
            uniratio.temp=wimean*log( (wimean+1)/(wjmean+1) )
            temp=data.frame(uniratio.temp=uniratio.temp, uniratio=uniratio)
            uniratio=apply(temp, 1, min)
          }
        }
        summary(uniratio)
        cweight[i, ]=uniratio
      }
      ### get top100 genes for each program cluster
      topn=100
      topgenes=list()
      for(i in cid)
      { 
        temp=sort(cweight[i, ], decreasing=TRUE)
        topgenes[[i]]=names(temp[1:topn])
      }
      ### check whether there are overlapped top genes between program cluters
      u.topgenes=Reduce(union, topgenes)
      str(u.topgenes)
      # chr [1:2668] "OGN" "PTGFR" "SLITRK5" "CHRD" "OMD" "EDNRB" "NOS1" "DNAH8" ...

      # only keep top 100 genes' weights
      temp.cweight=cweight
      for(i in cid)
      { 
        temp.cweight[i, !colnames(temp.cweight)%in%topgenes[[i]]]=0
      }
      # remove all-zero genes
      temp.sum=colSums(temp.cweight)
      temp.cweight=temp.cweight[, temp.sum>0]
    }

    ### give gene-programCluster weight to seurat obj's @reduction
    {
      ### Create a DimReduc object
      ### refer to https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/CreateDimReducObject
      str(epi.tumor@reductions$pca)
      ## matrix decomposition
      decomp.final = nnlm(t(temp.cweight), epi.tumor@assays$SCT@scale.data[colnames(temp.cweight),])
      # check similarity between origin expression matrix and W*H
      temp.exp.m=t(temp.cweight)%*%decomp.final$coefficients
      temp.exp.m.2=epi.tumor@assays$SCT@scale.data[colnames(temp.cweight),]
      # table(temp.exp.m.2>=0)
      table(colnames(temp.exp.m.2)==colnames(temp.exp.m))
      temp.exp.m[1:4, 1:4]
      temp.exp.m.2[1:4, 1:4]
      cor(temp.exp.m[1, ], temp.exp.m.2[1, ])
      cor(temp.exp.m[, 1], temp.exp.m.2[, 1])
      m.loadings=t(temp.cweight)
      m.embeddings=t(decomp.final$coefficients)
      table(colnames(m.embeddings)==colnames(m.loadings))
      temp=CreateDimReducObject(
        embeddings = m.embeddings,
        loadings = m.loadings,
        projected = matrix(0, 0, 0),
        assay = "SCT",
        stdev = numeric(),
        key = "P_",
        global = FALSE,
        jackstraw = NULL,
        misc = list() )
      epi.tumor@assays=epi.tumor.bak@assays
      epi.tumor@reductions$nmfp=temp
    }

    saveRDS(epi.tumor, "epi.tumor.nmf.32c.rds")
    # epi.tumor=readRDS("epi.tumor.nmf.32c.rds")
  }


  # featureplot NMF programs
  png("temp.featureplot.png", height=10*100, width=14*100)
  FeaturePlot(epi.tumor, features = sprintf("P_%02d", 1:14), ncol = 4)
  dev.off()

  ### dotplot NMF programs
  png("temp.dotplot.phenotype.png", height=4*100, width=10*100)
  DotPlot(object = epi.tumor, features = sprintf("P_%02d", 1:46), group.by="fig.pheno")+RotatedAxis()
  dev.off()

  ### dotplot NMF programs
  png("temp.dotplot.celltype.png", height=4*100, width=10*100)
  DotPlot(object = epi.tumor, features = sprintf("P_%02d", 1:46), group.by="fig.celltype")+RotatedAxis()
  dev.off()

  ### choose clusters' top genes to heatmap
  nmf.weight=epi.tumor@reductions$nmfp@feature.loadings
  nmf.genes=list()
  nmf.gset=c()
  pid=c(1:32)
  # pid=sprintf("P_%02d", pid)
  for(i in colnames(nmf.weight)[pid])
  {
    temp=sort(nmf.weight[, i], decreasing=TRUE)
    nmf.genes[[i]]=names(temp)[1:100]
    nmf.gset=union(nmf.gset, names(temp)[1:10])
  }

  ### heatmap topN DEGs for every NMF feature
  topn <- nmf.gset
  # topn$cluster=factor(topn$cluster, levels=c("Clike", "TypeC", "Luminal", "Basal", "CellCycle"))
  # topn=topn[order(topn$cluster), ]
  temp=SetIdent(epi.tumor, value = "fig.pheno")
  temp=subset(temp, downsample=300)
  # levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
  # tocheck=c(topn$gene, "KLK3", "AR")
  temp.heatmap<-DoHeatmap( temp, features = topn, slot="scale.data", assay="SCT", angle = 0 , size=2)
  # , subset=(fig.cluster%in%temp.clusters)
  pdf(paste("heatmap.fib.cluster.topn.pdf", sep=""),  width=(length(topn)/4+5), height=length(topn)/5)
  print(temp.heatmap)
  dev.off()

  setwd("./epi.tumor.nmf.res/")

  # for specific program, do heatmap and enrichment
  for(i in pid)#
  {
    res=try(
    {
      ii=sprintf("P_%02d", i)
      topn=nmf.genes[[ii]][1:20]
      write.table(nmf.genes[[ii]], file=paste(i, "txt", sep="."), col.name=F, row.name=F, quote=F, sep='\t')

      temp=SetIdent(epi.tumor, value = "fig.pheno")
      temp=subset(temp, downsample=300)
      # levels(temp)=c("Clike", "TypeC", "Luminal", "Basal")
      # tocheck=c(topn$gene, "KLK3", "AR")
      temp.heatmap<-DoHeatmap( temp, features = topn, slot="scale.data", assay="SCT", angle = 0 , size=2)
      # , subset=(fig.cluster%in%temp.clusters)
      pdf(paste( i, ".heatmap.top50.NMF.", ".pdf", sep="" ),  width=(length(topn)/4+5), height=length(topn)/8+2)
      print(temp.heatmap)
      dev.off()

      enrichment.gene.set.plot(nmf.genes[[ii]], topNterm=10, prefix=paste(i, ".", sep=""))
    })
    if(inherits(res, "try-error"))
    {
      #error handling code, maybe just skip this iteration using
      next
    }
  }

  setwd("..")

  ### correlation between programs
  {
    cweight=aggregate(temp.exp.m.2=epi.tumor@assays$SCT@scale.data[colnames(temp.cweight),], by = list(pweight.cluster), FUN = mean, na.rm = TRUE)
  }

}


