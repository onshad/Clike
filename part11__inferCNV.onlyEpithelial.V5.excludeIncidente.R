
options(error = function() traceback(2))

require("Seurat")
require("infercnv")

setwd("/ ... /inferCNV")

### prepare cell type annotation for cells
epi=readRDS("../ ... /epi.fullsamples.V3.rds")
epitype=epi@meta.data$epitype
names(epitype)=rownames(epi@meta.data)


### write matrix to harddisk
saveRDS(epi@assays$SCT@counts, "epitype.V5.sct.counts.matrix.rds")

celltype=data.frame(cellname=rownames(epi@meta.data), celltype=epi@meta.data$epitype)
write.table(celltype, file="epitype.V5.celltype.annotation.txt", col.name=F, row.name=F, quote=F, sep='\t')
###  nL1, nL4, nL5 are not treated as control
ref.celltype=c("nL1", "nB2")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="epitype.V5.sct.counts.matrix.rds",
                                    annotations_file="epitype.V5.celltype.annotation.txt",
                                    delim="\t",
                                    gene_order_file="pos_Homo_sapiens.GRCh38.Ensembl91.createdByPyScriptFromInferCNVwebsite.txt",
                                    ref_group_names=ref.celltype)

out_dir="epitype.V5.sct.counts.out"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             plot_steps=FALSE,
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads=20
                             )



