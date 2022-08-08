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
setwd("")

### read primary PCa samples
### process from raw counts

