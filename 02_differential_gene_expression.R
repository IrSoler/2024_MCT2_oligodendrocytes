# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# Manuscript: Monocarboxylate transporter 2 is required for the maintenance of myelin and axonal 
# integrity by oligodendrocytes
#
# Code author: Irene Soler Sáez (isoler@cipf.es)
#
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

# Created by: Irene Soler Sáez

# ----------------------------------------------------------- #
# --- 1. Load packages -------------------------------------- #
# ----------------------------------------------------------- #

library(Seurat)
library(DESeq2)

# ----------------------------------------------------------- #
# --- 2. Load data ------------------------------------------ #
# ----------------------------------------------------------- #

data.oligos = readRDS("./data_oligos.rds")

# ----------------------------------------------------------- #
# --- 3. Pseudobulk ----------------------------------------- #
# ----------------------------------------------------------- #

## Pseudobulk aggregation
data.oligos.av = AggregateExpression(data.oligos, return.seurat = TRUE, group.by = c('variable2test','sample_id'), assays = "RNA")
Idents(data.oligos.av) = factor(data.oligos.av@meta.data$class)

## Parsing names & variables
p = colnames(data.oligos.av)
p = strsplit(p, "_")
p1= unlist(lapply(p, function(x){x[[2]]}))

data.oligos@meta.data$sample_id = gsub("_", "-", data.oligos@meta.data$sample_id)

v = c()

for (p_sample in p1){
  a = unique(data.oligos@meta.data[data.oligos@meta.data$sample_id==p_sample, "variable2test"])
  v = c(v, a)
}

data.oligos.av@meta.data$age = v
data.oligos.av@meta.data$condition = Idents(data.oligos.av)

data.oligos.av$condition = factor(data.oligos.av$condition)

### Adapt to the variable to test


# ----------------------------------------------------------- #
# --- 3. Differential gene expression ----------------------- #
# ----------------------------------------------------------- #

## Deseq2 object
dds = DESeqDataSetFromMatrix(data.oligos.av[["RNA"]]$counts, 
                             colData = data.oligos.av@meta.data, 
                             design = ~ 0 + variable2test + covariables2adjust)

## Prefiltering
smallestGroupSize = 2
keep = rowSums(counts(dds) >= 1) >= smallestGroupSize
dds = dds[keep,]

## DGE
dds = DESeq(dds)

## Extract results DGE
dif1 = list(c("group1"), c("group2"))

p = results(dds, contrast=dif1, pAdjustMethod = "BH")

p = as.data.frame(p)

p$genes = rownames(p)
p$upregulation = ifelse(p$log2FoldChange > 0, "group1", "group2")
