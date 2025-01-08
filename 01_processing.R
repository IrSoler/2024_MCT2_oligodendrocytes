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

# ----------------------------------------------------------- #
# --- 1. Load packages -------------------------------------- #
# ----------------------------------------------------------- #

library(Seurat)

# ----------------------------------------------------------- #
# --- 2. Load data ------------------------------------------ #
# ----------------------------------------------------------- #

data = readRDS("seurat_integrated.rds")
class(data)

DefaultAssay(data) = "RNA"

# ----------------------------------------------------------- #
# --- 3. Quality control ------------------------------------ #
# ----------------------------------------------------------- #

## Percentage mito genes
data[["ratio.mito"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "ratio.mito"), ncol = 3)

## Quality control filtering
data = subset(data, subset = nFeature_RNA >= 250 & nCount_RNA >= 400 & ratio.mito < 5)


# ----------------------------------------------------------- #
# --- 3. Normalization -------------------------------------- #
# ----------------------------------------------------------- #

data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

# ----------------------------------------------------------- #
# --- 4. Subset oligodendrocytes ---------------------------- #
# ----------------------------------------------------------- #

# Select oligos
oligos.assign = readxl::read_excel("oligos_metadata.xlsx")
rownames(oligos.assign) = oligos.assign$cells

# Filter oligodendrocytes
data.oligos = data[,rownames(oligos.assign)]

# Generate a new variable with the classification 2 (clusters), establishing it as indents
data.oligos$clusters = oligos.assign$ol_subtype
data.oligos$clusters = factor(data.oligos$clusters, levels = c("OL-SLC5A11", "OL-LINC01608", "OL-HSPA1A", "OL-SGCZ"))
Idents(data.oligos) = "clusters"


