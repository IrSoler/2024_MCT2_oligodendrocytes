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
library(ggplot2)

# ----------------------------------------------------------- #
# --- 2. Load data ------------------------------------------ #
# ----------------------------------------------------------- #

data.oligos = readRDS("./data_oligos.rds")
data.oligos.av = readRDS("./data_oligos_av.rds")

# ----------------------------------------------------------- #
# --- 3. Barplot number of cells ---------------------------- #
# ----------------------------------------------------------- #

path = "./"

png(paste0(path, "number_cells.png"), width = 6400, height = 5830, res = 900)

ggplot(data = p, aes(x=Var1, y=Freq, fill=Var2)) + geom_col(position = "dodge2") + ylim(c(0,3000)) +
  theme_classic() + xlab("Clusters") + ylab("Number of cells") + 
  theme(text = element_text(size = 14, colour = "black"),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = c("Control" = "#89b0ae", "MS" = "#eaac8b"), name = "Condition")

dev.off() 


# ----------------------------------------------------------- #
# --- 3. Boxplots ------------------------------------------- #
# ----------------------------------------------------------- #

path = "./"

png(paste0(path, "clusters.png"), width = 6400, height = 5830, res = 900)

VlnPlot(data.oligos.av, features = features, slot = "data", alpha = 0, pt.size = 0, assay = "RNA", y.max = 1.3) +
  geom_boxplot() + ylab(expression("log"[10]*" expression level")) +  ggtitle(NULL) + xlab("Clusters") +
  geom_jitter(size = 0.5, colour = "#6c757d") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_signif(comparisons=list(c("OL-SLC5A11", "OL-SGCZ")), annotations="*",
              y_position = 1.2, tip_length = 0.01, vjust=0.3) +
  geom_signif(comparisons=list(c("OL-SLC5A11", "OL-LINC01608")), annotations="*",
              y_position = 1, tip_length = 0.01, vjust=0.3) +
  geom_signif(comparisons=list(c("OL-SLC5A11", "OL-HSPA1A")), annotations="*",
              y_position = 1.1, tip_length = 0.01, vjust=0.3) +
  scale_fill_manual(values = c("OL-SLC5A11" = "#F8CC67", "OL-LINC01608" = "#E78957", "OL-HSPA1A" = "#B48AC7", "OL-SGCZ" = "#4C7395"))

dev.off()