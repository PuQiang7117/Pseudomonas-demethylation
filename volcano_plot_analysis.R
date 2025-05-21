if (!require("edgeR")) install.packages("BiocManager"); BiocManager::install("edgeR")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("dplyr")) install.packages("dplyr")

library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)

otu_counts <- read.csv("otu_table.csv", row.names = 1, check.names = FALSE)

metadata <- read.csv("metadata.csv")

otu_counts <- otu_counts[, metadata$SampleID]

group <- factor(metadata$Group)
dge <- DGEList(counts = otu_counts, group = group)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)

res <- topTags(lrt, n = Inf)$table
res$OTU <- rownames(res)

res$group <- "NS" 

res$group[res$FDR >= 0.05 & abs(res$logFC) < 2] <- "down"

res$group[res$FDR < 0.05 & abs(res$logFC) >= 2] <- "up"

ggplot(res, aes(x = logFC, y = -log10(FDR), color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("down" = "cyan", "NS" = "gray", "up" = "salmon")) +
  labs(title = "Volcano Plot (edgeR)", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = ifelse(FDR < 0.05 & abs(logFC) >= 2, OTU, "")), size = 3)

write.csv(res, "differential_OTUs_edgeR.csv", row.names = FALSE)
ggsave("volcano_plot_edgeR.png", width = 8, height = 6)
