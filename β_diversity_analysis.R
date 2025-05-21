if (!require("vegan")) install.packages("vegan", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("ggrepel")) install.packages("ggrepel", dependencies = TRUE)
if (!require("ellipse")) install.packages("ellipse", dependencies = TRUE)
library(vegan)
library(ggplot2)
library(ggrepel)
library(ellipse)

otu_table <- read.csv("otu_table.csv", row.names = 1, check.names = FALSE)

otu_rel <- decostand(otu_table, method = "total")

bray_dist <- vegdist(otu_rel, method = "bray")

nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100)

stressplot(nmds_result)

nmds_points <- as.data.frame(nmds_result$points)
nmds_points$SampleID <- rownames(nmds_points)

metadata <- read.csv("metadata.csv") 
nmds_points <- merge(nmds_points, metadata, by = "SampleID")

ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = Group, shape = Condition)) +
  geom_point(size = 3, alpha = 0.7) + 
  stat_ellipse(aes(color = Group), level = 0.95, linetype = 1, size = 1) + 
  scale_shape_manual(values = c(16, 17)) + 
  theme_minimal() +
  labs(title = "nMDS",
       x = "nMDS1", y = "nMDS2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right")

ggsave("nmds_plot.png", width = 6, height = 5)
