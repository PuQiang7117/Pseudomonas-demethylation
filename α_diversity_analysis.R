if (!require("vegan")) install.packages("vegan", dependencies = TRUE)
library(vegan)

otu_table <- read.csv("otu_table.csv", row.names = 1, check.names = FALSE)

otu_rel <- decostand(otu_table, method = "total")

simpson <- diversity(otu_rel, index = "simpson")

simpson_diversity <- data.frame(
  SampleID = rownames(otu_rel),
  Simpson = simpson
)

print(simpson_diversity)