library(tidyverse)
library(pheatmap)
library(igraph)

# Read in the mash distance data, with the file directories cut out from the transcriptome names for brevity
mash_data <- read_delim("mash_out3.txt", col_names = c("referenceID", "queryID", "distance", "p_val", "shared_hashes"))

# Use igraph to make the table into a pairwise matrix of distances
mash_matrix <- as_adjacency_matrix(graph.data.frame(mash_data, directed=FALSE), names=TRUE, sparse=FALSE, attr="distance", type='both')

# Make a heatmap, remove the dendrogram for rows by setting height to 0, since the same dendrogram is present for columns
mash_heatmap <- pheatmap(mash_matrix, treeheight_row = 0)

# A function to save the heatmap as a png file with 150 dpi. It ended up being low resolution so I exported it as pdf instead
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(mash_heatmap, "my_heatmap.png")
