require(svglite)
require(dplyr)
require(grid)
# utility functions used for plotting

save_pheatmap <- function(x, filename, width=8, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  svg(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


drop_genecol <- function(df){
  df <- as.data.frame(df)
  stopifnot("gene" %in% colnames(df))
  genes <- df$gene
  df <- df %>% select(-gene)
  rownames(df) <- genes
  return(df)
}

# clustering based on correlation as implemented by MaSigPro
clusterfun <- function(clusterdata, k){
  dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                    nrow(clusterdata), nrow(clusterdata)) - 
    cor(t(clusterdata), use = "pairwise.complete.obs")
  clust <- hclust(as.dist(dcorrel), method = "ward.D")
  
  cut <- cutree(clust, k = k)
  
  OUTPUT <- list(cut, clust$order)
  names(OUTPUT) <- c("cut", "order")
  OUTPUT
}