require(ggplot2)
require(dplyr)
require(tidyr)
require(readr)
require(pheatmap)
require(RColorBrewer)

heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(100)
heatmap_colors2 <- colorRampPalette(colors = c("blue", "white", "red"))(100)


prep_data <- function(df, trafo= "norm_day0"){
  if(trafo == "norm_day0"){
    mycol <- "avg_norm"
  } else if(trafo == "rtm"){
    mycol <- "avg_norm2"
  } else if(trafo == "rtm2"){
    mycol <- "avg_norm3"
  }
  mycols <- c("time", mycol, "module")
  # remove duplicates and make wide --> prepare for pheatmap
  df <- df[mycols]
  df <- unique(df)
  df <- df %>% pivot_wider(names_from = time, values_from = !!as.symbol(mycol))
  df <- as.data.frame(df)
  rownames(df) <- df$module
  df <- df[,2:ncol(df)]
  
  return(df)
}


plot_heatmap <- function(df, sname, trafo){
  cellwidth = 8
  cellheight = 8
  
  if(trafo == "norm_day0"){
    sdir <- "figures/tcell_module/heatmap_normd0/"
    breaks <- seq(-1,1, length.out = 101)
  } else if(trafo == "rtm"){
    sdir <- "figures/tcell_module/heatmap_rtm/"
    breaks <- seq(-1,1, length.out = 101)
  }else if(trafo == "rtm2"){
    sdir <- "figures/tcell_module/heatmap_rtm2/"
    breaks <- seq(-1,1, length.out = 101)   
  }
  
  sname <- paste0(sdir, "heatmap_module_kinetics_", sname, ".png")
  
  pheatmap(df,
           breaks = breaks,
           color = heatmap_colors,
           scale = "none",
           cluster_rows = F,
           cluster_cols = F,
           cellwidth = cellwidth,
           cellheight = cellheight,
           filename = sname, 
           fontsize = 8,
           width = 4,
           height = 2)
}

df <- read.csv("data/data_summary/data_rtm_gene_module.csv")

mymodules <- c("expert_list")

df <- df[!(df$module %in% mymodules),]

out <- df %>% group_by(ID)


out_list <- out %>% group_split()
out_keys <- out %>% group_keys()
out_keys <- out_keys$ID


mytrafos <- c("norm_day0", "rtm", "rtm2")

for(trafo in mytrafos){
  out_list2 <- lapply(out_list, prep_data, trafo)
  mapply(plot_heatmap, out_list2, out_keys, trafo) 
}


