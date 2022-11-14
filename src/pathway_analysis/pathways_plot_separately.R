# take multiple pathways and plot them together in heatmap
# store table as xlsx
# increase memory for java application
options(java.parameters = "- Xmx1024m")
require(readr)
require(tidyr)
require(ggplot2)
require(stringr)
require(dplyr)
require(forcats)
require(RColorBrewer)
require(pheatmap)
require(xlsx)
source("src/utils.R")

get_files <- function(category, use_gene_module){
  
  # list files
  if(use_gene_module){
    mypath <- "data/ORA_output_gene_module/"
  } else{
    mypath <- "data/ORA_output/"
  }
  
  filenames <- list.files(mypath, pattern = ".csv", full.names = T)
  filenames2 <- list.files(mypath, pattern = ".csv", full.names = F)
  
  filenames <- filenames[grepl(category, filenames)]
  filenames2 <- filenames2[grepl(category, filenames2)]
  
  filenames <- filenames[grepl("Nir|Peine|Proserpio", filenames)]
  filenames2 <- filenames2[grepl("Nir|Peine|Proserpio", filenames2)]
  
  # read files and combine
  files <- lapply(filenames, read.table, header = T)
  out <- mapply(myfun1, files, filenames2, SIMPLIFY = F)
  out <- bind_rows(out)
  
  filenames_xlsx <- str_sub(filenames2, 1, -5)
  
  # store as supplementary table, note that this is quite memory consuming
  store_table <- F
  if(store_table){
    print("storing table, this takes memory, set to false if not needed...")
    for(i in seq_along(files)){
      xlsx_file <- files[[i]]
      # reformat the name to make it short enough to fit on an xlsx sheet
      xlsx_name <- str_split(filenames_xlsx[[i]], pattern = "_")
      xlsx_name <- paste(xlsx_name[[1]][1], xlsx_name[[1]][3], xlsx_name[[1]][7], sep = "_")
      
      if(i==1){
        append <- FALSE
      } else {
        append <- TRUE
      }
      
      xlsx_file$ID <- NULL
      write.xlsx(xlsx_file, file = "tables/Burt_etal_Supplementary_Table_Pathway_Analysis.xlsx",
                   sheetName=xlsx_name, append=append, row.names = FALSE)
      }
  }

  return(out)
}

myfun1 <- function(df, fname){
  # add the filename as a column, transform filename first
  name <- str_split(fname, "_")
  
  name <- name[[1]][5]
  
  df$name <- fname
  return(df)
}

# remove nonsignif categories
apply_fdr_filter <- function(df, alpha, filter){
  df <- df %>% filter(qvalue <= alpha)
  counts <- df %>% group_by(ID) %>% count()
  myselec <- counts$ID[counts$n>=filter]
  df <- df %>% filter(ID %in% myselec)  
  return(df)
}

proc_files <- function(out, user_curated, savedir){
  # add database as column and add -log10 column
  out <- out %>% separate(ID, into = c("DB"), remove = F)
  out <- out %>% mutate(fdrlog10 = -log10(qvalue))
  
  # rename transcription factor targets because the annotation doesnt work in combined feature
  mycategories <- c("GOBP", "WP", "REACTOME", "HALLMARK")
  out$DB[!(out$DB %in% mycategories)] <- "TFT"
  
  # use manual annotation 
  fname <- paste0(savedir, "pathways_curated.csv")
  out <- proc_pathways_curated(out, user_curated, fname)
  
  return(out)  
}

proc_pathways_curated <- function(out, user_curated, fname){
  
  # check if manual curation exists 
  # if it exists and user_curated True, reduce out to manual curation subset and change pathway names and ordering according to file
  # otherwise if it exists and user_curated False, add IDs to full table
  # otherwise return original
  if(file.exists(fname)){
    sep = ";"
    print("current separator is")
    print(sep)
    pathways_curated <- read.csv(fname, sep = sep)
    stopifnot("New_ID" %in% colnames(pathways_curated))
    
    # show only pathways in manual curation file or alternatively add index of manual curation to full table    
    if(user_curated){
      # reorder according to provided manual order in csv file
      pathways_curated$myorder <- rownames(pathways_curated)
      
      out <- inner_join(out, pathways_curated)
      out <- out %>% mutate(New_ID = fct_reorder(New_ID, as.numeric(myorder), .desc = T), ID = New_ID)
    } 
    
  } else {
    print("no manual curation found, returning original")
  }  
  return(out)
}


plot_heatmap <- function(out, sname, mythres){
  
  #prepare data output format for pheatmap, make wider
  # not that this can cause troubles if there are not enough hits in both categories
  n_categories <- length(unique(out$name))
  stopifnot(n_categories > 1)
  
  test <- out %>% select(name, fdrlog10, ID) %>% pivot_wider(names_from = name, values_from = fdrlog10)
  
  # all NAs get 0 instead (log10(1) = 0)
  test[is.na(test)]= 0
  
  test <- as.data.frame(test)
  # make coarse grained heatmap only showing signif *,** and *** differences
  rownames(test) <- str_remove(test$ID, "REACTOME_|HALLMARK_|")
  colnames(test) <- str_remove(colnames(test), ".csv")
  colnames(test) <- str_remove(colnames(test), "_ORA_REACTOME")
  colnames(test) <- str_remove(colnames(test), "_invitro")
  
  test <- test[,2:ncol(test)]
  mynames <- colnames(test)
  
  annot <- as.data.frame(mynames)
  annot$Model <- "Expo_Other"
  annot$Model[grepl("gamma", annot$mynames)] <- "Gamma"
  rownames(annot) <- annot$mynames
  annot <- annot %>% select(Model)

  # 
  # keep all entries where there is at least one significant hit
  test <- test[rowSums(test >= -log10(mythres)) > 0,]
  thres1 <- -log10(0.1)
  thres2 <- -log10(0.05)
  thres3 <- -log10(0.01)
 
  test2 <- test
  test2[test<thres1] = 0
  test2[test>thres1] = 1
  test2[test>thres2] = 2
  test2[test>thres3] = 3

  # convert to log

  
  
  #annot <- annot[,colnames(annot) %in% colnames(test2)]
  
  cellwidth = 10
  cellheight = 8
  width = 3.5
  height = 2.5
  # get three distinct types of red
  mycolors = colorRampPalette(brewer.pal(n = 4, name ="Reds"))(4)  
  mycolors[1] <- "white"
  mybreaks <- seq(-log10(mythres),-log10(0.01), length.out = 101)
  
  ann_colors = list(Model = c(Expo_Other = "chartreuse3", Gamma = "deepskyblue3"))
  
  pheatmap(test2,
           color = mycolors,
           annotation_col = annot,
           annotation_colors = ann_colors,
           annotation_legend = T,
           cluster_rows = T,
           cluster_cols = T,
           fontsize_row = 8,
           fontsize_col = 8,
           cellwidth = cellwidth,
           cellheight = cellheight,
           width = width,
           height = height,
           legend = T,
           filename = sname)
  
  
  p <- pheatmap(test2,
                color = mycolors,
                annotation_col = annot,
                annotation_colors = ann_colors,
                annotation_legend = T,
                cluster_rows = T,
                cluster_cols = T,
                fontsize_row = 8,
                fontsize_col = 8,
                cellwidth = cellwidth,
                cellheight = cellheight,
                legend = T,
                filename = sname)
  
  sname2 <- str_sub(sname, 1, -5)
  sname2 <- paste0(sname2, ".svg")
  
  save_pheatmap(p, sname2, width = 10, height = 20)
}


pipeline <- function(alpha, filter, width, height, category, user_curated, sname, use_gene_module){
  
  if(use_gene_module){
    savedir <- "figures/pathway_results/tcell_universe/"
  } else {
    savedir <- "figures/pathway_results/full_universe/"
  }
  #savedir <- "figures/pathway_results/"
  # list files in directory, read them and combine into data frame
  out <- get_files(category, use_gene_module)
  
  # process df
  out <- proc_files(out, user_curated, savedir)
  
  # plot heatmap (before applying filter)
  if(user_curated){
    sname0 <- "curated"
  } else {
    sname0 <- "uncurated"
  }

  sname_heatmap <- paste0(savedir, category, "_", sname,"_", sname0, ".png")
  print(paste("saving to", sname_heatmap))
  
  # keep only signif categories
  #out <- apply_fdr_filter(out, alpha, filter)
  
  # filtering step is now applied here
  plot_heatmap(out, sname_heatmap, alpha)
  
  
  return(out)
}

###################################################################################################
###################################################################################################
###################################################################################################
# analysis starts here
use_gene_module <- F # list of t cell genes as baseline
myFDR <- "FDR01"
alpha <- 0.1
out1 <- pipeline(alpha = alpha, filter = 1, width = 20, height = 30, category = "REACTOME", user_curated = F, sname = myFDR, use_gene_module)
#out1 <- pipeline(alpha = alpha, filter = 1, width = 20, height = 30, category = "HALLMARK", user_curated = F, sname = myFDR, use_gene_module)
#out1 <- pipeline(alpha = alpha, filter = 1, width = 20, height = 30, category = "TFT", user_curated = F, sname = myFDR, use_gene_module)
#out1 <- pipeline(alpha = alpha, filter = 1, width = 20, height = 30, category = "GOBP", user_curated = F, sname = myFDR, use_gene_module)
#out1 <- pipeline(alpha = 1, filter = 1, width = 20, height = 30, category = "WIKIPATHWAYS", user_curated = F, sname = myFDR, use_gene_module)
#out1 <- pipeline(alpha = alpha, filter = 1, width = 20, height = 30, category = "IMMUNESIG", user_curated = F, sname = myFDR, use_gene_module)



