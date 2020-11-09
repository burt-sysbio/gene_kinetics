# identify kinetic genes in crawford data using masigpro
# load crawford CD4 Arm data naive, d6,8,15,30
#   Differential expression analysis with limma

require(dplyr)
require(stringr)
require(readr)
require(maSigPro)

pdata <- read.csv("output/crawford_data_features.csv")
df <- read.csv("output/crawford_data_median_threshold.csv")
rownames(df) <- df$name
df <- df[colnames(df) != "name"]

cols <- c("name", "time", "infection", "cell_type")
pdata <- pdata[cols]
# rename cells not a mistake that its cell and cells below
pdata$cell_type[pdata$cell_type == "Naive CD44Lo CD8+ T cells"] <- "H2-Db GP33-specific CD8+ T cells"
pdata$cell_type[pdata$cell_type == "Naive CD44Lo CD4+ T cell"] <- "H2-IAb GP66 specific CD4+ T cell"
pdata <- pdata %>% mutate(cell_type = str_extract(cell_type, "CD."))


# create design matrix for masigpro, day 6 has for 1 case only 3 replicates instead of 4
make_des_matrix <- function(df, cell, inf){
  
  df <- filter(df, cell_type == cell & (infection == inf | infection == "None")) 
  
  if(cell == "CD4" & inf == "LCMV-Arm"){
    n_d6 <- 3
  } else {
    n_d6 <- 4
  }
    
  rep1 <- rep(1, 4)
  rep2 <- rep(2, n_d6)
  rep3 <- rep(3:5, each = 4)
  time <- as.numeric(df$time)
  replicate <- c(rep1, rep2, rep3)
  group <- rep(1, nrow(df))
  
  out <- cbind(time, replicate, group)
  rownames(out) <- df$name
  return(out)
}

#design_CD4_arm <- make_des_matrix(pdata, "CD4", "LCMV-Arm")
#design_CD4_cl13 <- make_des_matrix(pdata, "CD4", "LCMV-Clone 13")
#design_CD8_arm <- make_des_matrix(pdata, "CD8", "LCMV-Arm")
#design_CD8_cl13 <- make_des_matrix(pdata, "CD8", "LCMV-Clone 13")

#df_cd4_arm <- df[rownames(design_CD4_arm)]

#design <- make.design.matrix(design_CD4_arm, degree = 2)

#fit <- p.vector(df_cd4_arm, design=design)
#t <- T.fit(fit, step.method = "backward", alfa = 0.05)

#ss.example <- maSigPro(df_cd4_arm, design_CD4_arm, vars="each")

# next steps: get genes and make tidy df --> compute mean and SD as in previous scripts
# keep only genes that have their max at last position

#sigs <- get.siggenes(t, rsq = 0.6, vars = "all")
#genes <- sigs$summary

pipeline <- function(df, pdata, cell, inf){
  
  #make design matrix
  desmat <- make_des_matrix(pdata, cell, inf)
  # create design object for masigpro
  design <- make.design.matrix(desmat, degree = 2)
  
  # subset df to have same columns as design matrix rows (req. for masigpro)
  df <- df[rownames(desmat)]
  
  # run masigpro
  fit <- p.vector(df, design = design)
  t <- T.fit(fit, step.method = "backward", alfa = 0.05)
  
  # get signif genes and subset df
  sigs <- get.siggenes(t, rsq = 0.6, vars = "all")
  genes <- sigs$summary
  # subset
  df <- df[rownames(df) %in% genes,]
  

  return(df)
}

# run masigpro pipeline first
get_rtm_df <- function(df, pdata, cell, inf){

  # make tidy
  df$gene_name <- rownames(df)
  df <- pivot_longer(df, cols = -gene_name)
  df <- left_join(df, pdata, by = "name")
  df <- df %>% dplyr::select(-infection)
  
  # norm
  df <- norm_rtm(df)
  
  if(inf == "LCMV-Clone 13"){
    infstr <- "cl13"
  } else {
    infstr <- "arm"
  }
  write.csv(df, paste0("output/data_rtm_",cell, "_", infstr, ".csv"), row.names = F)
  
  return(df)
}


norm_rtm <- function(df){
  # compute avg per gene and time point
  df <- df %>% group_by(gene_name, time) %>% mutate(avg = mean(value)) %>% ungroup()
  df <- df %>% dplyr::select(-name)  
  
  # normalize to day 0
  df <- norm_d0(df)

  # change to response time dist get log and take abs values (downreg can also be resp. time prcess)
  df <- df %>% mutate(val_norm_rtm = abs(log2(val_norm)), avg_norm_rtm = abs(log2(avg_norm)))
  
  # get maximum val per gene over all times and divide by it
  df <- df %>% group_by(gene_name) %>%  mutate(max_gene = max(avg_norm_rtm)) %>% ungroup() 
  
  df <- df %>% mutate(avg_norm_rtm2 = avg_norm_rtm/max_gene, val_norm_rtm2 = val_norm_rtm/max_gene)
  
  # compute standard deviation and select columns needed for rtm
  df <- df %>% group_by(gene_name, time) %>% mutate(SD = sd(val_norm_rtm2)) %>% ungroup() %>%
    dplyr::select(gene_name, time, cell_type, SD, avg_norm_rtm2, val_norm_rtm2)
  
  return(df)
}

norm_d0 <- function(df){
  df_d0 <- df %>% filter(time == 0)
  df_d0 <- df_d0 %>% dplyr::select(-value, -time)
  df_d0 <- rename(df_d0, "avgd0" = "avg")
  df_d0 <- distinct(df_d0)
  df <- left_join(df, df_d0)
  df <- df %>% mutate(avg_norm = avg/avgd0, val_norm = value/avgd0)  
  return(df)
}

df_cd4_arm <- pipeline(df, pdata, cell = "CD4", inf = "LCMV-Arm")
df_cd4_cl13 <- pipeline(df, pdata, cell = "CD4", inf = "LCMV-Clone 13")
df_cd8_arm <- pipeline(df, pdata, cell = "CD8", inf = "LCMV-Arm")
df_cd8_cl13 <- pipeline(df, pdata, cell = "CD8", inf = "LCMV-Clone 13")


