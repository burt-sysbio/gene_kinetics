
# load data
require(readr)
require(dplyr)
require(tidyr)
# run masigpro pipeline first
get_rtm_df <- function(df, pdata, cell, inf){
  
  # make tidy
  df$gene_name <- rownames(df)
  df <- pivot_longer(df, cols = -gene_name)
  df <- left_join(df, pdata, by = "name")
  df <- df %>% dplyr::select(-infection)
  
  # norm
  df <- norm_rtm(df)

  write.csv(df, paste0("output/peine_data_rtm_",cell,".csv"), row.names = F)
  
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

# read data , note that rownames are in col X
df_peine <- read.csv("data/peine_data_median_threshold_annot.csv", row.names = "X")

timepoints <- rep(c(0,3,6,12,24,35,48,73,96,120), 4)
celltype <- rep(c("Th0", "Th1", "Th2", "Thm"), each = 10)

pdata <- data.frame(name = colnames(df_peine), time = timepoints, infection = "in vitro", cell_type = celltype)

df_th0 <- df_peine[,1:10]
df_th1 <- df_peine[, 11:20]
df_th2 <- df_peine[, 21:30]
df_thm <- df_peine[, 31:40]

# get kinetic genes
th1_kinetic <- na.omit(read.csv("gene_sets/peine_kinetic/out_table_Th1.csv"))
th2_kinetic <- na.omit(read.csv("gene_sets/peine_kinetic/out_table_Th2.csv"))
th0_kinetic <- na.omit(read.csv("gene_sets/peine_kinetic/out_table_Th0.csv"))
thm_kinetic <- na.omit(read.csv("gene_sets/peine_kinetic/out_table_ThMix.csv"))

fun1 <- function(df, pdata, celltype){
  file <- paste0("gene_sets/peine_kinetic/out_table_", celltype, ".csv")
  df_kinetic <- na.omit(read.csv(file))
  df <- df[rownames(df) %in% df_kinetic[["Symbol"]],]
  df_rtm <- get_rtm_df(df, pdata, celltype, "in vitro")
}

df_th1_rtm <- fun1(df_th1, pdata, "Th1")
df_th2_rtm <- fun1(df_th2, pdata, "Th2")
df_th0_rtm <- fun1(df_th0, pdata, "Th0")
df_thm_rtm <- fun1(df_thm, pdata, "ThMix")