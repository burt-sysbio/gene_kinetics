require(dplyr)
require(tidyr)
require(stringr)

# prepare processed data frame (data_processed) to run masigpro
prep_data <- function(df){
  df <- as.data.frame(df)
  df <- df[!is.na(df$gene),]
  genes <- make.unique(df$gene)
  df <- df %>% dplyr::select(-gene)
  rownames(df) <- genes
  df
}

# prepare design matrix from processed data frame
get_des_mat <- function(df){
  print("double check shared starting points, not adjusted for crawford")
  names <- colnames(df)
  
  # split column titles into identifiers and convert to wide format
  # fill missing values with 0, required for masigpro format
  df <- tibble("title" = colnames(df))
  df <- df %>% separate("title", 
                        c("study", "group", "Time", "Replicate"), 
                        sep = "_",remove = F)
  
  df$Replicate <- as.numeric(df$Replicate)
  df$Time <- as.numeric(df$Time)
  df <- df %>% dplyr::select(-study)
  df <- df %>% pivot_wider(id_cols = c("Time", "title"), 
                           names_from = group, 
                           values_from = Replicate,
                           values_fill = 0)
  
  df <- df %>% mutate(title2 = str_sub(title, 1, -3))
  df$Replicate <- 1
  
  # add replicate information, always check if title2 (info wo replicate)
  #is the same as in row before
  j <-1
  myseq <- 2:nrow(df)
  for (i in myseq){
    if (df[i, "title2"]!= df[i-1, "title2"]){
      j <- j+1
    }
    df[i, "Replicate"] <- j
  }
  
  # reformat for design matrix style
  df <- df %>% dplyr::select(-c(title, title2)) %>% 
    dplyr::select(Time, Replicate, everything())
  rownames(df) <- names
  
  return(df)
}


# take output from masigpro kinetic genes and make tidy
make_tidy <- function(df){
  df <- df %>% pivot_longer(-gene)
  df <- df %>% separate(name, into = c(NA, "cell_type", "time", NA), sep ="_")
  df$time <- as.numeric(df$time)
  df
}


norm_rtm <- function(df){
  # provide tidy data frame from make tidy fun output
  
  # compute avg per gene and time point
  df <- df %>% group_by(gene, time) %>% mutate(avg = mean(value)) %>% ungroup()
  
  # normalize to day 0
  df <- norm_d0(df)
  
  # change to response time dist get log and take abs values (downreg can also be resp. time prcess)
  df <- df %>% mutate(val_norm_rtm = abs(log2(val_norm)), 
                      avg_norm_rtm = abs(log2(avg_norm)))
  
  # get maximum val per gene over all times and divide by it
  df <- df %>% group_by(gene) %>%  
    mutate(max_gene = max(avg_norm_rtm)) %>% 
    ungroup() 
  
  df <- df %>% mutate(avg_norm_rtm2 = avg_norm_rtm/max_gene, 
                      val_norm_rtm2 = val_norm_rtm/max_gene)
  
  # compute standard deviation and select columns needed for rtm
  df <- df %>% group_by(gene, time) %>% 
    mutate(SD = sd(val_norm_rtm2)) %>% ungroup() %>%
    dplyr::select(gene, time, cell_type, SD, avg_norm_rtm2, val_norm_rtm2)
  
  return(df)
}

# normalize data to day0
norm_d0 <- function(df){
  df_d0 <- df %>% filter(time == min(time))
  df_d0 <- df_d0 %>% dplyr::select(-value, -time)
  df_d0 <- dplyr::rename(df_d0, "avgd0" = "avg")
  df_d0 <- distinct(df_d0)
  df <- left_join(df, df_d0)
  df <- df %>% mutate(avg_norm = avg/avgd0, val_norm = value/avgd0)  
  return(df)
}
