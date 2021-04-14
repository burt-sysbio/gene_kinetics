require(dplyr)
require(tidyr)
require(stringr)



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
  
  # change to response time dist and get log
  df <- df %>% mutate(val_norm_rtm = log2(val_norm), 
                      avg_norm_rtm = log2(avg_norm))
  
  # get maximum or minimum (new) to include downreg val per gene over all times and divide by it
  df <- df %>% group_by(gene) %>%  
    mutate(max_gene = ifelse(max(avg_norm_rtm) > abs(min(avg_norm_rtm)), max(avg_norm_rtm), min(avg_norm_rtm))) %>% 
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
