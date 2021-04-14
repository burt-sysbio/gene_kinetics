require(dplyr)
require(tidyr)
require(readr)
require(stringr)

# prepare design matrix from processed data frame
get_des_mat <- function(df){
  print("double check shared starting points, not adjusted for crawford")

  
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
  #df <- df %>% arrange(Time)
  # add replicate information, always check if title2 (info wo replicate)
  #is the same as in row before
  j <-1
  myseq <- 2:nrow(df)
  df <- df %>% arrange(Time)
  print(df)
  for (i in myseq){
    if (df[i, "title2"]!= df[i-1, "title2"]){
      j <- j+1
    }
    df[i, "Replicate"] <- j
  }
  
  # reformat for design matrix style
  rownames <- df$title
  df <- df %>% dplyr::select(-c(title, title2)) %>% dplyr::select(Time, Replicate, everything())
  rownames(df) <- rownames
  
  return(df)
}


study <- "Crawford"
savedir <- c("output/design_matrices/")
df <- read_csv(paste0("data/data_processed/", study, "_processed.csv"))

df <- df %>% select(-gene)
des_mat <- get_des_mat(df)
rownames <- rownames(des_mat)

des_mat_idx <- des_mat[, 3:ncol(des_mat)] != 0
des_mat[,3:ncol(des_mat)][des_mat_idx] <- 1

if(study == "Crawford"){
  
  grep_col <- grepl("None", colnames(des_mat))
  des_mat = des_mat[,!grep_col]
  des_mat[1:4,3:4] = 1
  des_mat[5:8,5:6] = 1
}

if(study == "Proserpio"){
  des_mat$Th2 <- 1
}
  

rownames(des_mat) <- rownames

write.table(des_mat, paste0(savedir, "design_matrix_", study, ".txt"))
print(des_mat)
