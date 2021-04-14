# take processed data and run masigpro pipeline to identify kinetic genes
require(readr)
require(dplyr)
require(maSigPro)
#source("code/data_pipeline/utils.R")

prep_data <- function(df){
  # check that one column is named gene
  colnames(df)[grepl("X1", colnames(df))] <- "gene"
  colnames(df)[grepl("gene_name", colnames(df))] <- "gene"
  stopifnot("gene" %in% colnames(df))
  
  df <- as.data.frame(df)
  df <- df[!is.na(df$gene),]

  genes <- make.unique(df$gene)
  df <- df %>% dplyr::select(-gene)
  
  rownames(df) <- genes
  df
}


study <- "Crawford"

df <- read_csv(paste0("data/data_processed/", study, "_processed.csv"))
des_mat <- read.table(paste0("output/design_matrices/design_matrix_", study, ".txt"))

df2 <- prep_data(df)

stopifnot(length(rownames(des_mat) %in% colnames(df2))==ncol(df2))
design <- make.design.matrix(des_mat, degree = 2)

# run masigpro

# p vector compares intercept model with regression model. Tfit compares differences between fits
fit <- p.vector(df2, design = design)
genes <- rownames(fit$SELEC)
df_kinetic <- df2[rownames(df2) %in% genes,]
df_kinetic$gene <- rownames(df_kinetic)
write_csv(df_kinetic, paste0("data/data_kinetic/", study, "_kinetic.csv"))

#t <- T.fit(fit, step.method = "backward", alfa = 0.05)

# get signif genes and subset df
#sigs <- get.siggenes(t, vars = "all")
#genes <- sigs$summary

# subset

### using pvector output as kinetic genes now!!! date 03.02.2021

