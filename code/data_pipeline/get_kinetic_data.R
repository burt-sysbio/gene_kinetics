# take processed data and run masigpro pipeline to identify kinetic genes
require(readr)
require(dplyr)
require(maSigPro)
source("code/data_pipeline/utils.R")

df <- read_csv("data/data_processed/nir_processed.csv")
df <- prep_data(df)
des_mat <- get_des_mat(df)

design <- make.design.matrix(des_mat, degree = 2)

# run masigpro
fit <- p.vector(df, design = design)
t <- T.fit(fit, step.method = "backward", alfa = 0.05)

# get signif genes and subset df
sigs <- get.siggenes(t, vars = "all")
genes <- sigs$summary

# subset
df_kinetic <- df[rownames(df) %in% genes,]

df_kinetic$gene <- rownames(df_kinetic)
write_csv(df_kinetic, "data/data_kinetic/nir_kinetic.csv")