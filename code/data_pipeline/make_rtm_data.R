# postprocess kinetic genes to generate data for rtm fitting
require(dplyr)
require(readr)
source("code/data_pipeline/utils.R")


### read data
study <- "nir"
df <- read.csv(paste0("data/data_kinetic/", study, "_kinetic.csv"))

# split into groups
df_tidy <- make_tidy(df)
df_list <- split(df_tidy, df_tidy$cell_type)

# transform for RTM
df_list_rtm <- lapply(df_list, norm_rtm)

# save output
groups <- names(df_list_rtm)
groups <- paste(study, "rtm", groups, sep = "_")
savedir <- "data/data_rtm/"

for(i in seq_along(df_list_rtm)){
  group <- groups[i]
  df_out <- df_list_rtm[[i]]
  file <- paste0(savedir, group, ".csv")
  write_csv(df_out, file)
}

