# postprocess kinetic genes to generate data for rtm fitting
require(dplyr)
require(readr)
source("code/data_pipeline/utils.R")


### read data
study <- "Proserpio"

infection <- "parasite"

df <- read.csv(paste0("data/data_kinetic/", study, "_kinetic.csv"))

# split into groups
df_tidy <- make_tidy(df)

if(study == "Crawford"){
  d0_cd4_cl13 <- df_tidy %>% filter(cell_type == "CD4.None")
  d0_cd8_cl13 <- df_tidy %>% filter(cell_type == "CD8.None")
  
  d0_cd4_cl13$cell_type <- "CD4.cl13"
  d0_cd8_cl13$cell_type <- "CD8.cl13"
  
  # rename existing none column to arm and then add the duplicated cl13 columns
  df_tidy$cell_type <- str_replace(df_tidy$cell_type, "CD4.None", "CD4.arm")
  df_tidy$cell_type <- str_replace(df_tidy$cell_type, "CD8.None", "CD8.arm")
  
  df_tidy <- bind_rows(d0_cd4_cl13, d0_cd8_cl13, df_tidy)
}

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
  file <- paste0(savedir, group, "_", infection, ".csv")
  write_csv(df_out, file)
}

