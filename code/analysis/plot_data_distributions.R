# plot distribution of normalized data sets before response-time processing

require(readr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(stringr)

mypath <- "data/data_processed/"
readnames <- list.files(path = mypath, pattern = "processed", full.names = T)

readfun <- function(readname){
  df <- read_csv(readname)
  df <- df %>% pivot_longer(-gene)
  df$study <- str_sub(readname, 21, -15)
  return(df)
}

df_list <- lapply(readnames, readfun)

df <- bind_rows(df_list)

g <- ggplot(data = df, aes(x = study, y = value)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(ylim = c(0,18)) + labs(y = "expr. norm.")
print(g)
ggsave("figures/data_distributions/data_processed_norm_log2.png", bg = "white", width =  3, height = 3)

