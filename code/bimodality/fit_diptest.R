require(msigdbr)
require(dplyr)
require(tidyr)
require(readr)
require(stringr)
require(diptest)
require(RColorBrewer)


# load own files
fit_summary <- read.csv("output/gamma_fits/feb2021/fit_summary_all.csv")
fit_summary <- fit_summary %>% select(c("gene", "study", "best_fit"))


# load own files
filenames <- list.files("data/data_rtm/", pattern="*.csv", full.names=TRUE)
idx <- grepl("_rtm_", filenames)
names <- filenames[idx]
files <- lapply(names, read.csv)

names2 <- str_sub(names, 15, -5)

# add study name to data frame, then concatenate rows

myfun <- function(data, name){
  data$study <- name
  data <- data %>% select("study", "gene", "time", "avg_norm_rtm2") %>% unique()
  data <- data %>% pivot_wider(names_from = time, values_from = avg_norm_rtm2)
  
  crit <- ((data[,3] < data[,4]) & (data[,4] > data[,5])) | ((data[,3] > data[,4]) & (data[,4] < data[,5]))
  
  # get diference between first and snd and snd and thrd timepoitn
  data$dip <- crit[,]
  
  delta1 <- abs(data[,4] - data[,3])
  delta2 <- abs(data[,5] - data[,4])
  df <- bind_cols(delta1, delta2)
  df$deltamin <- apply(df, 1, min)
  # assign dip value for each gene where a dip was found
  data$dip_idx <- df$deltamin
  data$dip_idx[data$dip == FALSE] <- NaN
  #data <- data %>% select("gene", "study", "dip")
  data
}


data_list <- mapply(myfun, files, names2)

data_list2 <- lapply(data_list, function(x){x %>% select("gene", "study", "dip", "dip_idx")})

data_list2 <- bind_rows(data_list2)

test <- left_join(fit_summary, data_list2, by = c("gene", "study"))

# how many dips for each fit assignment category across studies?
test2 <- test %>% group_by(best_fit) %>% summarize(ndips = sum(dip))

test3 <- test %>% filter(dip == T)

require(ggplot2)

pal1 <- brewer.pal(8, "Dark2")
pal2 <- brewer.pal(8, "Paired")

p <- ggplot(test3, aes(best_fit, dip_idx, fill = best_fit))
p1 <- p + geom_boxplot(outlier.shape=NA) + 
  xlab("") +
  ylab("dip index") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c(pal1[1], pal2[2], pal1[2], "grey"))

print(p1)

ggsave("figures/bimodality.png", width = 4, height = 3)
#for(i in seq_along(files)){
#  data <- files[[i]]
#  data$study <- names2[i]

  # only keep averages
#  data <- data %>% select("study", "gene", "time", "avg_norm_rtm2") %>% unique()
#  data <- data %>% pivot_wider(names_from = time, values_from = avg_norm_rtm2)

#  crit <- ((data[,3] < data[,4]) & (data[,4] > data[,5])) | ((data[,3] > data[,4]) & (data[,4] < data[,5]))
#  data$dip <- crit[,]
#}
