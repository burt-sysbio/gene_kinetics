require(ggplot2)
require(tidyr)
require(dplyr)

df_ann <- read.csv("output/crawford_CD4_anndata.csv")

gene_names <- c("Ifng", "Tbx21", "Gata3", "Il2", "Il21", "Il10", "Tnf")
test <- df_ann %>% filter(Gene.symbol %in% gene_names)
test <- pivot_longer(test, -c("Gene.title", "Gene.symbol", "ID", "Gene.ID"))
test <- separate(test, col = name, into = c("celltype", "day", "infection", "rep"), sep = "_")
test$day <- as.numeric(test$day)

avg <- test %>% group_by(day, Gene.symbol) %>% summarize(avg = mean(value))
test <- left_join(test, avg, by= c("day", "Gene.symbol"))

d0 <- test %>% filter(day == 0)
d0 <- d0 %>% select(-day, -value)
d0 <- rename(d0, c("avg_d0" = "avg"))

test <- left_join(test, d0)
test <- test %>% mutate(avg_norm = avg / avg_d0, val_norm = value / avg_d0)

p1 <- ggplot(data = test, aes(day, val_norm))
p1 + 
  geom_point(aes(color = Gene.symbol)) + 
  geom_line(aes(day, avg_norm, color = Gene.symbol))
