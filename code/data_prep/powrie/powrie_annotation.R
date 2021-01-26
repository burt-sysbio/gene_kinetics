#   Differential expression analysis with limma
library(data.table)
library(readr)
require(dplyr)
require(tidyr)
require(stringr)
require(ggplot2)
## read and annotate data

illumina_mapping_file <- "annotation/MouseWG-6_V2_0_R3_11278593_A_mod.txt"
#anno_map <- read.table(illumina_mapping_file,header=TRUE)
anno_map <- fread(illumina_mapping_file, select = c("Array_Address_Id","Probe_Id","RefSeq_ID","Entrez_Gene_ID","Symbol"))
anno_map <- subset(anno_map,!is.na(Entrez_Gene_ID))
row.names(anno_map) <- anno_map$Array_Address_Id

source_file <- paste("data","data_raw","powrie", "powrie_time_course.matrix",  sep="/") 
data <- read.table(source_file,header=TRUE)


# change column names (except id column)
names <- colnames(data)
names <- names[1:length(names)-1]
names <- str_sub(names, 2, -1)
names <- str_replace(names, ".R", "_")
names <- paste("Powrie", "innate", names, sep = "_")
names <- c(names, "ids")

colnames(data) <- names

# check if data is log transformed
#illumina_ids <- paste("ILMN_",data$ids,sep="")
#mapped_probes <- mappedkeys(illuminaMousev2GENENAME)
#xx <- as.list(x[mapped_probes])
#gene_symbols <- matrix(unlist(xx))

probes <- intersect(data$ids, anno_map$Array_Address_Id)
probe_idx <- match(probes,data$ids)
map_idx <- match(probes,anno_map$Array_Address_Id)

t_data <- data[probe_idx,1:ncol(data)-1]
t_data_annot <- cbind(data[probe_idx,],anno_map[map_idx,])

row.names(t_data) <- t_data_annot$ids
row.names(t_data_annot) <- t_data_annot$ids

t_data_annot <- t_data_annot %>% dplyr::select(-ids, -Entrez_Gene_ID, -RefSeq_ID, -Array_Address_Id)
cols <- colnames(t_data)
df_tidy <- pivot_longer(t_data_annot, cols <- cols, values_to = "expr", names_to = "sample")

# get median out of genes that match to multiple probe ids
df_tidy <- df_tidy %>% group_by(Symbol, sample) %>% 
  mutate(expr = median(expr)) %>% ungroup() 

df_tidy <- df_tidy %>% dplyr::select(-Probe_Id)
df_tidy <- distinct(df_tidy)

df_tidy <- df_tidy %>% mutate(sample = str_sub(sample, 2, -4)) %>% rename(time = sample, gene_name = Symbol)

g <- ggplot(data = df_tidy, aes(x = sample, y = expr))
g + geom_boxplot()

# make to long form again
df_proc <- df_tidy %>% pivot_wider(id_cols = Symbol, names_from = sample, values_from = expr)
colnames(df_proc)[1] <- "gene"

write.csv(df_proc, "data/data_processed/powrie_processed.csv", row.names = F)
# check if some probes map to multiple genes
