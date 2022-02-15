require(dplyr)
require(stringr)
require(DESeq2)
require(biomaRt)
require(ggplot2)
require(tidyr)

print("DESeq2 and biomart not installed at ITB computers")

# load files
data_dir <- "data/data_raw/proserpio/E-MTAB-3543.processed.1/"
files <- list.files(data_dir)
files2 <- paste0(data_dir, files)

myfiles = lapply(files2, read.delim, header = FALSE, col.names = c("gene", "count"))

# only keep count data in myfiles, store gene info extra
gene_names <- myfiles[[1]]$gene
myfiles = lapply(myfiles, function(x) { x["gene"] <- NULL; x })
df_raw <- bind_cols(myfiles)

# format column names for processed output
times <- str_sub(files, start = 5, end = -12)
reps <- rep(c("A", "B"), 6)
cell <- "Th2"

colnames <- paste("Proserpio", cell, times, reps, sep = "_")
colnames(df_raw) <- colnames

# filter data set with condition > median epxression
countmat <- as.matrix(df_raw)
rownames(countmat) <- gene_names


# keep all genes where expr is greater 1 in one condition
keep <- rowSums(countmat) > 1
countmat <- countmat[keep,]
print(dim(countmat))

genes_filtered <- rownames(countmat)


counts_tidy <- pivot_longer(as.data.frame(countmat), cols = everything())
counts_tidy <- counts_tidy %>% mutate(val_norm = log2(value+1))
g <- ggplot(data = counts_tidy, aes(x = name, y = val_norm)) +
  geom_boxplot() + labs(x = "", y = "expression")
print(g)

# variance transform for rnaseq data

trafo <- "log2"

if(trafo == "log2"){
  vst_data = log2(countmat + 1)
  sname <- "log2trafo"
} else {
  vst_data <- rlog(countmat)  
  sname <- "rlogtrafo"
}

# annotation
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
res <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), mart = mouse)
colnames(res) <- c("ensembl_id", "gene")

df_proc <- as.data.frame(vst_data)
df_proc$ensembl_id <- genes_filtered

df_proc <- left_join(df_proc, res, by = "ensembl_id")
df_proc <- df_proc %>% dplyr::select(-ensembl_id)
df_proc <- na.omit(df_proc)

df_proc <- df_proc %>% dplyr::select(gene, everything())

savename <- paste0("data/data_processed/Proserpio_processed_", sname, ".csv")
write.csv(df_proc, file = savename, row.names = F)