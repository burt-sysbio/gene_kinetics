# load crawford CD4 Arm data naive, d6,8,15,30
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
require(tidyr)
require(dplyr)
require(stringr)
require(ggplot2)
require(readr)
require(biomaRt)

# load series and platform data from GEO
gset <- getGEO("GSE41870", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# remove last 8 columns / samples because not kinetic data
gset <- gset[,1:(ncol(exprs(gset))-8)]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# merge feature data and expression data
ex <- as.data.frame(ex)
ex$ID <- as.numeric(rownames(ex))


# get annotation
ensembl <- useMart("ensembl")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
gene_ids <- ex$ID

foo <- getBM(attributes=c('affy_mogene_1_0_st_v1','external_gene_name', 'ensembl_gene_id', "entrezgene_id"),
             filters = 'affy_mogene_1_0_st_v1',
             values = gene_ids,
             mart = ensembl)

# kick out most non protein coding genes (should have no entrezgene id, at least for pseudogenes)
foo <- na.omit(foo)

# kick out genes that have duplicate gene annotations for same probe id
dups <- foo$'external_gene_name'[duplicated(foo[1])]

foo <- foo %>% filter(!(external_gene_name %in% dups))
colnames(foo) <- c("ID", "gene_name", "ensembl", "entrez")
foo <- foo %>% dplyr::select(ID, gene_name)

# add pheno data information
pdata <- pData(gset)
pdata <- rename(pdata, 
                c("name" = "geo_accession", 
                  "time" = "time (dpi):ch1", 
                  "infection" = "infection:ch1", 
                  "cell_type" = "cell type:ch1"))

pdata2 <- pdata %>% dplyr::select(name, time, infection, cell_type)


# subset of ex with corresponding probe IDs
ex_sub <- ex %>% dplyr::filter(ID %in% foo$ID)
ex_sub <- pivot_longer(ex_sub, cols = -ID)

# add feature info to ex data
ex_sub <- left_join(ex_sub, pdata2, by = "name")
ex_sub <- left_join(ex_sub, foo, by = "ID")

# take median for genes that match multiple probes for each sample (each timepoint/celltype/infection/replicate)
ex_sub <- ex_sub %>% group_by(name, gene_name) %>% mutate(value2 = median(value)) %>% ungroup()
# kick out unneccessary columns
ex_sub2 <- ex_sub[c("name", "gene_name", "value2")]
# keep unique rows
ex_sub2 <- distinct(ex_sub2)

# make wide again and keep only genes above median in at least one condition
ex_sub3 <- pivot_wider(ex_sub2, names_from = name, values_from = value2)
names <- ex_sub3$gene_name
ex_sub3 <- dplyr::select(ex_sub3, -gene_name)
rownames(ex_sub3) <- names

ex_sub4 <- apply(ex_sub3, 2, as.numeric)
rownames(ex_sub4) <- names

rowmax <- apply(ex_sub4, 1, max)
ex_sub4 <- ex_sub4[rowmax >= median(ex_sub4),]
ex_sub4 <- as.data.frame(ex_sub4)

# rename cells not a mistake that its cell and cells below
pdata2$cell_type[pdata2$cell_type == "Naive CD44Lo CD8+ T cells"] <- "H2-Db GP33-specific CD8+ T cells"
pdata2$cell_type[pdata2$cell_type == "Naive CD44Lo CD4+ T cell"] <- "H2-IAb GP66 specific CD4+ T cell"
pdata2 <- pdata2 %>% mutate(cell_type = str_extract(cell_type, "CD."))
pdata2$infection[pdata2$infection == "LCMV-Arm"] <- "arm"
pdata2$infection[pdata2$infection == "LCMV-Clone 13"] <- "cl13"

# group identifier separated by -, separate time and study and replicate by _
groups <- paste(pdata2$cell_type, pdata2$infection, sep = "-")
title <- paste("craw", groups, pdata2$time, sep = "_")

# add the replicate information
title <- make.unique(title)
idx <- !grepl("\\.", title)
title[idx] <- paste(title[idx], "1", sep = "_")
title <- sub("\\.1", "_2", title)
title <- sub("\\.2", "_3", title)
title <- sub("\\.3", "_4", title)
colnames(ex_sub4) <- title

# set gene name as column for csv
ex_sub4$gene <- rownames(ex_sub4)
ex_sub4 <- ex_sub4 %>% dplyr::select(gene, everything())
write.csv(ex_sub4, "data/data_processed/crawford_processed.csv", row.names = F)

