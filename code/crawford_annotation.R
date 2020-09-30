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

# take median for genes that match multiple probes
ex_sub <- ex_sub %>% group_by(cell_type, time, gene_name, name, infection) %>% mutate(value2 = median(value)) %>% ungroup()
# kick out unneccessary columns
ex_sub <- ex_sub %>% dplyr::select(-ID, -value, -name)
# keep unique rows
ex_sub <- distinct(ex_sub)


# need to have d0 as duplicate for arm and cl13
ex_sub$cell_type[(ex_sub$time == 0) & (ex_sub$cell_type == "Naive CD44Lo CD8+ T cells")] <- "H2-Db GP33-specific CD8+ T cells"
ex_sub$cell_type[(ex_sub$time == 0) & (ex_sub$cell_type == "Naive CD44Lo CD4+ T cell")] <- "H2-IAb GP66 specific CD4+ T cell"

# shorten cell type column
ex_sub <- ex_sub %>% mutate(cell_type = str_extract(cell_type, "CD."))
# add arm and cl13 infection for d0 cells even though they are naive
ex_sub$infection[ex_sub$time == 0] <- "LCMV-Arm"
ex_sub_d0 <- ex_sub %>% filter(time == 0)
ex_sub_d0$infection <- "LCMV-Clone 13"
ex_sub2 <- rbind(ex_sub, ex_sub_d0)

write.csv(ex_sub2, "output/crawford_annotation.csv", row.names = F)
