# get Yosef Nir Th17 diff. timecourse and make annotation
# filter expression values
library(GEOquery)
require(dplyr)
require(tidyr)
require(readr)

# load series and platform data from GEO
gset <- getGEO("GSE43970", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8321", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# get expression and make log transform
ex = exprs(gset)
boxplot(ex)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  ex <- log2(ex) 
}

# check value distribution
# box-and-whisker plot
#boxplot(ex)

# get phenoData
pdata <- pData(gset)
fdata <- fData(gset)

# exclude the following conditions from analysis
excl <- c("KO|WT|IL23|repeat")
title <- pdata$title
idx <- !grepl(excl, title)
p <- pdata[idx,]

# keep only those columns of expression set that are are filtered pdata
cols <- colnames(ex)
col_idx <- cols %in% rownames(p)
ex = ex[,col_idx]

# make column identifiers
time <- p$`time (hr):ch1`
group <- p$`treatment:ch1`
group[group == "Tgfb+Il6"] <- "Th17"
col_id <- paste("nir", group, time, "1", sep ="_")
colnames(ex) <- col_id


# format gene symbols
genes <- fdata %>% separate(`Gene Symbol`, into = c("gene", "sym2"), sep = "///")
# rename genes with different isoforms to get unique identifiers
genes <- genes %>% mutate(gene = make.unique(gene, sep = "_"))

# add gene names to expression set
g <- genes$gene
rownames(ex) <- g

# keep only genes where expr is above median in at least one condition
m <- median(ex)
crit <- apply(ex, 1, max)
# get rowmax
ex <- ex[crit > m,]
ex <- as.data.frame(ex)
ex$gene <- rownames(ex)
ex <- ex %>% dplyr::select(gene, everything())
write_csv(ex, "data/data_processed/nir_processed.csv")