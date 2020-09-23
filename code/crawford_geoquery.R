# load crawford CD4 Arm data naive, d6,8,15,30
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
require(tidyr)
require(dplyr)
# load series and platform data from GEO
# load series and platform data from GEO

gset <- getGEO("GSE41870", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("0000111222233334444XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }



fdata <- pData(featureData(gset))
fdata <- fdata %>% select(c("ID", "Gene.title", "Gene.symbol", "Gene.ID"))
adata <- exprs(gset)

# remove rows with no entries
fdata2 <- fdata %>% filter(Gene.symbol != "")

# merge feature data and expression data
adata <- as.data.frame(adata)

# assign colnames
n_naive <- 4
n_d6 <- 3
n_d8 <- 4 
n_d15 <- 4
n_d30 <- 4

g1 <- paste(rep("CD4_0", n_naive), "Arm", 1:n_naive, sep = "_")
g2 <- paste(rep("CD4_6", n_d6), "Arm", 1:n_d6, sep = "_")
g3 <- paste(rep("CD4_8", n_d8), "Arm", 1:n_d8, sep = "_")
g4 <- paste(rep("CD4_15", n_d15), "Arm", 1:n_d15, sep = "_")
g5 <- paste(rep("CD4_30", n_d30), "Arm", 1:n_d30, sep = "_")

colnames <- c(g1,g2,g3,g4,g5)
colnames(adata) <- colnames

# add ID column for merging
adata$ID <- rownames(adata)
adata$ID <- as.numeric(adata$ID)

df_ann <- inner_join(adata, fdata2, by = c("ID"))

write.csv(df_ann, "output/crawford_CD4_anndata.csv", row.names = F)


