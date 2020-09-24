# load crawford CD4 Arm data naive, d6,8,15,30
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
require(tidyr)
require(dplyr)
require(stringr)
require(ggplot2)

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

# filter feature data
fdata <- pData(featureData(gset))
fdata <- fdata %>% select(c("ID", "Gene title", "Gene symbol", "Gene ID"))

# remove rows with no entries
fdata2 <- fdata %>% filter(`Gene symbol` != "")

# choose gene subset
wei_th1 <- read_csv("gene_sets/references/wei_th1.CSV")
zhu_2012 <- read_csv("gene_sets/references/zhu_2012.CSV")
stubbington_th1 <- read_csv("gene_sets/references/stubbington_th1.CSV")
colnames(wei_th1)[1] <- "gene_symbol"
colnames(zhu_2012)[1] <- "gene_symbol"
colnames(stubbington_th1)[1] <- "gene_symbol"


gene_names <- wei_th1$gene_symbol
fdata_sub <- fdata2 %>% filter(`Gene symbol` %in% gene_names)

# subset of ex with corresponding probe IDs
ex_sub <- ex %>% filter(ID %in% fdata_sub$ID)
ex_sub <- pivot_longer(ex_sub, cols = -ID)

# add pheno data information
pdata <- pData(gset)
pdata <- rename(pdata, 
                c("name" = "geo_accession", 
                  "time" = "time (dpi):ch1", 
                  "infection" = "infection:ch1", 
                  "cell_type" = "cell type:ch1"))

pdata2 <- pdata %>% select(name, time, infection, cell_type)

# add feature info to ex data
ex_sub <- left_join(ex_sub, fdata_sub, by = "ID")
ex_sub <- left_join(ex_sub, pdata2, by = "name")

# need to have d0 as duplicate for arm and cl13
ex_sub$cell_type[(ex_sub$time == 0) & (ex_sub$cell_type == "Naive CD44Lo CD8+ T cells")] <- "H2-Db GP33-specific CD8+ T cells"
ex_sub$cell_type[(ex_sub$time == 0) & (ex_sub$cell_type == "Naive CD44Lo CD4+ T cell")] <- "H2-IAb GP66 specific CD4+ T cell"

# shorten cell type column
ex_sub <- ex_sub %>% mutate(cell_type = str_extract(cell_type, "CD."))

# add arm and cl13 infection for d0 cells even though they are naive
ex_sub$infection[ex_sub$time == 0] <- "LCMV-Arm"
ex_sub_d0 <- ex_sub %>% filter(time == 0)
ex_sub_d0$infection <- "LCMV-Clone 13"
ex_sub <- rbind(ex_sub, ex_sub_d0)

# compute avg expression per infection, cell type, time point and gene
ex_sub <- ex_sub %>% group_by(time, `Gene symbol`, infection, cell_type) %>% mutate(avg = mean(value)) %>% ungroup()

# normalize expression to day 0
ex_sub_d0 <- ex_sub %>% filter(time == 0)
ex_sub_d0 <- ex_sub_d0 %>% select(-value, -time, -name)
ex_sub_d0 <- rename(ex_sub_d0, "avgd0" = "avg")
# kick out duplicate rows that were induced because col removal
ex_sub_d0 <- distinct(ex_sub_d0)
ex_sub <- left_join(ex_sub, ex_sub_d0)

ex_sub <- ex_sub %>% mutate(avg_norm = avg / avgd0, val_norm = value / avgd0)

# change to response time dist subtract min
ex_sub <- ex_sub %>% mutate(val_norm_rtm = val_norm - 1.0, avg_norm_rtm = avg_norm - 1.0)

# divide by avg 
ex_sub <- ex_sub %>% group_by(`Gene symbol`, infection, cell_type) %>% mutate(val_norm_rtm2 = val_norm_rtm / max(avg_norm_rtm))
ex_sub <- ex_sub %>% group_by(`Gene symbol`, infection, cell_type) %>% mutate(avg_norm_rtm2 = avg_norm_rtm / max(avg_norm_rtm))

# add errors for resp time distr and norm vals
ex_sub <- ex_sub %>% group_by(cell_type, `Gene symbol`, infection, time) %>% mutate(err = sd(val_norm), err_rtm = sd(val_norm_rtm2))
#ex_sub <- ex_sub %>% filter(`Gene symbol` == "Ifng")

ex_sub$time <- as.numeric(ex_sub$time)

p1 <- ggplot(ex_sub, aes(time, val_norm))
p1 + 
  geom_point(aes(color = `Gene symbol`)) + 
  geom_line(aes(time, avg_norm, color = `Gene symbol`)) +
  facet_grid(infection~cell_type) +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave("figures/gene_kinetics.pdf")

p1 <- ggplot(ex_sub, aes(time, val_norm_rtm))
p1 + 
  geom_point(aes(color = `Gene symbol`)) + 
  geom_line(aes(time, avg_norm_rtm, color = `Gene symbol`)) +
  facet_grid(infection~cell_type) +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave("figures/gene_kinetics_rtm.pdf")

p1 <- ggplot(ex_sub, aes(time, val_norm_rtm2))
p1 + 
  geom_point(aes(color = `Gene symbol`)) + 
  geom_line(aes(time, avg_norm_rtm2, color = `Gene symbol`)) +
  facet_grid(infection~cell_type) +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave("figures/gene_kinetics_rtm.pdf")

# save a wide form df with times and averages for python processing 
python <- ex_sub %>% ungroup() %>% select(infection, cell_type, `Gene symbol`, time, avg_norm, err, avg_norm_rtm2, err_rtm)
python <- distinct(python)

#python2 <- split(python, list(python$infection, python$cell_type))
#python2 <- lapply(python2, pivot_wider, names_from = time, values_from = avg_norm)
#python2 <- bind_rows(python2)
write.csv(python, "output/avg_expression_norm.csv", row.names = F)
