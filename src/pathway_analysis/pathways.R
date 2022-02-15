require(ontologyIndex)
require(clusterProfiler)
require(msigdbr)
require(readr)
require(dplyr)


# needs to be dataframe with gene as column
mydata <- read.csv("output/gamma_fits/feb2021/category_assignment_winners.csv")

deg_gamma <- mydata %>% filter(winner == "gamma") %>% select(gene) %>% unique() 
deg_expo <- mydata %>% filter(winner == "expo") %>% select(gene) %>% unique()
universe <- as.character(unique(mydata$gene))

files <- c(deg_gamma, deg_expo)
filenames_save <- c("gamma.csv", "expo.csv")

# load ontologies
cat1 <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS")
cat2 <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
cat3 <- msigdbr(species = "Mus musculus", category = "H")
cat4 <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
cat5 <- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:TFT_Legacy")

categories <- list(cat1, cat2, cat3, cat4, cat5)
cat_names <- c("WIKIPATHWAYS", "REACTOME", "HALLMARK", "GOBP", "TFT")

all_sets <- bind_rows(categories)
all_sets_nogo <- all_sets[all_sets$gs_subcat != "GO:BP",]
# pathway analysis for each file and each data base

sdir <- "data/ORA_output/"

genesets1 <- as.data.frame(all_sets) %>% select(gs_name, gene_symbol)
genesets2 <- as.data.frame(all_sets_nogo) %>% select(gs_name, gene_symbol)

for(i in seq_along(files)){
  savename <- filenames_save[i]
  degs <- files[i]
  
  # run ora for all categories combined
  print("")
  print("not using universe at the moment")
  print("")
  
  test <- enricher(degs$gene, TERM2GENE = genesets1, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1)
  res <- test@result
  out <- paste0(sdir, "ORA_combined_", savename)
  write.table(res, out)
  
  test <- enricher(degs$gene, TERM2GENE = genesets2, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1)
  res <- test@result
  out <- paste0(sdir, "ORA_nogo_", savename)
  write.table(res, out)
  # run ORA for each category separately
  for(j in seq_along(categories)){
    genesets <- categories[[j]]

    catname <- cat_names[[j]]
    genesets <- as.data.frame(genesets) %>% select(gs_name, gene_symbol)
    
    test <- enricher(degs$gene, TERM2GENE = genesets, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1)
    res <- test@result
    
    out <- paste0(sdir, "ORA_", catname, "_", savename)
    write.table(res, out)
  }
}


