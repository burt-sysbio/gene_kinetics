require(ontologyIndex)
require(clusterProfiler)
require(msigdbr)
require(readr)
require(dplyr)
require(tidyr)
require(stringr)

myfiles <- list.files("output/gamma_fits/feb2021/", pattern = "fit_summary")
myfiles <- myfiles[grepl("Peine|Nir|Proserpio", myfiles)]
myfiles <- paste0("output/gamma_fits/feb2021/", myfiles)

df_list = c()
for (f in myfiles){
  df <- read.csv(f)
  df$study <- str_sub(f, 39, -5)
  
  if(str_detect(f, "Nir")){
    df <- df %>% separate(gene, into = c("gene", "iso1", "iso2"), sep = "_",remove = TRUE)
  } else {
    df <- df %>% separate(gene, into = c("gene", "iso1", "iso2"), sep = "\\.",remove = TRUE)
  }

  df$gene <- str_trim(df$gene)
  df <- df %>% dplyr::filter(gene != "")
  df <- df %>% dplyr::select(gene, best_fit, study) %>% unique()
  
  df_list <- c(df_list, list(df))
}


# load ontologies
cat1 <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS")
cat2 <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
cat3 <- msigdbr(species = "Mus musculus", category = "H")
cat4 <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
cat5 <- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:TFT_Legacy")

categories <- list(cat1, cat2, cat3, cat4, cat5)
cat_names <- c("WIKIPATHWAYS", "REACTOME", "HALLMARK", "GOBP", "TFT")

all_sets <- dplyr::bind_rows(categories)
all_sets_nogo <- all_sets[all_sets$gs_subcat != "GO:BP",]

# do analyisis with and wo GOBP
genesets1 <- as.data.frame(all_sets) %>% select(gs_name, gene_symbol)
genesets2 <- as.data.frame(all_sets_nogo) %>% select(gs_name, gene_symbol)


# for each study, check enrichment for gamma expo and other category
sdir <- "data/ORA_output/"

for (mydata in df_list){
  
  
  deg_gamma <- mydata %>% filter(best_fit == "gamma") %>% select(gene) %>% unique() 
  deg_expo <- mydata %>% filter(best_fit == "expo") %>% select(gene) %>% unique()
  deg_other <- mydata %>% filter(best_fit == "other") %>% select(gene) %>% unique()
  universe <- as.character(unique(mydata$gene))
  
  files <- c(deg_gamma, deg_expo, deg_other)
  filenames_save <- c("gamma.csv", "expo.csv", "other.csv")
  
  mystudy <- mydata$study[1]
  
  savedir <- paste0(sdir, mystudy, "_")
  for(i in seq_along(files)){
    savename <- filenames_save[i]
    degs <- files[i]
    
    # run ora for all categories combined
    test <- enricher(degs$gene, TERM2GENE = genesets1, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1, universe = universe)
    res <- test@result
    out <- paste0(savedir, "ORA_combined_", savename)
    write.table(res, out)
    
    test <- enricher(degs$gene, TERM2GENE = genesets2, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1, universe = universe)
    res <- test@result
    out <- paste0(savedir, "ORA_nogo_", savename)
    write.table(res, out)
    
    # run ORA for each category separately
    for(j in seq_along(categories)){
      genesets <- categories[[j]]
      
      catname <- cat_names[[j]]
      genesets <- as.data.frame(genesets) %>% select(gs_name, gene_symbol)
      
      test <- enricher(degs$gene, TERM2GENE = genesets, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1, universe = universe)
      res <- test@result
      
      
      out <- paste0(savedir, "ORA_", catname, "_", savename)
      write.table(res, out)
    }
  }
}

