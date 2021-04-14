require(msigdbr)
require(dplyr)
require(tidyr)
require(readr)


check_pathway <- function(pathway_genes, del_genes, nodel_genes){
  # double check that I dont have duplicates with changed names _1, _2, etc
  del_genes = del_genes %>% separate(gene, into = c("gene", "isoform"), sep = "_") %>% select(gene) %>% unique()
  nodel_genes = nodel_genes %>% separate(gene, into = c("gene", "isoform"), sep = "_") %>% select(gene) %>% unique()
  
  # get the genes that are in the pathway and those that are not in the pathway for both delay and nodelay genes
  n_del_in <- sum(del_genes$gene %in% pathway_genes)
  n_nodel_in <- sum(nodel_genes$gene %in% pathway_genes)
  n_del_out <- nrow(del_genes) - n_del_in
  n_nodel_out <- nrow(nodel_genes) - n_nodel_in 
  
  # run fishers exact test to test if odds ratio is greater than expected (enriched) in delay genes
  mat <- matrix(c(n_del_in, n_nodel_in, n_del_out, n_nodel_out), 2, 2)
  pval <- fisher.test(mat, alternative="greater")$p.value
  #if(pval<0.05){
  #  print(mat)
  #}
  pval
}

# load some target gene sets from data base
#s1 = msigdbr(species = "Mus musculus", category = "H")
#s2 = msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:TFT_Legacy")
#s3 = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
database <- "C7"
s4 = msigdbr(species = "Mus musculus", category = database)


#
patterns <- c("_T_CELL", "TCELL", "CD3", "CD28", "CD4", "CD8","TH1", "TH2", "TH17", "NAIVE", "EXHAUSTED", "TFH", "TREG",
              "ACTIVATION", "SIGNALING", "IMMUNE", "PROLIFERATION", "DIVISION", "METABOLISM", "GLYCOLISIS",
              "IL2", "IL4", "IL17", "IL6", "IL21", "IL10", "IFN", "TNF", "IL-2", "IL-4", "IL-17", "IL-6", "IL-21", "IL-10",
              "CYTOKINE", "INTERFERON", "EFF", "MYC", "APOPTOSIS", "TRANSCRIPTION",
              "TF", "INFLAMM", "BCL6", "TBET", "GATA", "ICOS", "PD1", "PD-1", "CXCR5")

pattern_cellfate <- c("TH1", "TH2", "TH17")

pattern_regex <- paste(patterns, collapse = "|")

# finalize df that I want to analyse (use only pathways that have a match for my patterns)
s4_pattern <- s4 %>% mutate(pat = str_detect(gs_name, pattern_regex)) %>% filter(pat == TRUE)


# load own files
filenames <- list.files("output/gamma_fits/feb2021/", pattern="*.csv", full.names=TRUE)
idx <- grepl("summary", filenames)
names <- filenames[idx]
files <- lapply(names, read.csv)


#files <- files[1:2]
# run pathway analysis for each set (nodel / del) for each study with pattern filtering
for(i in seq_along(files)){
  f <- files[[i]]
  set1 <- f %>% filter(best_fit == "gamma") %>% select(gene)
  set2 <- f %>% filter(best_fit == "expo") %>% select(gene)

  # run pathway analysis for current sets and perform fdr correction
  out <- s4_pattern %>% group_by(gs_name) %>% summarize(pval = check_pathway(gene_symbol, set1, set2))
  out$padj <- p.adjust(out$pval, method = "fdr")
  
  # store study name
  study <- str_sub(filenames[i], 35,-1)
  out$study <- study
  
  # store output
  fileout <- paste("pathway", study, sep ="_")
  
  savedir <- "output/pathways/"
  write.csv(out, paste0(savedir, fileout))
}

# conv to regex 