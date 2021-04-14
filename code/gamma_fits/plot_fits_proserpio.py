
import seaborn as sns
import pandas as pd
import numpy as np
from utils_plot import plot_single_fit
import matplotlib.pyplot as plt
sns.set(context = "poster", style = "ticks")


# load tfs, cytos, cytoR, th2 th1 genes
cells_dir = "../../gene_sets/references/effector_states/"
sign_dir = "../../gene_sets/references/signaling/"

th2_genes = list(pd.read_csv(cells_dir + "stubbington_th2.txt", sep=" ", header=None).values.flatten())
th1_genes = list(pd.read_csv(cells_dir + "stubbington_th1.txt", sep=" ", header=None).values.flatten())
th17_genes = list(pd.read_csv(cells_dir + "stubbington_th17.txt", sep=" ", header=None).values.flatten())

cyto_genes = list(pd.read_csv(sign_dir + "stubbington_cyto.csv").values.flatten())
cytor_genes = list(pd.read_csv(sign_dir + "stubbington_cytor.csv").values.flatten())
tf_genes = list(pd.read_csv(sign_dir + "stubbington_tf.csv").values.flatten())


# read in all files
month = "feb2021"
fitdir = "../../output/gamma_fits/" + month + "/"
datadir = "../../data/data_rtm/"

study = "Proserpio_rtm_Th2_parasite.csv"
fFit = "fit_res_" + study
fTest = "ftest_" + study
fGenes = "fit_summary_" + study

fit_res = pd.read_csv(fitdir + fFit)
data = pd.read_csv(datadir + study)
ftest = pd.read_csv(fitdir + fTest)
df_genes = pd.read_csv(fitdir + fGenes)


#if "peine" in study:
#    data = data.rename(columns={"gene_name": "gene"})

fit_res = fit_res.set_index(["gene", "model"])
data = data.set_index("gene")

show_fit = "other"
genes = df_genes.loc[df_genes["best_fit"] == show_fit, ["gene"]].values
genes = genes.flatten()

genes = np.random.choice(genes, 9)

fig, axes = plt.subplots(3,3, figsize = (12,10))
axes = axes.flatten()
for gene, ax in zip(genes,axes):
    plot_single_fit(gene, data, fit_res, ax, capsize = 10)

plt.tight_layout()
plt.show()


# check if I have anything to show in each category

th2_sig = df_genes.loc[df_genes.gene.isin(cyto_genes)]


genes = ["Tnfsf4", "Socs2", "Cntf", "Il4"]
fig, axes = plt.subplots(1,4, figsize = (20,4))
axes = axes.flatten()
for gene, ax in zip(genes,axes):
    plot_single_fit(gene, data, fit_res, ax, show_rmse = False)

plt.tight_layout()
plt.show()


for gene in genes:
    fig, ax = plt.subplots()
    plot_single_fit(gene, data, fit_res, ax, capsize = 10, show_rmse = False, show_title = True)
    plt.tight_layout()
    plt.show()
#fig.savefig("../../figures/example_fits_proserpio.pdf")
#fig.savefig("../../figures/example_fits_proserpio.svg")

# double-check f-test for specific gene
#gene = "Bub1b"
#n = len(df2.time)
#rss1 = fit_res.loc[gene, "expo"].rss
#rss2 = fit_res.loc[gene, "gamma"].rss

#pval = f_test(rss1, rss2, n)


test = data.loc["Emb"]
