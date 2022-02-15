
import seaborn as sns
import pandas as pd
import numpy as np
from code.analysis.utils_plot import plot_single_fit
import matplotlib.pyplot as plt
sns.set(context = "paper", style = "ticks")

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


fit_res = fit_res.set_index(["gene", "model"])
data = data.set_index("gene")

show_fit = "other"
genes = df_genes.loc[df_genes["best_fit"] == show_fit, ["gene"]].values
genes = genes.flatten()

genes = np.random.choice(genes, 9)

fig, axes = plt.subplots(3,3, figsize = (12,10))
axes = axes.flatten()
for gene, ax in zip(genes,axes):
    plot_single_fit(gene, data, fit_res, ax)

plt.tight_layout()
plt.show()

