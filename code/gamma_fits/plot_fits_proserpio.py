
import seaborn as sns
import pandas as pd
import numpy as np
from utils_plot import plot_single_fit
import matplotlib.pyplot as plt


def plot_multi_genes(genes, savename, data , fit_res_gamma , fit_res_gmm, ncol = 2, nrow =3, figsize = (4, 5),
                     style = "ticks", context = "paper", make_tight = False, **kwargs):
    sns.set(context=context, style=style)
    fig, axes = plt.subplots(nrow, ncol, figsize= figsize)
    axes = axes.flatten()
    for i in range(len(genes)):
        gene = genes[i]
        ax = axes[i]
        plot_single_fit(gene, data, fit_res_gamma, fit_res_gmm, ax, **kwargs)
        ax.set_xticks([0,50,100])
        #ax.set_yticks([0, 0.5, 1.0])
        if make_tight:
            if i != 0:
                ax.set_ylabel("")
                ax.set_yticklabels([])

    plt.tight_layout()
    plt.show()
    fig.savefig("../../figures/example_fits_proserpio/timecourse_fits_proserpio_" + savename +".pdf")
    fig.savefig("../../figures/example_fits_proserpio/timecourse_fits_proserpio_" + savename + ".svg")

# read in all files
month = "feb2021"
fitdir = "../../output/gamma_fits/" + month + "/"
datadir = "../../data/data_rtm/"

study = "Proserpio_rtm_Th2_parasite.csv"
study2 = study[:-4]

data = pd.read_csv(datadir + study)
fit_res_gmm = pd.read_csv(fitdir + "fit_GMM_" + study)
fit_res_gamma = pd.read_csv(fitdir + "fit_res_" + study)
df_categories = pd.read_csv("../../output/fit_summary/category_assignment.csv")

fit_res_gamma["study"] = study2
fit_res_gamma = pd.merge(fit_res_gamma, df_categories, on = ["gene", "study"], how = "left")

fit_res_gamma = fit_res_gamma.set_index(["gene", "model"])
fit_res_gmm = fit_res_gmm.set_index("gene")
data = data.set_index("gene")

# check if I have anything to show in each category
#genes = ["Tnfsf4", "Socs2", "Cntf", "Il4"]
#genes = ["Sdf2l1", "Alox12", "Olfr56", "Ktn1", "Upp1"]
genes = ["Ifng", "Tbx21", "Eomes", "Gata3", "Il4"]
candidate_categories = [["expo"], ["gamma"], ["longtail"], ["bimodal"], ["expo"]]
# other expo candidates: "Gm15558"
fig, axes = plt.subplots(1,5,figsize = (8,2))
axes = axes.flatten()
for i in range(len(genes)):
    gene = genes[i]
    ax = axes[i]
    my_cats = candidate_categories[i]
    plot_single_fit(gene, data, fit_res_gamma, fit_res_gmm, ax, my_categories= my_cats, show_rmse=False, lw = 3)
    ax.set_xticks([0,50,100])
    ax.set_yticks([0, 0.5, 1.0])

    if i != 0:
        ax.set_ylabel("")
        ax.set_yticklabels([])

plt.tight_layout()
plt.show()
#fig.savefig("../../figures/example_fits_proserpio/timecourse_fits_proserpio_candidate_genes.pdf")
#fig.savefig("../../figures/example_fits_proserpio/timecourse_fits_proserpio_candidate_genes.svg")

# plot some genes from the other category
fittypes = ["bimodal", "gamma", "expo", "longtail", "other"]

# sample some genes for each category and plot it
# for fit in fittypes:
#     df = df_categories.loc[(df_categories["best_fit"] == fit) & (df_categories["study"] == study2), :]
#     n_genes = 6
#     genes = df["gene"].drop_duplicates().sample(n_genes).tolist()
#
#     l1 = fit_res_gmm.index.get_level_values(0).tolist()
#     l2 = fit_res_gamma.index.get_level_values(0).tolist()
#
#     genes = [x for x in genes if (x in l1) & (x in l2)]
#     plot_multi_genes(genes, fit, data, fit_res_gamma, fit_res_gmm, my_categories = ["gamma", "expo", "longtail", "bimodal"],
#                      show_rmse = False, lw = 2)


#genes_bad = ["1700016C15Rik", "Akap6", "Aqp4", "Ccl7", "Cda", "Mok", "Pkp1", "Osbpl10"]
#genes_no_success = ["Parp12", "Pgls", "Rpl9", "Ass1", "Zfp131", "Poglut2"]
#genes_bimodal = ["Pam16", "Gatb", "Psrc1", "Hey1", "Jag1", "Pja1", "Lpar1"]

#plot_multi_genes(genes_bad, "high_error", data , fit_res_gamma , fit_res_gmm)
#plot_multi_genes(genes_no_success, "fit_unsuccessful", data , fit_res_gamma , fit_res_gmm)
#plot_multi_genes(genes_bimodal, "bimodal_pos", data , fit_res_gamma , fit_res_gmm)
