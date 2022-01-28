
import numpy as np
import pandas as pd
import seaborn as sns
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

sns.set(style = "ticks", context = "talk")


month = "feb2021"
readdir = "../../output/gamma_fits/" + month + "/"

filenames = os.listdir(readdir)
pattern = "fit_summary"
filenames = [f for f in filenames if pattern in f]
filepaths = [readdir + f for f in filenames]

fits = []
for f,n in zip(filepaths,filenames):
    #
    if "fit_summary_all" not in n:
        fit = pd.read_csv(f)
        fit["study"] = n[12:-4]
        fits.append(fit)

df_fits = pd.concat(fits)

# kick out crawford data set for this analysis because the kinetics are so late?
#df_fits = df_fits[df_fits['study'].str.contains('Crawford')==False]


# get the number of times a genes is found to be kinetic across samples (studies and celltypes)
df_fits["n_obs"] = df_fits.groupby(["gene"])["study"].transform("count")

# keep only genes that are kinetic in at least 3 studies

use_crit = True
fit_type = "gamma"

if use_crit:
    df_fits = df_fits.loc[df_fits["n_obs"] >= 3, :]
    critname = "_withcrit"
else:
    critname = "_default"
# reduce to only relevant hits for gamma best fit
df_fits["fit_type"] = 0
df_fits.loc[df_fits["best_fit"] == fit_type, "fit_type"] = 1
df_fits["n_obs_gamma"] = df_fits.groupby(["gene"])["fit_type"].transform("sum")

# get percentage of studies where gamma was best fit
df_fits["rel_gamma"] = 100 * df_fits["n_obs_gamma"] / df_fits["n_obs"]

# kick out information from which study the data is?
df_fits2 = df_fits[["gene", "n_obs", "n_obs_gamma", "rel_gamma"]]
df_fits2 = df_fits2.drop_duplicates()

# for some reason all gamma genes do not have "_" as sep
#df_fits2[['gene']] = df_fits2['gene'].str.split("_", expand=True)
df_fits2[['gene','split2', "split3"]] = df_fits2['gene'].str.split(r"\.", expand=True)

# kick out unused separator columns, note that this df may still contain duplicates bc isoforms occured with different obs
df_fits2 = df_fits2[["gene", "n_obs", "n_obs_gamma", "rel_gamma"]].drop_duplicates()

# another sanity check: for those genes where diff. isoforms had uneven n obs keep the higher one
df_fits2 = df_fits2.sort_values('n_obs_gamma').drop_duplicates('gene', keep='last')
df_fits2 = df_fits2.sort_values("rel_gamma")

# assign deg list based on threshold, remove leftover duplicates
df_deg = df_fits2.loc[df_fits2["rel_gamma"] > 80]["gene"].drop_duplicates()
df_deg.to_csv("../../output/delayed_genes/delayed_top80.txt", index = False, header = False)

df_deg = df_fits2.loc[df_fits2["rel_gamma"] == 0]["gene"].drop_duplicates()
df_deg.to_csv("../../output/delayed_genes/delayed_min10.txt", index = False, header = False)

df_gsea = df_fits2[["gene", "rel_gamma"]]
df_gsea = df_gsea.drop_duplicates(subset = "gene")
df_gsea.to_csv("../../output/delayed_genes/gsealist_ranked_relgamma.txt", index = False, header = False, sep = "\t")
# also save as rnk file
df_gsea.to_csv("../../output/delayed_genes/gsealist_ranked_relgamma.rnk", index = False, header = False, sep = "\t")

# visualize no of kinetic genes across studies and

# load tfs, cytos, cytoR, th2 th1 genes
sign_dir = "../../../genesets_literature/references/signaling/"
tf_genes = pd.read_csv(sign_dir + "stubbington_tf.csv", names= ["gene"])

cyto_genes = pd.read_csv(sign_dir + "stubbington_cyto.csv", names= ["gene"])
cytor_genes = pd.read_csv(sign_dir + "stubbington_cytor.csv", names= ["gene"])

#df_genes = pd.concat([tf_genes, cyto_genes, cytor_genes])
#df_genes = df_genes.drop_duplicates()
df_genes = pd.read_csv("../../../genesets_literature/references/genes_caro.csv", names = ["gene"])

df_th1 = pd.read_csv("../../../genesets_literature/references/effector_states/stubbington_th1.txt", names = ["gene"])
df_th2 = pd.read_csv("../../../genesets_literature/references/effector_states/stubbington_th2.txt", names = ["gene"])
df_th17 = pd.read_csv("../../../genesets_literature/references/effector_states/stubbington_th17.txt", names = ["gene"])

# subset fit summary

def barplot_delay(df, geneset, fname, show_top = 30):
    df = df.loc[df.gene.isin(geneset)]
    # count hits per gene across studies
    df = df.reset_index(drop =True)
    df = df.sort_values("rel_gamma", ascending= False)
    # only look at top 50 genes that are delayed across studies
    n = df.shape[0] ### replace this with n = x for subset of genes
    if show_top is not None:
        n = show_top
    df_top = df.iloc[:n,:]

    fig,ax = plt.subplots(figsize = (4,9))
    yvals = np.arange(1,(n+1))
    ax.barh(y = yvals, width = df_top["rel_gamma"].values, color = "grey")
    plt.yticks(yvals, df_top.gene.values)
    plt.ylim([0,(n+1)])
    plt.xlabel("delayed (% of studies)")
    sns.despine()
    plt.show()

    fig.savefig("../../figures/barplot_genelists_delays/delays_barplot_" + fname + critname + fit_type +".pdf")

    fig, ax = plt.subplots()
    sns.histplot(data = df, x = "rel_gamma", stat = "probability")
    ax.set_xlabel("delayed in % of studies")
    ax.set_ylabel("freq." + fname)
    ax.set_ylim([0,0.3])
    plt.show()
    fig.savefig("../../figures/barplot_genelists_delays/delays_histplot_" + fname + critname + fit_type +".pdf")

#barplot_delay(df_fits2, df_genes["gene"], "CD4_expert_list")
#barplot_delay(df_fits2, tf_genes["gene"], "TF_genes")
#barplot_delay(df_fits2, cyto_genes["gene"], "cyto_genes")
#barplot_delay(df_fits2, cytor_genes["gene"], "cytoR_genes")
#barplot_delay(df_fits2, df_th1["gene"], "th1_genes")
#barplot_delay(df_fits2, df_th2["gene"], "th2_genes")
#barplot_delay(df_fits2, df_th17["gene"], "th17_genes")


fig, ax = plt.subplots()
sns.histplot(df_fits2,x = "n_obs_gamma", bins= 10)
ax.set_xlabel("delayed in n studies")
ax.set_ylabel("n genes")

if use_crit:
    title = str(df_fits2.shape[0]) + " kinetic genes in at least 3 studies"
else:
    title = str(df_fits2.shape[0]) + " kinetic genes (no filter)"
ax.set_title(title)
plt.show()

fig.savefig("../../figures/histo_relgamma" + critname + fit_type + ".pdf")

