
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

df_fits = df_fits.loc[df_fits.best_fit == "gamma", ["gene", "study"]]

# load some gene lists

# load tfs, cytos, cytoR, th2 th1 genes
sign_dir = "../../gene_sets/references/signaling/"
tf_genes = pd.read_csv(sign_dir + "stubbington_tf.csv", names= ["gene"])

cyto_genes = pd.read_csv(sign_dir + "stubbington_cyto.csv", names= ["gene"])
cytor_genes = pd.read_csv(sign_dir + "stubbington_cytor.csv", names= ["gene"])

df_genes = pd.concat([tf_genes, cyto_genes, cytor_genes])
df_genes = df_genes.drop_duplicates()

# subset fit summary
df = df_fits.loc[df_fits.gene.isin(df_genes.gene)]

# count hits per gene across studies
df2 = df.groupby("gene").count()
df2 = df2.sort_values("study", ascending = False)
df2 = df2.reset_index()
df2["category"] = "TF"
df2.loc[df2.gene.isin(cyto_genes.gene), "category"] = "Cytokine"
df2.loc[df2.gene.isin(cytor_genes.gene), "category"] = "CytoR"


# only look at top 50 genes that are delayed across studies
n = 30
df2 = df2.iloc[:n,:]
df2["ynorm"] = (df2["study"] / 10) * 100

df2.to_csv("../../output/gamma_fits/feb2021/delays_allstudies.csv")

fig,ax = plt.subplots(figsize = (4,9))
yvals = np.arange(1,(n+1))
ax.barh(y = yvals, width = df2.ynorm.values, color = "grey")
plt.yticks(yvals, df2.gene.values)
plt.ylim([0,(n+1)])
plt.xlabel("delayed (% of studies)")
sns.despine()
plt.show()

fig.savefig("../../figures/delays_all_studies.pdf")

# now also plot pvalues of genes that I found that occur in many studies
# read in ftest results
month = "feb2021"
readdir = "../../output/gamma_fits/" + month + "/"
genesetdir = "../../output/genesets/" + month + "/"
filenames = os.listdir(readdir)

pat1 = "ftest"
files1 = [f for f in filenames if pat1 in f]
ftest_all = [pd.read_csv(readdir+ f1) for f1 in files1]
fnames = [f2[6:-4] for f2 in files1]

for f,n in zip(ftest_all,fnames):
    f["study"] = n
    f = f.loc[f.comp == "expo_gamma"]

df_ftests = pd.concat(ftest_all)

# subset genes
df_subset = df2.loc[df2.study >= 5]
df3 = pd.merge(df_ftests, df_subset, on = "gene")
df4 = df3.groupby("gene").min("p")
df4 = df4.reset_index()
df4["log10p"] = -np.log10(df4.p.values)

# make barplot
df4 = df4.sort_values("log10p")
g = sns.catplot(data = df4, x = "log10p", y = "gene", kind = "bar",
                color = "grey", aspect= 0.5)
g.set(xlabel = "-log10(pval)")
plt.show()

g.savefig("../../figures/delays_all_studies_pval.pdf")
