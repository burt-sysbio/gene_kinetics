import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(context = "poster", style = "ticks")

dir_pw = "../../gene_sets/references/stubbington_"
tfs = pd.read_csv(dir_pw + "tfs.csv")
cyto = pd.read_csv(dir_pw + "cytokines.csv")
cytoR = pd.read_csv(dir_pw + "receptors.csv")

dir_fit = "../../output/gamma_fits/jan2021/fit_summary_all.csv"
df = pd.read_csv(dir_fit)

df["pathway"] = "other"

idx_tf = df.gene.isin(tfs.gene_name)
idx_cyto = df.gene.isin(cyto.gene_name)
idx_cytoR = df.gene.isin(cytoR.gene_name)

df.loc[idx_tf, ["pathway"]] = "TF"
df.loc[idx_cyto, ["pathway"]] = "Cyto"
df.loc[idx_cytoR, ["pathway"]] = "CytoR"

order = ["TF", "Cyto", "CytoR"]

df_red = df[(df.pathway != "other") & (df.keep_fit == True)]

g = sns.catplot(data = df_red, x = "pathway", hue = "f-test", col = "study",
                col_wrap= 3, kind = "count", order = order, hue_order= ["ns", "sig"])
g.set_titles("{col_name}")
plt.show()

g.savefig("../../figures/barplots_pathways.pdf")

df_sum = df.groupby(["study", "pathway", "f-test"])["gene"].count()
df_sum2 = df.groupby(["study", "f-test"])["gene"].count()

df_sum = df_sum.to_frame()
df_sum = df_sum.reset_index()
