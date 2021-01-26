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

df_red = df[df.pathway != "other"]
g = sns.catplot(data = df_red, x = "pathway", hue = "f-test", col = "study",
                col_wrap= 3, kind = "count", hue_order= ["ns", "sig"])
g.set_titles("{col_name}")
plt.show()

g.savefig("../../figures/barplots_pathways.pdf")

df_sum = df.groupby(["study", "pathway", "f-test"])["gene"].count()
df_sum2 = df.groupby(["study", "f-test"])["gene"].count()

df_sum = df_sum.to_frame()
df_sum = df_sum.reset_index()

df_sum2 = df_sum2.to_frame()
df_sum2 = df_sum2.reset_index()

df_merge = pd.merge(df_sum, df_sum2, on = ["study", "f-test"])
df_merge["h_pathway"] = df_merge["gene_x"] / df_merge["gene_y"]


df_merge_red = df_merge[df_merge.pathway != "other"]
g = sns.catplot(data = df_merge_red, x = "pathway", y = "h_pathway",
                hue = "f-test", col = "study", col_wrap= 3, kind = "bar",
                hue_order= ["ns", "sig"], facet_kws= {"sharey" : False})
g.set_titles("{col_name}")
plt.show()

g.savefig("../../figures/barplots_rel_freq_pathways.pdf")