#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use("../paper_theme_python.mplstyle")

"""
take summary of fit categories and do some postprocessing
first, clean up isoforms
second, filter genes that are expressed at least n times across studies
for those genes, assign "winner" based on how often gene was classified as gamma, expo, bimodal etc
store output
"""

month = "feb2021"
readdir = "../../output/fit_summary/"
df_fits = pd.read_csv(readdir + "category_assignment.csv")

# kick out crawford data set for this analysis because the kinetics are so late?
#df_fits = df_fits[df_fits['study'].str.contains('Crawford')==False]

df_fits[["gene2", "iso1", "iso2"]] = df_fits["gene"].str.split(r"_|\.", expand = True)
df_fits = df_fits.sort_values(by = "gene")
df_fits = df_fits[["gene", "best_fit", "study", "gene2"]]

df_fits.loc[df_fits.best_fit == "longtail", "best_fit"] = "expo"

# get the number of times a genes is found to be kinetic across samples
# note that lambda expr is used coz now multiple isoforms are reduced to gene name
# so study could show up more than once fur each gene
# "nunique" as string also works
df_fits["n_obs"] = df_fits.groupby(["gene2"])["study"].transform(lambda s: s.nunique())

# keep only genes with at least n hits across studies
# this filter seems to affect mainly number of kinetic genes but not the ratios (gamma/expo/etc)
df_fits = df_fits.loc[df_fits["n_obs"] >= 5, :]

# for each gene, find out the highest number of best fit hits
# i.e Eomes --> gamma 4 hits, expo 1 hit --> assign gamma
# problem arises for equal number of hits!
df_win = df_fits.groupby(["gene2", "best_fit"])["study"].nunique()
df_win = df_win.reset_index()

# out = df_win.copy()
# #out["total"] = out.groupby(["gene2"]).transform("sum")


# mygenes = ["Tbx21", "Eomes", "Ifng", "Ifngr", "Gata3",
#            "Ccl5", "Bcl6", "Cxcr5", "Pdcd1", "Rorgt"]

# out = out.loc[out["gene2"].isin(mygenes)]

# out = df_fits.loc[df_fits["gene2"].isin(mygenes)]
# myorder = ["Ccl5", "Bcl6", "Cxcr5", "Eomes", "Gata3",
#            "Ifng", "Tbx21"]

# test = out.groupby(["gene2", "best_fit"]).size().reset_index().pivot
# fig, ax = plt.subplots(figsize = (1.7,1.3))
# g = sns.countplot(y = 'gene2', 
#                   hue = 'best_fit',
#                 data=out,
#                  palette = palette, alpha = 1,
#                  hue_order = ["gamma", "expo", "other", "bimodal"],
#                  order = myorder,
#                  #multiple = 'stack',
#                  #shrink = 0.75
#                  )

# #ax.set_yticklabels(style = "italic")
# ax.set_xlabel("best fit across data sets")
# #ax.set_ylim([0,100])
# #plt.xticks(rotation = 90)

# plt.show()

colors = sns.color_palette("deep", 10)
palette = [colors[0], colors[2], colors[7], "purple"]

# src taken instead of idxmax from stack overflow
df_win = df_win.pivot(columns = "best_fit", values = "study", index = "gene2")

test = df_win.copy()
test = df_win.reset_index()
#%%
df_modules = pd.read_csv("../../genesets_literature/gene_module_summary.csv")

mygenes = ["Ifng", "Il4", "Il21", 
           "Il10","Tnfsf13b", "Eomes", "Ly6c1", "Cxcr3", 
            "Il12rb2", "Ifngr1","Il6ra", "Il2ra",
           "Hlx", "Hopx",  "Ccr4","Cxcr5", "Ccr6", "Ccl5",
          "Tbx21", "Foxp3", "Bcl6", "Rorc", "Gata3"]
test.index = test["gene2"]
test = test.loc[mygenes, ["gamma", "expo", "other", "bimodal"]]

#%%
fig, ax = plt.subplots(figsize = (1.2,3.6))

test.plot(kind = "barh", stacked = True, ax = ax,
          color = palette, width = 0.85, rot =0, legend= False)

ax.set_yticklabels(mygenes, style = "italic")
ax.set_xlabel("no of best fit assignments")
ax.set_ylabel("")
plt.show()
fig.savefig("barplot_categories_signature_genes.svg")
#%%
a = df_win.copy()
df_win['winner'] = a.idxmax(axis=1)
s = a.eq(a.max(axis=1), axis=0).sum(axis=1)
df_win['winner'] = df_win['winner'].mask(s > 1, 'Equality')

# kick out the category that is unclear
df_win = df_win.loc[df_win["winner"] != "Equality"]
df_win = df_win.reset_index()

# save the fit summary
df_win = df_win.rename(columns = {"gene2" : "gene"})
df_win.to_csv("../../output/fit_summary/category_assignment_winners.csv")
