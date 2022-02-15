"""
load fit results and match with gene modules
needs thorough cleaning of gene annotations to match with pathways

take fit information from individual genes and compute avg+SEM per module and study
"""

from scipy.special import gammainc
import numpy as np

def gamma_cdf(t, alpha, beta):
    dummy = t * beta
    return gammainc(alpha, dummy)

import pandas as pd
import seaborn as sns
import os
from scipy.stats import gamma
import matplotlib.pyplot as plt
import warnings
sns.set(context = "paper", style = "ticks")
warnings.warn("time unit might not be adjusted between data sets")

df_modules = pd.read_csv("../../genesets_literature/gene_module_summary.csv")

# kick out some categories
rm_module = ["Th17", "Differentiation"]
df_modules = df_modules.loc[~df_modules["module"].isin(rm_module)]
# only keep genes with good gamma fits
rmse_thres = 0.2

ylim = [None,None]
xlim = [None,None]
capsize = 4
linewidth = 1.5
xlabel = "avg (h)"
ylabel = "SD (h)"
modelfit = "gamma"


df = pd.read_csv("../../output/fit_summary/fit_res_all.csv")

# convert beta because unit is days for Crawford and Powrie
df.loc[df["study"].str.contains("Craw|Powrie", regex =True), "beta"] = df["beta"]/24

df = df.loc[df["rmse"] < rmse_thres]
df = df.loc[df["model"] == modelfit]
#df = df.loc[df["alpha_err"]<10]
#df = df.loc[df["beta_err"]<10]

# clean gene names
df[["gene2", "iso1", "iso2"]] = df["gene"].str.split(r"_|\.", expand = True)

df.drop(columns = ["gene"], inplace = True)
df.rename(columns = {"gene2" : "gene"}, inplace = True)

df = pd.merge(df, df_modules, how = "inner", on = "gene")

df["mu"] = df["alpha"] / df["beta"]
df["SD"] = np.sqrt(df["alpha"] / df["beta"]**2)

df_exvivo = df.loc[df["study"].str.contains("Peine_rtm_Th2|Pros")]

g = sns.catplot(data = df_exvivo, x = "study", y = "mu", kind = "bar", ci = 68)
plt.show()
g = sns.catplot(data = df_exvivo, x = "study", y = "SD", kind = "bar", ci =68)
plt.show()

#df_exvivo = df_exvivo.loc[df_exvivo["module"]!="Th1"]
# kick out Craford and Powrie because of timescale

# only focus on Th1 study and compare with Th2 study
df = df.loc[df["study"].str.contains("Peine")]
df = df.loc[df["study"].str.contains("Th1|Th2", regex = True)]

# compute mean and sem for genes across studies and modules
out = df.groupby(["module", "study"]).agg({'mu': ['mean', 'sem'], 'SD': ["mean", "sem"],
                                           'alpha': ['mean', 'sem'], 'beta': ["mean", "sem"]})
out.reset_index(inplace = True)
out.columns = ["module", "study", "avg", "avg_err", "SD", "SD_err", "alpha", "alpha_err", "beta", "beta_err"]

g = sns.FacetGrid(data = out, col = "study", hue = "module", legend_out=True, aspect = 0.7)
g.map_dataframe(plt.errorbar, x= "avg", y = "SD", yerr = "SD_err", xerr= "avg_err",elinewidth=linewidth, capsize=capsize)
g.set_titles("{col_name}")
g.add_legend()
g.set(xlim = xlim, ylim = ylim, xlabel = xlabel, ylabel = ylabel)
plt.show()


#g.savefig("../../figures/module_quantification/th1_th2_condition.pdf")
#g.savefig("../../figures/module_quantification/th1_th2_condition.svg")


out_th1 = out.loc[out["study"] == "Peine_rtm_Th1_invitro"]
mymodules = ["Cytokine Receptor", "Cytokines", "Transcription Factor", "Proliferation", "Th1"]
out_th1 = out_th1.loc[out_th1["module"].isin(mymodules)]
g = sns.FacetGrid(data = out_th1, hue = "module",  aspect = 1.1, legend_out=False,
                  hue_order= ["Th1", "Proliferation","Transcription Factor","Cytokines","Cytokine Receptor"],
                  height = 2.1)
g.map_dataframe(plt.errorbar, x= "avg", y = "SD", yerr = "SD_err", xerr= "avg_err",elinewidth=linewidth, capsize=capsize)
g.add_legend()
g.set(xlim = xlim, ylim = ylim, xlabel = xlabel, ylabel = ylabel)
sns.despine(top = False, right = False)
plt.show()


g.savefig("../../figures/module_quantification/peine_th1cell_modules.pdf")
g.savefig("../../figures/module_quantification/peine_th1cell_modules.svg")


#
df_th1 = df.loc[(df["study"] == "Peine_rtm_Th1_invitro") & (df["module"]=="Th1")]

g = sns.relplot(data = df_th1, x = "alpha", y = "beta",
                kind = "scatter")
g.set(xlim = [0,100], ylim = [0,2])
plt.show()

mydf = out_th1.loc[out_th1["module"] == "Th1"]

time = np.linspace(0,96,100)

alpha = mydf["alpha"].values[0]
alpha_err = mydf["alpha_err"].values[0]
beta = mydf["beta"].values[0]
beta_err = mydf["beta_err"].values[0]

y1 = gamma_cdf(time, alpha-alpha_err, beta-beta_err)
y2 = gamma_cdf(time, alpha-alpha_err, beta+beta_err)
y3 = gamma_cdf(time, alpha+alpha_err, beta-beta_err)
y4 = gamma_cdf(time, alpha+alpha_err, beta+beta_err)
y5 = gamma_cdf(time, alpha, beta)

mylist = [y2,y3]

fig, ax = plt.subplots(figsize = (2.5,2.2))
for y in mylist:
    ax.plot(time, y, color="tab:blue", lw=0.5)
ax.plot(time, y5, color = "tab:blue", lw = 2)

ax.fill_between(time, y2, y3, color = "tab:blue", alpha = 0.1)
ax.set_xlim([0,time[-1]])
ax.set_xlabel("time (h)")
ax.set_ylabel("expr. norm.")
ax.set_xticks([0,24,48,72,96])
ax.set_title("Th1 differentiation onset")
ax.set_ylim([0,1])
plt.show()

fig.savefig("../../figures/module_quantification/th1_module_timecourse.svg")

#########################################################################
#########################################################################
#########################################################################
# plot avg and SD distr as barplot + SEM (ci=68)
# for data in [df, df_exvivo]:
#     for y in ["mu", "SD"]:
#         g = sns.catplot(data = data, x = "module", y = y, aspect = 0.8,
#                         hue = "study", kind = "bar", ci = 68, capsize = 0.1)
#         g.set(xlabel = "")
#         g.set_titles("{col_name}")
#         g.set_xticklabels(rotation = 90)
#         sns.despine(top = False, right = False)
#         plt.show()
#


# plot the actual gamma distributions
# alpha_values = out["alpha"].values
# beta_values = out["beta"].values
# theta_values = 1/beta_values
# x = np.linspace(0, 96, 200)
#
# gamma_list = []
# for alpha, theta in zip(alpha_values, theta_values):
#     dist = gamma(alpha, 0, theta)
#     gamma_vals = dist.pdf(x)
#     gamma_df = pd.DataFrame({"time": x, "value" : gamma_vals})
#     gamma_df["alpha"] = alpha
#     gamma_list.append(gamma_df)
#
# gamma_df = pd.concat(gamma_list)
# gamma_df = pd.merge(gamma_df, out, how = "left", on = "alpha")
#
# g = sns.relplot(data = gamma_df, x = "time", y = "value", hue = "module", col = "study",
#                 kind = "line")
# plt.show()