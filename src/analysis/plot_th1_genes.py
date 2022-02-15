import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import gammainc
import numpy as np

def gamma_cdf(t, alpha, beta):
    dummy = t * beta
    return gammainc(alpha, dummy)

# read in all files
month = "feb2021"
fitdir = "../../output/gamma_fits/" + month + "/"
datadir = "../../data/data_rtm/"

study = "Peine_rtm_Th1_invitro.csv"
study2 = study[:-4]

data = pd.read_csv(datadir + study)

fit_res = pd.read_csv(fitdir + "fit_res_" + study)
fit_res = fit_res.loc[fit_res["model"]=="gamma"]
fit_res[["gene2", "iso1", "iso2"]] = fit_res["gene"].str.split(r"_|\.", expand = True)

df_categories = pd.read_csv("../../output/fit_summary/category_assignment.csv")

df_categories = df_categories.loc[df_categories["best_fit"] == "gamma"]
df_categories = df_categories.loc[df_categories["study"] == study2]

df_modules = pd.read_csv("../../genesets_literature/gene_module_summary.csv")
df_modules = df_modules.loc[df_modules["module"] == "Th1"]

fig, ax = plt.subplots()
time = np.linspace(0,120,120)
mygenes = ["Ifng", "Tbx21", "Eomes", "Ly6c1"]

for gene in df_modules["gene"].values:
    mydf = fit_res.loc[fit_res["gene2"] == gene]

    if mydf.shape[0] > 0:
        mydf = mydf.loc[mydf["rmse"] == mydf["rmse"].min()]

        if mydf["gene"].values[0] in df_categories["gene"].values:
            color = "k"
            lw = 0.1
        if mydf["gene"].values[0] in mygenes:
            color = "tab:blue"
            lw = 3
        else:
            color = "k"
            lw = 0.1
        yarr = gamma_cdf(time, mydf["alpha"].values, mydf["beta"].values)
        ax.plot(time, yarr, color = color, lw = lw)

plt.show()

