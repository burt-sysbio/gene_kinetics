import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from utils import gamma_cdf
"""
plot fit results for any gene module and study together with actual data
current script uses SEM, can be replaced with SD
"""

sns.set(context = "poster", style = "ticks")
data = pd.read_csv("../../data/data_summary/data_rtm_gene_module.csv")

ftest = pd.read_csv("../../data/data_summary/module_ftest_results.csv")
fit_results = pd.read_csv("../../data/data_summary/module_fit_results.csv")

mymodule = "th1_genes"
myID = "Peine_rtm_Th1_invitro"

data = data.loc[(data["ID"] == myID) & (data["module"] == mymodule)]
ftest = ftest.loc[(ftest["ID"] == myID) & (ftest["gene"] == mymodule)]
fit_results = fit_results.loc[(fit_results["ID"] == myID) & (fit_results["gene"] == mymodule)]


data["SEM"] = data.groupby(["time"])["val_norm"].transform("sem")
data = data[["time", "avg_norm2", "SEM"]].drop_duplicates()

x_idx = data.avg_norm2.values
x_idx = np.where(x_idx == 1)[0][0]
data = data.iloc[:(x_idx + 1), :]

alphas = fit_results["alpha"].values
betas = fit_results["beta"].values

# get gamma cdf vals
data_x = data.time.drop_duplicates().values
x_sim = np.linspace(0, 120, 100)
print("xsim max was set manually to 120 for Peine data")
y_sim = [gamma_cdf(x_sim, alpha, beta) for alpha, beta in zip(alphas, betas)]

y_sim_gamma = y_sim[1]

fig, ax = plt.subplots()
#sns.scatterplot(data=data, x="time", y="avg_norm2", ax=ax, color="k")
ax.plot(x_sim, y_sim_gamma, color = "tab:blue", ls = "--", lw = 3)
ax.errorbar(data.time, data["avg_norm2"], yerr=data["SEM"], ecolor="tab:blue", fmt=".", capsize=6)
ax.set_xlim([x_sim.min(), x_sim.max()])
ax.set_xticks([0,24,48,72,96,120])
ax.set_xlabel("time (h)")
ax.set_ylabel("% of maximum")
plt.show()

#fig.savefig("../../figures/tcell_module/peine_th1_module_timecourse.pdf")
#fig.savefig("../../figures/tcell_module/peine_th1_module_timecourse.svg")
