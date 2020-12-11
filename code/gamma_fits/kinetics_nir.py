import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils_plot import *
from utils import *

sns.set(context = "talk", style = "ticks")
fit = pd.read_csv("../../output/gamma_fits/dec2020/gamma_fit_nir_Th0.csv")
data = pd.read_csv("../../data/data_rtm/nir_rtm_Th0.csv")

# use gene name as index - this should be changed when output is generated
fit = fit.set_index("gene", drop = True)
data = data.set_index("gene", drop = True)

# look at residual squares distribution
ylim = 700
r = (0,0.5)
fig,ax = plt.subplots(1,2, figsize = (6,3))
ax[0].hist(fit.rss_x, range = r)
ax[0].set_ylim(0,ylim)
ax[0].set_xlabel("fit alpha=1")
ax[0].set_ylabel("fit error")
ax[1].hist(fit.rss_y, range = r)
ax[1].set_xlabel("fit alpha!=1")
ax[1].set_ylim(0,ylim)
plt.show()

# look at error distribution
r_err = (0,2)
fig,ax = plt.subplots(1,3, figsize = (7,2))
ax[0].hist(fit.beta_err_x, range = r_err)
ax[0].set_xlabel("err beta x")
ax[1].hist(fit.beta_err_y, range = r_err)
ax[1].set_xlabel("err beta y")
ax[2].hist(fit.alpha_err_y, range = r_err)
ax[2].set_xlabel("err alpha y")
plt.show()


# filter for low error and residual squares
fit_o = fit.loc[fit.sig == "other"]
fit = fit.loc[((fit.rss_x < 1) | (fit.rss_y < 1)) & (fit.alpha_err_y < 2)]

# make volcona plot of pval vs alpha
fit1 = fit.loc[fit.sig != "other"]
fit1["log10padj"] = -np.log10(fit1.padj.values)
g = sns.relplot(data = fit1, x = "alpha_y", y = "log10padj", hue = "sig")
plt.show()

# plot distribution of mean arrivals
fit_g = fit.loc[fit.sig == "Non-exponential"]
fit_e = fit.loc[fit.sig == "Exponential"]

fit_g["mu"] = fit_g.alpha_y.values / fit_g.beta_y.values
fit_e["mu"] = fit_e.alpha_x.values / fit_e.beta_x.values
fit_o["mu"] = np.nan


r_mu = (0,100)
fig, ax = plt.subplots(1, 2)
ax[0].hist(fit_g["mu"].values)
ax[0].set_xlabel("mu x")
ax[1].hist(fit_e["mu"].values)
ax[1].set_xlabel("mu y")
plt.show()

# plot distribution of alpha values for non exponential fit
r_alpha = (0,10)
fig, ax = plt.subplots()
ax.hist(fit_g.alpha_y.values)
ax.set_xlabel("alpha non Exponential")
plt.show()

# filter non exponential fits for high alphas and focus on distinct mean arrival range
mu_min = 30
mu_max = 120
alpha_min = 2
fit_g = fit_g.loc[fit_g.alpha_y > alpha_min]
fit_mu = pd.concat([fit_e, fit_g, fit_o])
fit_mu = fit_mu.loc[((fit_mu.mu > mu_min) & (fit_mu.mu < mu_max)) | fit_mu.mu.isnull()]

# plot timecourse of gene expression for exponential non exponential and other category
df = data.merge(fit_mu, how = "left", left_index= True, right_index=True)
df = df.reset_index()
g = sns.relplot(data = df, x = "time", y = "avg_norm_rtm2", ci = "sd",
                col = "sig", kind = "line",  legend = False)
g.set(ylabel = "expr. norm.")
plt.show()

#plot_fit("Fbxo36", data, fit)