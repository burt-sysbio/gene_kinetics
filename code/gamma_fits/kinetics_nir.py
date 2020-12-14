import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils_plot import *
from utils import *

sns.set(context = "talk", style = "ticks")
fit = pd.read_csv("../../output/gamma_fits/dec2020/gamma_fit_nir_Th0.csv")
data = pd.read_csv("../../data/data_rtm/nir_rtm_Th0.csv")


# look at residual squares distribution
g = sns.displot(data = fit, x = "rss", col = "model")
plt.show()

# look err distribution
#g = sns.displot(data = fit, x = "beta_err", col = "model")
#g.set(xscale = "log")
#plt.show()


# look err distribution
#fit_gamma = fit.loc[fit.model == "gamma"]
#g= sns.displot(data = fit_gamma, x = "alpha_err")
#g.set(xscale = "log")
#plt.show()

# filter for low error and residual squares
fit = fit.loc[(fit.rss < 0.5) &
              ((fit["alpha_err"] < 2) | (np.isnan(fit["alpha_err"]))) &
              (fit["beta_err"] < 2)]

# make volcona plot of pval vs alpha
fit["log10padj"] = -np.log10(fit.padj.values)

# plot distribution of mean arrivals
fit["mean_arrival"] = fit.alpha.values / fit.beta.values
#sns.displot(data = fit, x = "mean_arrival", hue = "model")
#plt.show()


# plot volcano
fit_gamma = fit.loc[fit.model == "gamma"]
g = sns.relplot(data = fit_gamma, x = "alpha", y = "log10padj", hue = "f-test")
g.set(xscale = "log")
plt.show()

# plot alpha distribution
#g = sns.displot(data = fit_gamma, x = "alpha")
#g.set(xscale = "log")
#plt.show()

# merge fit results with data and plot normalized arrival times for expo and gamma fits
fit_merge = fit[["gene", "model", "f-test", "mean_arrival"]]
df = fit_merge.merge(data, how = "left", on = "gene")
df["time_norm"] = df.time.values / df.mean_arrival.values

data_gamma = df.loc[(df.model != "expo") & (df["f-test"] == "sig")]
data_expo = df.loc[(df.model == "expo") & (df["f-test"] == "ns")]
df = pd.concat([data_expo, data_gamma])
#g = sns.relplot(data = df, x = "time_norm",y = "avg_norm_rtm2",
#                col = "model", legend = False)
#g.set(xlim = (0,10))
#plt.show()

df2 = df.loc[(df.gene == "Nudt1") | (df.gene == "Actn1")]
g = sns.relplot(data = df2, x = "time_norm", y = "avg_norm_rtm2",
                col = "model", hue = "gene", legend = False)
plt.show()

# get time when max is reached
time_idx = df.loc[df.avg_norm_rtm2 == 1.0][["time", "gene"]]
time_idx = time_idx.rename(columns = {"time" : "timemax"})

df = df.merge(time_idx, how = "left", on = "gene")
df = df.loc[df.time <= df.timemax]
df["time_norm2"] = df.time / df.timemax

g = sns.relplot(data = df, x = "time_norm2", y = "avg_norm_rtm2",
                kind = "line", col = "model", legend = False)
plt.show()

# normalize time to mean arrival time
#df = data.merge(fit, how = "left", left_index= True, right_index=True)
#df = df.reset_index()
#g = sns.relplot(data = df, x = "time", y = "avg_norm_rtm2", ci = "sd",
#                col = "sig", kind = "line",  legend = False)
#g.set(ylabel = "expr. norm.")
#plt.show()

#plot_fit("Fbxo36", data, fit)

# t = np.linspace(0,10, 100)
# x1 = gamma_cdf(t, 1,1)
# x2 = gamma_cdf(t, 1,2)
# x3 = gamma_cdf(t, 1, 0.5)
#
# fig,ax = plt.subplots()
# ax.plot(t, x1)
# ax.plot(t, x2)
# ax.plot(t, x3)
# plt.show()
#
# y1 = t
# y2 = t/0.5
# y3 = t/2
#
#
# fig,ax = plt.subplots()
# ax.plot(y1, x1)
# ax.plot(y2, x2)
# ax.plot(y3, x3)
# plt.show()
