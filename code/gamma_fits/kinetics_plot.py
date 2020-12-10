import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils_plot import *
from utils import *

sns.set(context = "paper", style = "ticks")
fit = pd.read_csv("../../output/gamma_fits/dec2020/gamma_fit_peine_Th0.csv")
data = pd.read_csv("../../data/data_rtm/peine_rtm_Th0_invitro.csv")

# use gene name as index - this should be changed when output is generated
fit = fit.set_index("gene", drop = True)

# look at error distribution
fig,ax = plt.subplots(1,2, figsize = (5,2))
ax[0].hist(fit.rss_x, range = (0,1))
ax[1].hist(fit.rss_y, range = (0,1))
plt.show()

# look at error distribution
fig,ax = plt.subplots(1,3, figsize = (7,2))
ax[0].hist(fit.beta_err_x, range = (0,2))
ax[1].hist(fit.beta_err_y, range = (0,2))
ax[2].hist(fit.alpha_err_y, range = (0,2))
plt.show()


fit_o = fit.loc[fit.sig == "other"]
fit = fit.loc[((fit.rss_x < 1) | (fit.rss_y < 1)) & (fit.alpha_err_y < 2)]

data = data.set_index("gene_name", drop = True)

fit1 = fit.loc[fit.sig != "other"]
fit1["log10padj"] = -np.log10(fit1.padj.values)

g = sns.relplot(data = fit1, x = "alpha_y", y = "log10padj", hue = "sig")
plt.show()

# plot ran
thres = 1.0
fit_g = fit.loc[fit.sig == "Non-exponential"]
fit_e = fit.loc[fit.sig == "Exponential"]

fit_g["mu"] = fit_g.alpha_y.values / fit_g.beta_y.values
fit_e["mu"] = fit_e.alpha_x.values / fit_e.beta_x.values
fit_o["mu"] = np.nan

fig, ax = plt.subplots()
ax.hist(fit_g["mu"].values, range(0,100))
plt.show()

fig, ax = plt.subplots()
ax.hist(fit_e["mu"].values, range(0,100))
#plt.show()

fig, ax = plt.subplots()
ax.hist(fit_g.alpha_y.values, range = (0,10), bins = 30)
#plt.show()

fit_g = fit_g.loc[fit_g.alpha_y > 2]
fit_mu = pd.concat([fit_e, fit_g, fit_o])
fit_mu = fit_mu.loc[((fit_mu.mu > 30) & (fit_mu.mu < 120)) | fit_mu.mu.isnull()]

df = data.merge(fit_mu, how = "left", left_index= True, right_index=True)
df = df.reset_index()

g = sns.relplot(data = df, x = "time", y = "avg_norm_rtm2", ci = "sd",
                col = "sig", kind = "line",  legend = False)
plt.show()


def plot_fit(gene, data, fit):

    # get fit vals but only until 1 is reached
    exp = data.loc[gene]
    x_idx = exp.avg_norm_rtm2.values
    x_idx = np.where(x_idx == 1)[0][0]
    exp = exp.iloc[:(x_idx+1), :]

    alpha_2 = fit.loc[gene].alpha_y
    beta_1 = fit.loc[gene].beta_x
    beta_2 = fit.loc[gene].beta_y

    # get gamma cdf vals
    x_sim = exp.time.values
    x_sim = np.linspace(min(x_sim), max(x_sim), 100)
    y_sim1 = gamma_cdf1(x_sim, beta_1)
    y_sim2 = gamma_cdf(x_sim, alpha_2, beta_2)

    fig, ax = plt.subplots()
    sns.scatterplot(data = exp, x = "time", y = "val_norm_rtm2", ax = ax)
    ax.plot(x_sim, y_sim1)
    ax.plot(x_sim, y_sim2)
    plt.tight_layout()
    plt.show()

#plot_fit("Fbxo36", data, fit)