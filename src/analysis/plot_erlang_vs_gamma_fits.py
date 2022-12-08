
import seaborn as sns
import pandas as pd
from src.analysis.utils_plot import plot_single_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gammainc
import random
plt.style.use("../paper_theme_python.mplstyle")


def gamma_cdf(t, alpha, beta):
    dummy = t * beta
    return gammainc(alpha, dummy)

def gamma_cdf1(t, beta):
    return gamma_cdf(t, 1, beta)

def fit_gamma(x, y, sigma):
    # dummy output in case fit does not work (e.g. gene kinetics are bimodal)
    # run fit catch runtime exception for bad fit
    fit_val, fit_err = curve_fit(f=gamma_cdf, xdata=x, ydata=y, sigma=sigma,
                                 absolute_sigma=True)
    # error of fitted params (not used atm)
    err = np.sqrt(np.diag(fit_err))

    # assign fit values depending on gamma_fun (exponential fit only fits beta not alpha)
    alpha_fit, beta_fit = fit_val
    alpha_err, beta_err = err

    # compute chi sqr and get residuals
    # need to change pipeline function in R to write sigma as 1 instead of NA
    nexp = gamma_cdf(x, *fit_val)
    rss = np.sum((y-nexp)**2)
    # root mean squared
    rmse = np.sqrt(rss/len(y))

    out = [alpha_fit, beta_fit, rss, rmse, alpha_err, beta_err]
    return out

# read in all files
month = "feb2021"
fitdir = "../../output/gamma_fits/" + month + "/"
datadir = "../../data/data_rtm/"

df_modules = pd.read_csv("../../genesets_literature/gene_module_summary.csv")
#df_modules = df_modules.loc[df_modules["module"].isin(["Cytokines", "Cytokine Receptor"])]

study = "Proserpio_rtm_Th2_parasite.csv"
study2 = study[:-4]

fit_res_gamma = pd.read_csv(fitdir + "fit_res_" + study)
fit_res_gamma = fit_res_gamma.loc[fit_res_gamma["model"] == "gamma"]
fit_res_gamma = fit_res_gamma.loc[fit_res_gamma["gene"].isin(df_modules["gene"])]
fit_res_gamma = fit_res_gamma.loc[fit_res_gamma["rmse"]<0.3]
fit_res_gamma = fit_res_gamma.loc[fit_res_gamma["alpha"]<30]
fit_res_gamma = fit_res_gamma.loc[fit_res_gamma["alpha"]>=2]
fit_res_gamma = fit_res_gamma.loc[fit_res_gamma["alpha_err"]<1]
fit_res_gamma = fit_res_gamma.loc[fit_res_gamma["beta_err"]<1]



data = pd.read_csv(datadir + study)


def myfit(df):
    df = df[["time", "avg_norm_rtm2", "SD"]].drop_duplicates().copy()
    x = df["time"].values
    y = df["avg_norm_rtm2"].values
    sigma = df["SD"].values

    out = fit_gamma(x,y, sigma)
    return out

def plot_fit(data, gene, ax):
    df = data.loc[data.gene == gene]

    # run fit
    alpha_fit, beta_fit, rss, rmse, alpha_err, beta_err = myfit(df)

    # based on fit, compute corresponding Erlang alpha and beta
    mu_fit = alpha_fit/beta_fit
    alpha_round = round(alpha_fit)
    beta_round = alpha_round / mu_fit

    x = df["time"].values
    y = df["avg_norm_rtm2"].values
    SD = df["SD"].values

    nexp = gamma_cdf(x, alpha_round, beta_round)
    rss_round = np.sum((y-nexp)**2)
    # root mean squared
    rmse_round = np.sqrt(rss_round/len(y))

    x_arr = np.linspace(x[0],x[-1],100)
    y_arr_gamma = gamma_cdf(x_arr, alpha_fit, beta_fit)
    y_arr_erlang = gamma_cdf(x_arr, alpha_round, beta_round)

    ax.plot(x_arr, y_arr_gamma, c = "k", lw = 1.8)
    ax.plot(x_arr, y_arr_erlang, c = "grey", lw = 1.8)

    ax.text(20,0.1,f"E({alpha_round},{round(beta_round,2)}), rmse={round(rmse_round,2)}",size = 6)
    ax.text(20,0.,f"g({round(alpha_fit,2)},{round(beta_fit,2)}), rmse={round(rmse,2)}", size = 6)
    ax.text(0,0.9, gene, size = 8, style = "italic")
    ax.errorbar(x,y, yerr =SD, fmt = ".")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("")
    ax.set_yticks([0,0.5,1])
    ax.set_yticklabels(["","",""])
    ax.set_xticks([0,24,48,72,96])
    ax.set_ylim([-0.1,1.1])

fig, axes = plt.subplots(1,6, figsize = (8.5,1.2))
gene_list = list(fit_res_gamma.gene.sample(10))

gene_selection = ["Nkg7", "Gata3", "Il4", "Il21", "Icos", "Ly6c2"]


for gene, ax in zip(gene_selection, axes):
    plot_fit(data, gene, ax)

axes[0].set_ylabel("expr. norm.")
axes[0].set_yticklabels(["0", "0.5", "1.0"])
#plt.tight_layout()
plt.show()

fig.savefig("../../figures/gamma_vs_erlang/fit_proserpio_gamma_vs_erlang.pdf", transparent = True)
fig.savefig("../../figures/gamma_vs_erlang/fit_proserpio_gamma_vs_erlang.svg", transparent = True)

