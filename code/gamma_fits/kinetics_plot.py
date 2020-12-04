import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils_plot import *

sns.set(context = "poster", style = "ticks")
fit = pd.read_csv("../../output/gamma_fits/dec2020/gamma_fit_peine_Th0.csv")
data = pd.read_csv("../../output/rtm_data/peine_rtm_Th0_invitro.csv")

# use gene name as index - this should be changed when output is generated
fit = fit.set_index("gene", drop = True)
data = data.set_index("gene_name", drop = True)

# subset data with transcription factors
tfs = pd.read_csv("../../gene_sets/references/stubbington_tfs.csv")

# make scatter plot for alpha and beta fits and indicate significance
fit = fit.loc[fit["sig"] == True]
fit = fit.loc[fit["alpha_err_y"] < 3]

tf_list = list(tfs.gene_name.values)
fit_tf = fit.loc[fit.index.isin(tf_list)]


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

plot_fit("Fbxo36", data, fit)