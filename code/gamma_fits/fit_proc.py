
import pandas as pd
from utils import run_f_test
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
sns.set(style = "ticks", context = "poster")

def filter_error(df):
    thres_beta = df.beta_err.quantile(0.9)
    thres_rss = df.rss.quantile(0.9)

    df = df.loc[(df.rss < thres_rss) & (df.beta_err < thres_beta)]
    if not np.isnan(df.alpha_err).any():
        thres_alpha = df.alpha_err.quantile(0.9)
        df = df.loc[df.alpha_err < thres_alpha]

    return df


study = "peine"
groups = ["Th0"]
activations = ["invitro"]
month = "jan2021"

readdir = "../../output/gamma_fits/" + month + "/"

for gr in groups:
    for act in activations:
        dummy = study + "_" + gr + "_" + act + ".csv"
        file_a = readdir + "fit_expo_" + dummy
        file_b = readdir + "fit_gamma_" + dummy
        file_c = readdir + "ftest_" + dummy

        # load data
        df_exp = pd.read_csv(file_a)
        df_gamma = pd.read_csv(file_b)
        df_fit = pd.read_csv(file_c)
        df = pd.concat([df_exp, df_gamma])

        df = df.dropna()
        # # look at residual squares distribution

        # plot error distributions after filtering
        df_exp_f = filter_error(df_exp)
        df_gamma_f = filter_error(df_gamma)

        df_filt = pd.concat([df_exp_f, df_gamma_f])
        #g = sns.displot(data = df_filt, x = "rss", hue = "model")
        #plt.show()

        #g = sns.displot(data = df_filt, x = "beta_err", hue = "model")
        #plt.show()

        #g = sns.displot(data = df_filt, x = "alpha_err", hue = "model")
        #plt.show()

        # get genes that have good fit parameters for gamma or exponential fit
        genes_good_fit = df_filt.gene.drop_duplicates()

        df_fit_good = df_fit[df_fit.gene.isin(genes_good_fit)]
        df_fit_good.to_csv(readdir + "fit_summary_" + dummy)


