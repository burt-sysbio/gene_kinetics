
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
sns.set(style = "ticks", context = "poster")
from statsmodels.stats.multitest import fdrcorrection
import os

def filter_error(df):
    thres = 0.9
    thres_beta = df.beta_err.quantile(thres)
    thres_rss = df.rss.quantile(thres)

    df["keep_fit"] = False
    df.loc[(df.rss < thres_rss) & (df.beta_err < thres_beta), ["keep_fit"]] = True

    # only for gamma fit there is an error for alpha so use this but only then
    if not np.isnan(df.alpha_err).any():
        thres_alpha = df.alpha_err.quantile(thres)
        df.loc[df.alpha_err < thres_alpha, ["keep_fit"]]  = True

    return df


# read in all files
month = "jan2021"
readdir = "../../output/gamma_fits/" + month + "/"
filenames = os.listdir(readdir)

pattern1 = "fit_expo"
pattern2 = "fit_gamma"
pattern3 = "ftest"
files1 = [f for f in filenames if pattern1 in f]
files2 = [f for f in filenames if pattern2 in f]
files3 = [f for f in filenames if pattern3 in f]

readdir = "../../output/gamma_fits/" + month + "/"

for f1,f2,f3 in zip(files1, files2, files3):
    # load data
    df_exp = pd.read_csv(readdir + f1)
    df_gamma = pd.read_csv(readdir + f2)
    df_fit = pd.read_csv(readdir + f3)


    #df = pd.concat([df_exp, df_gamma])

    #df = df.dropna()
    # # look at residual squares distribution

    # plot error distributions after filtering
    df_exp = filter_error(df_exp)
    df_gamma = filter_error(df_gamma)

    #df_filt = pd.concat([df_exp, df_gamma])
    #g = sns.displot(data = df_filt, x = "rss", hue = "model")
    #plt.show()

    #g = sns.displot(data = df_filt, x = "beta_err", hue = "model")
    #plt.show()

    #g = sns.displot(data = df_filt, x = "alpha_err", hue = "model")
    #plt.show()

    # get genes that have good fit parameters for gamma or exponential fit
    good_fit_idx = df_exp["keep_fit"] & df_gamma["keep_fit"]
    df_fit["keep_fit"] = False
    df_fit.loc[good_fit_idx, ["keep_fit"]] = True

    # get pvals for the good fits
    pvals = df_fit.loc[good_fit_idx, ["p"]].values

    # add fdr correction but only for those genes that have a good fit statistic
    df_fit["padj"] = np.nan
    df_fit["sig_delay"] = False
    padj = fdrcorrection(pvals.flatten())[1]
    df_fit.loc[good_fit_idx, ["padj"]] = padj
    df_fit.loc[df_fit["padj"] <= 0.05, ["sig_delay"]] = True

    # save output, use f3 and subtract prefix words
    df_fit.to_csv(readdir + "fit_summary_" + f3[6:])


