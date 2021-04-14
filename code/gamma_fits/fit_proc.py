import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set(style="ticks", context="talk")
from statsmodels.stats.multitest import fdrcorrection
import os

# def myfun(fit_res, ftest):
#     fit1 = filter_error(fit1)
#
#     good_fit = fit1["keep_fit"] & fit2["keep_fit"]
#
#     ftest["keep_fit"] = False
#     ftest["padj"] = np.nan
#     ftest["sig_padj"] = "ns"
#     ftest.loc[good_fit, ["keep_fit"]] = True
#     pvals = ftest.loc[good_fit, ["p"]].values
#     ftest.loc[good_fit, ["padj"]] = fdrcorrection(pvals.flatten())[1]
#     ftest.loc[ftest["padj"] <= 0.05, ["sig_padj"]] = "sig"
#
#     return fit1, fit2, ftest

# def filter_error(df):
#     thres = 0.9
#     thres_beta = df.beta_err.quantile(thres)
#     thres_rss = df.rss.quantile(thres)
#
#     df["keep_fit"] = False
#     df.loc[(df.rss < thres_rss) & (df.beta_err < thres_beta), ["keep_fit"]] = True
#
#     # only for gamma fit there is an error for alpha so use this but only then
#     if not np.isnan(df.alpha_err).any():
#         thres_alpha = df.alpha_err.quantile(thres)
#         df.loc[df.alpha_err < thres_alpha, ["keep_fit"]] = True
#
#     return df


def filter_error(df):
    rmse_thres = 0.2
    df["keep_fit"] = False
    df.loc[df["rmse"] < rmse_thres, ["keep_fit"]] = True
    return df


def get_genesets(df_fit):
    # take dataframe after filtering for errors and return gene sets for pathway analysis
    col = "f-test"
    val1 = "sig"

    # get genes where fit was bad
    df0 = df_fit.loc[df_fit["keep_fit"] == False, ["gene"]].drop_duplicates()
    df0["best_fit"] = "other"

    # get genes where fit was good
    df_fit = df_fit.loc[df_fit["keep_fit"] == True]

    # get number of genes where delay was significant
    df1 = df_fit.loc[(df_fit[col] == val1) & (df_fit["comp"] == "expo_gamma"), ["gene"]]
    df2 = df_fit.loc[(df_fit[col] == val1) & (df_fit["comp"] == "expo_longtail"), ["gene"]]

    # expo are all genes with good fit except those in df1 or df2
    df3 = df_fit.loc[~(df_fit.gene.isin(df1.gene) | df_fit.gene.isin(df2.gene)), ["gene"]]
    df3 = df3.drop_duplicates()

    df1["best_fit"] = "gamma"
    df2["best_fit"] = "longtail"
    df3["best_fit"] = "expo"

    df = pd.concat([df0, df1, df2, df3])
    return df


# read in all files
month = "feb2021"
readdir = "../../output/gamma_fits/" + month + "/"
genesetdir = "../../output/genesets/" + month + "/"
filenames = os.listdir(readdir)

pat1 = "fit_res"
pat2 = "ftest"

files1 = [f for f in filenames if pat1 in f]
files2 = [f for f in filenames if pat2 in f]

#files1 = [files1[0]]
#files2 = [files2[0]]

fits_all = [pd.read_csv(readdir+ f1) for f1 in files1]
ftest_all = [pd.read_csv(readdir+ f2) for f2 in files2]
fnames = [f2[6:-4] for f2 in files2]

for f,n in zip(fits_all,fnames):
    f["study"] = n

df_dist = pd.concat(fits_all)

#g = sns.displot(data = df_dist, x = "rmse", hue = "model", row = "study",
#                alpha = 0.5,
#                aspect = 1)
#plt.show()

xlabels = ["[1]CD4-Arm", "[1]CD4-Cl13", "[1]CD8-Arm", "[1]CD8-Cl13",
           "[2]Th0", "[2]Th17",
           "[3]Th0", "[3]Th1", "[3]Th2", "[3]Th1/2",
           "[4]Innate",
           "[5]Th2"]

palette = sns.color_palette()
palette_reordered = [palette[0], palette[2], palette[1]]
g = sns.catplot(data = df_dist, x = "study", y = "rmse", hue = "model", kind = "box", aspect= 2,
                fliersize = 0, hue_order= ["gamma", "expo", "longtail"], whis=[5, 95], legend_out=True,
                palette = palette_reordered)

g.set(ylim = (0,0.7), ylabel = "RMSE (a.u.)", xlabel = "")
g.set_xticklabels(rotation=90, labels= xlabels)
g.savefig("../../figures/rmse_dist.pdf")
g.savefig("../../figures/rmse_dist.svg")

plt.tight_layout()
plt.show()


g = sns.displot(data = df_dist, x = "rmse", hue = "model", palette = palette_reordered,
                hue_order= ["gamma", "expo", "longtail"], aspect= 1.5)
g.set(xlim = (0,0.8), ylim = (0,1700), xlabel = "fit error", ylabel = "n obs")

g.savefig("../../figures/rmse_dist_comb.pdf")
g.savefig("../../figures/rmse_dist_comb.svg")
plt.show()



for fit_res, ftest, fname in zip(fits_all, ftest_all, fnames):

    fit_res = filter_error(fit_res)

    # get genes that have good fit parameters for gamma  longtail or exponential fit
    genes_good_fit = fit_res.loc[fit_res["keep_fit"], ["gene"]]
    ftest["keep_fit"] = False
    ftest.loc[ftest.gene.isin(genes_good_fit.gene), "keep_fit"] = True

    # save gene sets with delayed and nondelayed genes separately for pathway analysis
    df = get_genesets(ftest)
    df.to_csv(readdir + "fit_summary_" + fname + ".csv")
