
import numpy as np
import pandas as pd
import seaborn as sns
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

sns.set(style = "ticks", context = "paper")

month = "jan2021"
readdir = "../../output/gamma_fits/" + month + "/"

filenames = os.listdir(readdir)
pattern = "fit_summary"
filenames = [f for f in filenames if pattern in f]
filepaths = [readdir + f for f in filenames]

fits = []
for f,n in zip(filepaths,filenames):
    #
    if "fit_summary_all" not in n:
        fit = pd.read_csv(f)
        fit["study"] = n[12:-4]
        fits.append(fit)

df_fits = pd.concat(fits)


order = ["crawford_CD4_arm",
         "crawford_CD8_arm",
         "crawford_CD4_cl13",
         "crawford_CD8_cl13",
         "powrie_innate_colitis",
         "proserpio_Th2_parasite",
         "peine_Th0_invitro",
         "peine_Th1_invitro",
         "peine_Th2_invitro",
         "peine_ThMix_invitro",
         "nir_Th0_invitro",
         "nir_Th17_invitro"
         ]

use_fdr = False
if use_fdr:
    col = "sig_delay"
    val1 = True
    val2 = False
else:
    col = "f-test"
    val1 = "sig"
    val2 = "ns"

# get total number of kinetic genes
n_kin = df_fits.groupby(["study"])["gene"].count()

# get number of genes where fit was good
df_good_fit = df_fits.loc[df_fits["keep_fit"] == True]
n_good_fit = df_good_fit.groupby(["study"])["gene"].count()

# get number of genes where delay was significant (with fdr correction)
df_sig_del = df_good_fit.loc[df_good_fit[col] == val1]
n_sig_del = df_sig_del.groupby(["study"])["gene"].count()

# get n genes where delay was not
df_ns_del = df_good_fit.loc[df_good_fit[col] == val2]
n_ns_del = df_ns_del.groupby(["study"])["gene"].count()

df_sum = pd.concat([n_kin, n_good_fit, n_sig_del, n_ns_del], axis = 1)
df_sum.columns = ["kinetic_genes", "good_fits", "sig_delay", "ns_delay"]
df_sum["bad_fits"] = df_sum.kinetic_genes - df_sum.good_fits

df_sum = df_sum.drop(columns = ["kinetic_genes", "good_fits"])
ax = df_sum.plot.bar(stacked = True)
plt.show()

# normalize data frame to sum everything up to 1, to do this divide by all kinetic genes
df_norm = df_sum.values / np.reshape(n_kin.values, (n_kin.shape[0],1))
df_norm = pd.DataFrame(df_norm, columns = df_sum.columns)
df_norm.index = df_sum.index

ax = df_norm.plot.bar(stacked = True)
plt.show()

#fig.savefig("../../figures/fit_summary_datasets.pdf")

df_fits.to_csv(readdir + "fit_summary_all.csv")