
import numpy as np
import pandas as pd
import seaborn as sns
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

sns.set(style = "ticks", context = "poster")

month = "jan2021"
readdir = "../../output/gamma_fits/" + month + "/"

filenames = os.listdir(readdir)
pattern = "fit_summary"
filenames = [f for f in filenames if pattern in f]
filepaths = [readdir + f for f in filenames]

fits = []
for f,n in zip(filepaths,filenames):
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

fig, ax = plt.subplots(figsize = (12,12))
sns.countplot(data = df_fits, y = "study", hue = "f-test", ax = ax, order = order)
plt.show()

fig.savefig("../../figures/fit_summary_datasets.pdf")

df_fits.to_csv(readdir + "fit_summary_all.csv")