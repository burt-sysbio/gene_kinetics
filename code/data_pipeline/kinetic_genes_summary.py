
import numpy as np
import pandas as pd
import seaborn as sns
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

sns.set(style = "ticks", context = "talk")


month = "feb2021"
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
df_fits = df_fits.sort_values("gene")

# kick out crawford data set for this analysis because the kinetics are so late?
df_fits = df_fits[df_fits['study'].str.contains('Crawford')==False]

# get the number of times a genes is found to be kinetic across samples (studies and celltypes)
df_fits["n_obs"] = df_fits.groupby(["gene"])["study"].transform("count")

# keep only genes that are kinetic in at least 3 studies
df_fits = df_fits.loc[df_fits["n_obs"] >= 3, :]

# refactor gene names, at the moment there are different gene isoforms both indicated by "." and by "_"
#df_fits[['gene1','split1']] = df_fits['gene'].str.split("_", expand=True)
df_fits[['gene2','split2', "split3"]] = df_fits['gene'].str.split(r"\.", expand=True)


df_out = df_fits["gene2"].drop_duplicates()
df_out.to_csv("../../output/delayed_genes/background_kinetic_all_studies.txt", index = False, header = False)

