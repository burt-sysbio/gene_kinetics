
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


order = ["crawford_CD4.arm",
         "crawford_CD4.cl13",
         "crawford_CD8.arm",
         "crawford_CD8.cl13",
         "powrie_innate_colitis",
         "proserpio_Th2_parasite",
         "peine_Th0_invitro",
         "peine_Th1_invitro",
         "peine_Th2_invitro",
         "peine_ThMix_invitro",
         "nir_Th0_invitro",
         "nir_Th17_invitro"
         ]

#xlabels = ["[1]CD4-Arm", "[1]CD4-Cl13", "[1]CD8-Arm", "[1]CD8-Cl13",
#           "[2]Th0", "[2]Th17",
#           "[3]Th0", "[3]Th1", "[3]Th2", "[3]Th1/2",
#           "[4]Innate",
#           "[5]Th2"]

xlabels = ["Craw. 1", "Craw. 2", "Craw. 3", "Craw. 4",
           "Nir 1", "Nir 2",
           "Peine 1", "Peine 2", "Peine 3", "Peine 4",
           "Ilott 1",
           "Pros. 1"]

# get total number of kinetic genes
df = df_fits.groupby(["study", "best_fit"])["gene"].count()
df = df.reset_index()
df = df.pivot(index = "study", columns = "best_fit", values = "gene")

total = df.sum(axis = 1)
df_norm = df.values / np.reshape(total.values, (df.shape[0],1))

df_norm = pd.DataFrame(df_norm, columns = df.columns)
df_norm.index = df.index


ax = df.plot.bar(stacked = True)
ax.set_ylabel("n kinetic genes fitted")
ax.set_xticklabels(xlabels)
plt.show()

df_norm = df_norm * 100

colors = sns.color_palette()
colors_reordered = [colors[2], colors[0], colors[1], "tab:grey"]


df_norm.loc[:,"other"] = df_norm.longtail + df_norm.other

plot_longtail = False
if not plot_longtail:
    print("longtail dist not shown?")
    df_norm = df_norm[["expo", "gamma", "other"]]
    savename = "fit_assignment_reduced"
    colors_reordered = [colors[2], colors[0], "tab:grey"]
else:
    savename = "fit_assignment"

ax = df_norm.plot.bar(stacked = True, color = colors_reordered)
ax.set_ylabel("fit assignment (% kinetic genes)")
ax.set_xticklabels(xlabels)
ax.set_ylim(0,100)
ax.set_xlabel("")
plt.show()

fig = ax.get_figure()
fig.savefig("../../figures/"+ savename + ".pdf")
fig.savefig("../../figures/"+ savename + ".svg")


# save the fit summary
df_fits.to_csv("../../output/gamma_fits/feb2021/fit_summary_all.csv")
