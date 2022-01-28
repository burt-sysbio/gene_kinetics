import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt
sns.set(context = "poster", style = "ticks")
# read in all files
month = "feb2021"
readdir = "../../output/gamma_fits/" + month + "/"
genesetdir = "../../output/genesets/" + month + "/"
filenames = os.listdir(readdir)

pat1 = "fit_res"
pat2 = "ftest"

files1 = [f for f in filenames if pat1 in f]
files2 = [f for f in filenames if pat2 in f]

fits_all = [pd.read_csv(readdir+ f1) for f1 in files1]
ftest_all = [pd.read_csv(readdir+ f2) for f2 in files2]
fnames = [f2[6:-4] for f2 in files2]

for f,n in zip(fits_all,fnames):
    f["study"] = n

df_dist = pd.concat(fits_all)

xlabels = ["Craw. 1", "Craw. 2", "Craw. 3", "Craw. 4",
           "Nir 1", "Nir 2",
           "Peine 1", "Peine 2", "Peine 3", "Peine 4",
           "Ilott 1",
           "Pros. 1"]

palette = sns.color_palette()
palette_reordered = [palette[0], palette[2], palette[1]]
g = sns.catplot(data = df_dist, x = "study", y = "rmse", hue = "model", kind = "box", aspect= 2,
                fliersize = 0, hue_order= ["gamma", "expo", "longtail"], whis=[5, 95], legend_out=True,
                palette = palette_reordered)

g.set(ylim = (0,0.7), ylabel = "RMSE (a.u.)", xlabel = "")
g.set_xticklabels(rotation=90, labels= xlabels)

savedir = "../../figures/fit_error_distributions/"
g.savefig(savedir + "rmse_dist.pdf")
g.savefig(savedir + "rmse_dist.svg")

plt.tight_layout()
plt.show()


g = sns.displot(data = df_dist, x = "rmse", hue = "model", palette = palette_reordered,
                hue_order= ["gamma", "expo", "longtail"], legend = False, alpha = 0.5, aspect = 1.2)
g.set(xlim = (0,0.6), ylim = (0,1700), xlabel = "fit error (RMSE)", ylabel = "n genes")
sns.despine(top = False, right = False)
g.savefig(savedir + "rmse_dist_comb.pdf")
g.savefig(savedir + "rmse_dist_comb.svg")
plt.show()