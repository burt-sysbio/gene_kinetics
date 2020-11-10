#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:40:08 2020

@author: burt
"""

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sns.set(context = "poster", style = "ticks")

df_fit_res = pd.read_csv("../output/kinetics_gamma_fit.csv")
df_fit_res = df_fit_res.dropna()

df_kinetics = pd.read_csv("../output/kinetics_tidy.csv")

# load gene sets for filtering
df_th1 = pd.read_csv("../gene_sets/references/stubbington_th1.csv")
df_receptor = pd.read_csv("../gene_sets/references/stubbington_receptors.csv")
df_cytos = pd.read_csv("../gene_sets/references/stubbington_cytokines.csv")
df_tf = pd.read_csv("../gene_sets/references/stubbington_tfs.csv")

df_list = [df_kinetics, df_th1, df_tf, df_cytos, df_receptor]
savenames = ["all", "th1", "tfs", "cytokines", "receptors"]

# do analysis with and wo error threshold
err_list = [np.inf, 1]

df_out = [[],[]]
err_dist = [[],[]]

for err, out, err_dis in zip(err_list, df_out, err_dist):
    # max error is either 1 or inf
    max_err = err
      
    for df, name in zip(df_list, savenames):
        df_fit = df_fit_res[df_fit_res["gene_name"].isin(df["gene_name"])]
        df_kin = df_kinetics[df_kinetics["gene_name"].isin(df["gene_name"])]
        
        # visualize error distribution
        df_err = df_fit.melt(id_vars = ["infection", "cell_type", "gene_name"])
        df_err["chisqr"] = np.log10(df_err.value)
        
     
        df_err["name"] = name
        err_dis.append(df_err)
     
        crit1 = df_fit.alpha1 < max_err
        crit2 = df_fit.alpha2 < max_err
        crit3 = df_fit.alpha10 < max_err
        crit = crit1 | crit2 | crit3
        
        df_fit = df_fit[crit]
        
        # find number of genes where respective alpha is best fit
        n_alpha1 = sum((df_fit.alpha1 < df_fit.alpha2) & (df_fit.alpha1 < df_fit.alpha10))
        n_alpha2 = sum((df_fit.alpha2 < df_fit.alpha1) & (df_fit.alpha2 < df_fit.alpha10))
        n_alpha10 = sum((df_fit.alpha10 < df_fit.alpha1) & (df_fit.alpha10 < df_fit.alpha2))
        
        x = [r"$\alpha=1$", r"$\alpha=2$", r"$\alpha=10$"]
        y = [n_alpha1, n_alpha2, n_alpha10]
        
        df = pd.DataFrame({"best_fit":x, "n_genes":y})
        

        df["name"] = name
        out.append(df)

        
        if name != "all":   
            g = sns.relplot(data = df_kin, 
                            x = "time", 
                            hue = "gene_name", 
                            y = "log2FC_mean", 
                            row = "infection",
                            col = "cell_type", 
                            kind = "line", 
                            facet_kws= {"margin_titles" : True})
            
            g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
            g.savefig("../figures/kinetic_genes_"+name+".png")
        
df1 = pd.concat(df_out[0])
df2 = pd.concat(df_out[1])

df = pd.concat(err_dist[0])


# compute relative contributions of alphas to total number of kinetic genes
df1["n"] = df1.groupby(["name"])["n_genes"].transform("sum")
df1["n_rel"] = df1.n_genes/df1.n

df2["n"] = df2.groupby(["name"])["n_genes"].transform("sum")
df2["n_rel"] = df2.n_genes/df2.n


g = sns.catplot(data = df1, x = "best_fit", y = "n_rel", hue = "name", kind = "bar", aspect = 1.5)
#g.legend.set_title("")
g.savefig("../figures/gamma_fits2.png")

g = sns.catplot(data = df2, x = "best_fit", y = "n_rel", hue = "name", kind = "bar", aspect = 1.5)
#g.legend.set_title("")

g.savefig("../figures/gamma_fits_err2.png")


g = sns.catplot(data = df1, x = "best_fit", y = "n_genes", hue = "name", kind = "bar", aspect = 1.5)
#g.legend.set_title("")
g.savefig("../figures/gamma_fits_total2.png")

g = sns.catplot(data = df2, x = "best_fit", y = "n_genes", hue = "name", kind = "bar", aspect = 1.5)
#g.legend.set_title("")

g.savefig("../figures/gamma_fits_err_total2.png")


g = sns.catplot(data = df, x = "variable", y = "chisqr", hue = "name", kind = "box")
g.set(xlabel = "", ylabel = "log10(chisqr)")
#g.legend.set_title("")
g.savefig("../figures/gamma_fits_err_dist2.png")

