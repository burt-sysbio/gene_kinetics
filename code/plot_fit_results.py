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

df_fit_res = pd.read_csv("../output/kinetics_gamma_fit.csv")
df_fit_res = df_fit_res.dropna()

df_kinetics = pd.read_csv("../output/kinetics_tidy.csv")

# load gene sets for filtering
df_th1 = pd.read_csv("../gene_sets/references/stubbington_th1.csv")
df_receptor = pd.read_csv("../gene_sets/references/stubbington_receptors.csv")
df_cytos = pd.read_csv("../gene_sets/references/stubbington_cytokines.csv")

df_fit_res = df_fit_res[df_fit_res["gene_name"].isin(df_receptor["gene_name"])]
df_kinetics = df_kinetics[df_kinetics["gene_name"].isin(df_receptor["gene_name"])]


# visualize error distribution
df_err = df_fit_res.melt(id_vars = ["infection", "cell_type", "gene_name"])
df_err["errlog"] = np.log10(df_err.value)
g = sns.catplot(data = df_err, x = "variable", y = "errlog", kind = "box")
g.set(ylabel = "log(chisqr)", xlabel = "")
# overall error
vals = df_err.value.values


# kick out all values where error is not below maxerr in one condition
maxerr = np.inf

crit1 = df_fit_res.alpha1 < maxerr
crit2 = df_fit_res.alpha2 < maxerr
crit3 = df_fit_res.alpha10 < maxerr
crit = crit1 | crit2 | crit3

df_fit_res = df_fit_res[crit]

# find number of genes where respective alpha is best fit
n_alpha1 = sum((df_fit_res.alpha1 < df_fit_res.alpha2) & (df_fit_res.alpha1 < df_fit_res.alpha10))
n_alpha2 = sum((df_fit_res.alpha2 < df_fit_res.alpha1) & (df_fit_res.alpha2 < df_fit_res.alpha10))
n_alpha10 = sum((df_fit_res.alpha10 < df_fit_res.alpha1) & (df_fit_res.alpha10 < df_fit_res.alpha2))

x = [r"$\alpha=1$", r"$\alpha=2$", r"$\alpha=10$"]
y = [n_alpha1, n_alpha2, n_alpha10]

df = pd.DataFrame({"best_fit":x, "n_genes":y})

fig, ax = plt.subplots()
ax = sns.barplot(data = df, x = "best_fit", y = "n_genes", color = "tab:blue")
plt.tight_layout()


g = sns.relplot(data = df_kinetics, 
                x = "time", 
                y = "expr_norm", 
                hue = "gene_name", 
                row = "infection",
                col = "cell_type", 
                kind = "line", 
                facet_kws= {"margin_titles" : True}, 
                legend= False)

g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.savefig("../figures/kinetic_genes_all_FC.pdf")

g = sns.relplot(data = df_kinetics, 
                x = "time", 
                hue = "gene_name", 
                y = "log2FC_mean", 
                row = "infection",
                col = "cell_type", 
                kind = "line", 
                facet_kws= {"margin_titles" : True}, 
                legend= False)

g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
#g.savefig("../figures/kinetic_genes_all_gamma_dist.pdf")