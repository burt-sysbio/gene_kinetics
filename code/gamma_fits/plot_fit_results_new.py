#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:46:15 2020

@author: burt
"""
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sns.set(context = "poster", style = "ticks")
savedir = "../../figures/gamma_fits/"

from scipy.special import gammainc


def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)



def get_data(data, cell, inf):
    """
    load data and add info which alpha has least error
    """
    file = "../../output/gamma_fits/gamma_fit_"+data+"_"+cell+"_"+inf+".csv"

    df_fit_res = pd.read_csv(file)
    df_fit_res = df_fit_res.dropna()
    
    df_alpha = df_fit_res[["alpha1", "alpha2", "alpha10"]]
    best_alpha = df_alpha.apply(find_best_alpha, axis = 1)
    
    err_alpha = df_alpha.min(axis = 1)
    
    df_fit_res["best_alpha"] = best_alpha
    df_fit_res["err_alpha"] = err_alpha
    df_fit_res["infection"] = inf
    df_fit_res["study"] = data
    
    # find minimum of alpha2 and alpha10 err and add this as column to compare vs expo alpha = 1    
    df_dummy = df_fit_res[["alpha2", "alpha10"]]
    minval = df_dummy.min(axis=1).values
    df_fit_res["alpha_2_10"] = minval
    
    # add difference between errors to capture those where differences are large in scatterplot
    err_diff = np.log2(df_fit_res["alpha1"].values / df_fit_res["alpha10"].values)
    df_fit_res["diff_err"] = np.abs(err_diff)
    return(df_fit_res)


def find_best_alpha(row):
    """
    helper fun for get data, check index of minimum in row and add alpha 
    corresponding to it (1,2 or 10)
    """
    min_idx = np.argmin(row)
    alpha = [1,2,10]
    best_alpha = alpha[min_idx]
    return(best_alpha)


# load TFs cytos and receptors
ref_file = "../../gene_sets/references/"
df_tf = pd.read_csv(ref_file+"stubbington_tfs.csv")
df_cyto = pd.read_csv(ref_file+"stubbington_cytokines.csv")
df_receptor = pd.read_csv(ref_file+"stubbington_receptors.csv")


peine1 = ["peine", "Th1", "invitro"]
peine2 = ["peine", "Th2", "invitro"]
peine3 = ["peine", "Th0", "invitro"]
peine4 = ["peine", "ThMix", "invitro"]


crawford1 = ["crawford", "CD4", "arm"]
crawford2 = ["crawford", "CD8", "arm"]
crawford3 = ["crawford", "CD4", "cl13"]
crawford4 = ["crawford", "CD8", "cl13"]


df_list = [peine1, peine2, peine3, peine4, crawford1, crawford2, crawford3, crawford4]
df_list = [get_data(*df) for df in df_list]

# note that when I concatenate here, NAns are induced for crawford data which has fewer timepoints
df_all = pd.concat(df_list)
# visualize the error distribution
#df1 = get_data(*peine1)
df_all["category"] = "other"
df_all["category"][df_all.gene_name.isin(df_tf.gene_name)]="TF"
df_all["category"][df_all.gene_name.isin(df_cyto.gene_name)]="Cytokine"
df_all["category"][df_all.gene_name.isin(df_receptor.gene_name)]="Cyto. Rec."


# split up
df_peine = df_all.loc[df_all.study == "peine"]
df_craw = df_all.loc[df_all.study == "crawford"]


# plotting params
alpha = 0.8
s = 60
lower = 1e-2
upper = 1e1


# visualize error distributions per study
g = sns.relplot(data = df_all, x = "alpha1", y = "alpha10", col = "study", 
                kind = "scatter", alpha = alpha, s = s)
g.set(yscale="log", xscale = "log")
for ax in g.axes.flat:
    ax.plot([0,100], [0,100], '--k') 
g.savefig(savedir+"scatter_fit.pdf")


# zoom in
g = sns.relplot(data = df_all, x = "alpha1", y = "alpha10", col = "study", 
                kind = "scatter", alpha = alpha, s = s)
g.set(yscale="log", xscale = "log", xlim = (lower,upper), ylim = (lower, upper))
for ax in g.axes.flat:
    ax.plot([0,100], [0,100], '--k') 
g.savefig(savedir+"scatter_fit_zoom.pdf")


# kick out genes with no annotation
df3 = df_all.loc[df_all.category != "other"]
g = sns.relplot(data = df3, x = "alpha1", y = "alpha10", col = "category", kind = "scatter")
g.set(yscale="log", xscale = "log", xlim = (lower,upper), ylim = (lower, upper))
for ax in g.axes.flat:
    ax.plot([0,100], [0,100], '--k') 
g.savefig(savedir+"scatter_fit_gene_sets.pdf")


# compare arm vs cl13
df4 = df_craw.loc[df_craw.cell_type == "CD4"]
df4 = df4.loc[(df4.infection == "arm") | (df4.infection == "cl13")]

g = sns.relplot(data = df3, x = "alpha1", y = "alpha10", hue = "category", 
                kind = "scatter", col = "infection")
g.set(yscale="log", xscale = "log", xlim = (lower,upper), ylim = (lower, upper))
for ax in g.axes.flat:
    ax.plot([0,100], [0,100], '--k') 
g.savefig(savedir+"scatter_fit_infection.pdf")


def tidy_df(df, timepoints, alpha, n_genes = 3):
    # keep cols where desired alpha is best fit and of those only good fits   
    df = df[df.best_alpha == alpha]
    thres = 0.1
    df = df[df.err_alpha < thres]
    
    # sort by bad fits for other alpha
    sort_by = "alpha10" if alpha == 1 else "alpha1"
    df = df.sort_values(sort_by, ascending = False)
    df = df.iloc[:n_genes,:]
    # convert to long format for timepoints

    id_vars = ["gene_name", "beta_1", "beta_10"]
    df = df[id_vars + timepoints]
    df = df.melt(id_vars = id_vars, var_name = "time", value_name = "exp_data")
    df = df.dropna()
    df["time"] = pd.to_numeric(df["time"])
    

    return(df)


def plot_rtm_genes(df, n_genes, alpha, savename):
    
    # check that df is only from one study and add corresponding timepoints
    study = df.study.unique()[0]

    assert((study == "crawford") | (study == "peine"))
    if study == "peine":
        timepoints = [0,3,6,12,24,35,48,73,96,120]
    else: 
        timepoints = [0,6,8,15,30]
        
    timepoints = [str(val) for val in timepoints]

    df = tidy_df(df, timepoints, alpha = alpha, n_genes = n_genes)
    
    # plot data
    g = sns.relplot(data = df, x = "time", y = "exp_data", col = "gene_name",
                    kind = "scatter")
    
    # get unique genes and add model fits
    genes = df.gene_name.drop_duplicates()
    genes = genes.values
    for (ax, gene) in zip(g.axes.flat, genes):
        p = df.loc[df.gene_name == gene]
        betas = ["beta_1", "beta_10"]
        alphas = [1, 10]
        colors = ["tab:blue", "tab:red"]
        for alpha, beta, c in zip(alphas, betas, colors):
            # create x array spanning experiment time
            x = np.linspace(p.time.min(), p.time.max(), 100)
            # get one beta value, should all be the same
            b = p[beta].values[0]
    
            y = gamma_cdf(x, alpha, b)
            ax.plot(x,y, ls = "--", color = c)
    
    savedir = "../../figures/gamma_fits/"
    g.savefig(savedir+savename+".pdf")
# focus on specific cell type and infection
my_df = df_craw[(df_craw.infection == "arm") & (df_craw.cell_type == "CD4")]

plot_rtm_genes(my_df, alpha = 10, n_genes = 3, savename = "crawford_rtm")
plot_rtm_genes(my_df, alpha = 1, n_genes = 3, savename = "crawford_expo")

peine_df = df_peine[df_peine.cell_type == "Th1"]

plot_rtm_genes(peine_df, alpha = 10, n_genes = 3, savename = "peine_rtm")
plot_rtm_genes(peine_df, alpha = 1, n_genes = 3, savename = "peine_expo")
