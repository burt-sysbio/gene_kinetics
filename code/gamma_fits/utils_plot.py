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
from utils import *


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


def plot_fit(gene : str, data, fit):
    """
    purpose: plot gene expression time course from data with corresponding fits
    data: df from data_rtm processed data
    fit : df from gamma fits generated from data
    """
    # get fit vals but only until 1 is reached
    exp = data.loc[gene]
    x_idx = exp.avg_norm_rtm2.values
    x_idx = np.where(x_idx == 1)[0][0]
    exp = exp.iloc[:(x_idx+1), :]

    alpha_2 = fit.loc[gene].alpha_y
    beta_1 = fit.loc[gene].beta_x
    beta_2 = fit.loc[gene].beta_y

    # get gamma cdf vals
    x_sim = exp.time.values
    x_sim = np.linspace(min(x_sim), max(x_sim), 100)
    y_sim1 = gamma_cdf1(x_sim, beta_1)
    y_sim2 = gamma_cdf(x_sim, alpha_2, beta_2)

    fig, ax = plt.subplots()
    sns.scatterplot(data = exp, x = "time", y = "val_norm_rtm2", ax = ax)
    ax.plot(x_sim, y_sim1)
    ax.plot(x_sim, y_sim2)
    plt.tight_layout()
    plt.show()