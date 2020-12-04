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
