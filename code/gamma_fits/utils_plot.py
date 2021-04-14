#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:46:15 2020

@author: burt
"""
import seaborn as sns
import pandas as pd
import numpy as np
from utils import gamma_cdf


def prep_fit(gene : str, data, fit_res):
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

    alphas = fit_res.loc[gene].alpha.values
    betas = fit_res.loc[gene].beta.values

    # get gamma cdf vals
    x_sim = exp.time.values
    x_sim = np.linspace(0, max(x_sim)+20, 100)
    y_sim = [gamma_cdf(x_sim, alpha, beta) for alpha, beta in zip(alphas, betas)]

    df_sim = pd.DataFrame({"time" : x_sim, "fit1" : y_sim[0], "fit2" : y_sim[1], "fit3" : y_sim[2]})
    return df_sim, exp


def plot_single_fit(gene, data, fit_res, ax, show_rmse = True, show_title = False, capsize = 10):
    # get simulation values for the gene
    df1, df2 = prep_fit(gene, data, fit_res)

    # plot data and fim values
    palette = sns.color_palette()
    palette_reordered = [palette[1], palette[0], palette[2]]
    sns.scatterplot(data = df2, x = "time", y = "avg_norm_rtm2", ax = ax, color = "k")
    # add some error bars if available
    if not np.isnan(df2.SD).any():
        ax.errorbar(df2.time, df2.avg_norm_rtm2, yerr = df2.SD, ecolor = "k", fmt = "none", capsize = capsize)


    sns.lineplot(data = df1, x = "time", y = "fit1", ax = ax, color = palette[1])
    sns.lineplot(data = df1, x = "time", y = "fit2", ax = ax, color = palette[0])
    sns.lineplot(data = df1, x = "time", y = "fit3", ax = ax, color = palette[2])

    ax.set_ylabel("expr. norm.")
    ax.set_xlabel("time (h)")

    if show_rmse:
        rmse = fit_res.loc[gene]["rmse"]
        names = ["exp", "del", "fat"]
        rmse = [str(np.round(val,2))+n for val, n in zip(rmse, names)]
        title = gene + " " + " ".join(rmse)
        ax.set_title(title)

    if show_title:
        assert show_rmse == False
        ax.set_title(gene)