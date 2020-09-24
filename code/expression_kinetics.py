#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:48:39 2020

@author: burt
"""
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.special import gammainc

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)


sns.set(context = "poster", style = "ticks")
data = pd.read_csv("../output/avg_expression_norm.csv")

# =============================================================================
# split and to wide format
# =============================================================================
df_raw_val = data[["infection", "cell_type", "Gene symbol", "time", "avg_norm"]]
df_raw_err = data[["infection", "cell_type", "Gene symbol", "time", "err"]]
df_rtm_val = data[["infection", "cell_type", "Gene symbol", "time", "avg_norm_rtm2"]]
df_rtm_err = data[["infection", "cell_type", "Gene symbol", "time", "err_rtm"]]

df_list = [df_raw_val, df_raw_err, df_rtm_val, df_rtm_err]
df_list = [pd.pivot_table(data = df, 
                          index = ["infection", "cell_type", "Gene symbol"], 
                          columns = "time") for df in df_list]
df_list = [df.reset_index() for df in df_list]

# =============================================================================
# raw data find genes that change at all
# =============================================================================
df = df_list[0]
df_vals = df.iloc[:,-5:]
arr = df_vals.values

# normalization doncition: day 0
norm = arr[:,0]

# define thresold crit how much vals should be greater than norm condition
crit = 1.3
norm = norm*crit

# add extra dimension to normalization conditino(d0) because then I can broadcast to whole arr
norm = norm[:, None]

# compare array with itself shifted by 1 column, and get bitwise AND&, this gives true for two adjacent trues
arr2 = arr > norm
arr3 = arr2[:,1:] & arr2[:,:-1]
# check if there is at least one pair of adjacent true vals
arr3 = np.sum(arr3, axis = 1)
adj_pairs = 1
arr3 = arr3 > adj_pairs

# subset df for genes that change
df2 = df[arr3]
df2 = df2.reset_index(drop = True)

# =============================================================================
# raw data add Nans after plateau phase
# =============================================================================
# kick out values that decrease after plateau
df3 = df2.iloc[:,-5:].values
# for each row check if row is monotonic, if not find first deviating element
for i in range(len(df3)):
    row = df3[i,:]
    # check if row is monotonic
    desc_thres = 0.1
    idx = np.where((row[1:]-row[:-1]) < -desc_thres)
    idx = idx[0]
    # check if element is empty if not assign everything to nan
    if idx.size != 0:
        idx = idx[0]
        df3[i, (idx+1):] = np.nan

# update nans in original df
df4 = df2.copy()
df4.iloc[:,-5:] = df3

# =============================================================================
# visualize
# =============================================================================
df4 = pd.melt(df4, id_vars = ["infection", "cell_type", "Gene symbol"])
df4 = df4[["infection", "cell_type", "Gene symbol", "time", "value"]]

df4["time"] = pd.to_numeric(df4["time"])

g = sns.relplot(data = df4, x = "time", y = "value", 
                hue = "Gene symbol", 
                col = "infection",
                row = "cell_type",
                kind = "line",
                facet_kws = {"margin_titles":True})

g.set_titles(row_template="{row_name}", col_template = "{col_name}")


# =============================================================================
# use response time normed data subset from meaningful raw data
# =============================================================================
# filter for values where sth happens at all and only for plateau phase (all in df4)
df5 = df4.dropna()
df5 = df5.drop(["value"], axis = 1)

df_rtm_val2 = pd.merge(df5, df_rtm_val, how = "left")
df_rtm_err2 = pd.merge(df5, df_rtm_err, how = "left")


g = sns.relplot(data = df_rtm_val2, x = "time", y = "avg_norm_rtm2", 
                hue = "Gene symbol", 
                col = "infection",
                row = "cell_type",
                kind = "line",
                facet_kws = {"margin_titles":True})

g.set_titles(row_template="{row_name}", col_template = "{col_name}")


# need to make wide again for fit
df_list = [df_rtm_val2, df_rtm_err2]
df_list = [pd.pivot_table(data = df, 
                          index = ["infection", "cell_type", "Gene symbol"], 
                          columns = "time") for df in df_list]
df_list = [df.reset_index() for df in df_list]

# extract data and errors
vals = df_list[0].iloc[:,-5:].values
errs = df_list[1].iloc[:,-5:].values

# xdata for gamma dist fit (time series from crawford et al)
x = np.array([0, 6.0, 8.0, 15.0, 30.0])
alpha_list = []
cov_list = []

# fit gamma dist for each gene
for i in range(len(vals)):
    y = vals[i,:]
    # kick out nans in y and later in sigma
    y = y[~np.isnan(y)]
    
    alpha_fit = np.nan
    alpha_err = np.nan
    # only do fit procedure if at least 3 vals in y to fit
    if len(y) > 3:
        xdata = x[:len(y)]
        sigma = errs[i,:]
        sigma = sigma[~np.isnan(sigma)]
        fit_val, fit_err = curve_fit(f = gamma_cdf, xdata= xdata, ydata = y, sigma = sigma,
                                     absolute_sigma = True)

        # according to docs this is the error of covariance of fitted params
        # I am only interested in the identifiability of alpha, so take only first param        
        alpha_err =  np.sqrt(np.diag(fit_err))[0]
        alpha_fit = fit_val[0]
        
    # add alpha fit and err
    alpha_list.append(alpha_fit)
    cov_list.append(alpha_err)    


#bounds=([1.0,0], np.inf)

df_final = df_list[0]
df_final["alpha"] = alpha_list
df_final["fit_err"] = cov_list

