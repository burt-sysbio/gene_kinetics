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
import matplotlib.pyplot as plt

def gamma_cdf(t, alpha, beta):
    dummy = t*beta
    return gammainc(alpha,dummy)


def gamma_cdf1(t, beta):
    dummy = t*beta
    return gammainc(1.0,dummy)


def gamma_cdf2(t, beta):
    dummy = t*beta
    return gammainc(2,dummy)


def gamma_cdf3(t, beta):
    dummy = t*beta
    return gammainc(10.0,dummy)


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


a = df_raw_val["Gene symbol"].values
print(a, len(np.unique(a)))
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
    desc_thres = 0.05
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
                facet_kws = {"margin_titles":True},
                legend = False)

g.set(ylabel = "expr. norm.")
g.set_titles(row_template="{row_name}")


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
                facet_kws = {"margin_titles":True},
                legend = False)

g.set_titles(row_template="{row_name}")
g.set(ylabel = "expr. norm.")

#g.savefig("../figures/gene_kinetics_norm.pdf")

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


gamma_list = [gamma_cdf1, gamma_cdf2, gamma_cdf3]
err_list = [[], [], []]

for fun, err in zip(gamma_list, err_list):
    # fit gamma dist for each gene
    for i in range(len(vals)):
        y = vals[i,:]
        # kick out nans in y and later in sigma
        y = y[~np.isnan(y)]
        
        #alpha_fit = np.nan
        chisq = np.nan
        # only do fit procedure if at least 3 vals in y to fit
        if len(y) > 2:
            xdata = x[:len(y)]
            sigma = errs[i,:]
            sigma = sigma[~np.isnan(sigma)]
            
            try:
                # run fit catch runtime exception for bad fit
                # restrain alpha within 1 and 100
                fit_val, fit_err = curve_fit(f = fun, 
                                             xdata= xdata, 
                                             ydata = y, 
                                             sigma = sigma,
                                             absolute_sigma = True)
        
                # according to docs this is the error of covariance of fitted params
                alpha_err =  np.sqrt(np.diag(fit_err))[0]
                
                # compute chi sqr
                nexp = fun(xdata, *fit_val)
                # get residuals
                r = y - nexp
                # only reassign chisq for low error in identified parameter              
                if alpha_err < 1.0:
                    chisq = np.sum((r/sigma)**2)
            except RuntimeError:
                print("max calls reached no good fit")
        # add alpha fit and err
        #alpha_list.append(alpha_fit)
        #cov_list.append(alpha_err)    
        err.append(chisq)

#bounds=([1.0,0], np.inf)

df_final = df_list[0]
df_final["alpha1"] = err_list[0]
df_final["alpha2"] = err_list[1]
df_final["alpha3"] = err_list[2]

df_fit_res = pd.DataFrame({"alpha1": err_list[0], 
                           "alpha2" : err_list[1], 
                           "alpha10" : err_list[2]})


df_fit_res = df_fit_res.dropna()

n_alpha1 = sum((df_fit_res.alpha1 < df_fit_res.alpha2) & (df_fit_res.alpha1 < df_fit_res.alpha10))
n_alpha2 = sum((df_fit_res.alpha2 < df_fit_res.alpha1) & (df_fit_res.alpha2 < df_fit_res.alpha10))
n_alpha10 = sum((df_fit_res.alpha10 < df_fit_res.alpha1) & (df_fit_res.alpha10 < df_fit_res.alpha2))

x = [r"$\alpha=1$", r"$\alpha=2$", r"$\alpha=10$"]
y = [n_alpha1, n_alpha2, n_alpha10]

df = pd.DataFrame({"best_fit":x, "n_genes":y})

fig, ax = plt.subplots()
ax = sns.barplot(data = df, x = "best_fit", y = "n_genes", color = "tab:blue")
plt.tight_layout()
#fig.savefig("../figures/gene_kinetics_fit_results.pdf")

#df_fit_res["crit"] = crit

#n_hits = sum(df_fit_res.crit)
#n_total = len(df_fit_res.crit)
#print(n_hits, n_total)
