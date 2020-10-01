#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:54:25 2020

@author: burt
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import gammainc


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


df_kinetic = pd.read_csv("../output/kinetics_tidy.csv")

df_kinetic = df_kinetic[df_kinetic.kinetic == True]

df_err = df_kinetic[["infection", "cell_type", "gene_name", "time", "sd"]]
df_err = df_err.drop_duplicates()
df_val = df_kinetic[["infection", "cell_type", "gene_name", "time", "expr_norm"]]
df_val = df_val.drop_duplicates()


# need to make wide again for fit
df_list = [df_val, df_err]
df_list = [pd.pivot_table(data = df, 
                          index = ["infection", "cell_type", "gene_name"], 
                          columns = "time") for df in df_list]
df_list = [df.reset_index() for df in df_list]

# extract data and errors
n_timepoints = len(df_kinetic.time.drop_duplicates())

vals = df_list[0].iloc[:,-n_timepoints:].values
errs = df_list[1].iloc[:,-n_timepoints:].values

# xdata for gamma dist fit (time series from crawford et al)
x = np.array([0, 6.0, 8.0, 15.0, 30.0])


gamma_list = [gamma_cdf1, gamma_cdf2, gamma_cdf3]
err_list = [[], [], []]

for fun, err in zip(gamma_list, err_list):
    # fit gamma dist for each gene
    for i in range(len(vals)):
        y = vals[i,:]
        
        assert np.any(y==1.0)
        max_idx = np.where(y==1.0)[0][0]
        y = y[:max_idx]
        # kick out nans in y and later in sigma
        #y = y[~np.isnan(y)]
        
        #alpha_fit = np.nan
        chisq = np.nan
        # only do fit procedure if at least 3 vals in y to fit
        if len(y) > 2:
            xdata = x[:len(y)]
            sigma = errs[i,:]
            #sigma = sigma[~np.isnan(sigma)]
            sigma = sigma[:max_idx]   
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
df_final["alpha10"] = err_list[2]


# nas were generated during fit if not enough data was available
df_final = df_final.dropna()
df_final = df_final[["infection", "cell_type", "gene_name", "alpha1", "alpha2", "alpha10"]]
df_final = df_final.reset_index(drop=True)
df_final.to_csv("../output/kinetics_gamma_fit.csv", index = False)
