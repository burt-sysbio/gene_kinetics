#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:24:52 2020

@author: burt
new script for output generated from R should work on all kinetic data together
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


def prep_data(df):

    df_err = df[["cell_type", "gene_name", "time", "SD"]]
    df_err = df_err.drop_duplicates()

    df_val = df[["cell_type", "gene_name", "time", "avg_norm_rtm2"]]
    df_val = df_val.drop_duplicates()
    
    
    # need to make wide again for fit
    df_list = [df_val, df_err]
    df_list = [pd.pivot_table(data = df, 
                              index = ["cell_type", "gene_name"], 
                              columns = "time") for df in df_list]
    df_list = [df.reset_index() for df in df_list]
    
    # extract data and errors
    timepoints = df.time.drop_duplicates().values
    n_timepoints = len(timepoints)
    
    vals = df_list[0].iloc[:,-n_timepoints:].values
    errs = df_list[1].iloc[:,-n_timepoints:].values    

    return(df_list, vals, errs, timepoints)



def fit_kinetic(df, rep_available, gamma_list = [gamma_cdf1, gamma_cdf2, gamma_cdf3]):

    df_list, vals, errs, x = prep_data(df)

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
                
                if rep_available == True:
                    sigma = errs[i,:]
                    #sigma = sigma[~np.isnan(sigma)]
                    sigma = sigma[:max_idx]   
                    absolute_sigma = True
                else:
                    sigma = None
                    absolute_sigma = False
                try:
                    # run fit catch runtime exception for bad fit
                    # restrain alpha within 1 and 100
                    fit_val, fit_err = curve_fit(f = fun, 
                                                 xdata= xdata, 
                                                 ydata = y, 
                                                 sigma = sigma,
                                                 absolute_sigma = absolute_sigma)
            
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
    
    return(df_list, err_list)


def run_fit(df, data, cell, inf, rep_available, gamma_list = [gamma_cdf1, gamma_cdf2, gamma_cdf3]):
    
    df_list, err_list = fit_kinetic(df, rep_available = rep_available, gamma_list = gamma_list)
    df_final = df_list[0]
    df_final["alpha1"] = err_list[0]
    df_final["alpha2"] = err_list[1]
    df_final["alpha10"] = err_list[2]
    
    
    # nas were generated during fit if not enough data was available
    df_final = df_final.dropna()
    df_final = df_final[["cell_type", "gene_name", "alpha1", "alpha2", "alpha10"]]
    df_final = df_final.reset_index(drop=True)
    
    filename = "../output/gamma_fit_"+data+"_"+cell+"_"+inf+".csv"
    df_final.to_csv(filename, index = False)
#bounds=([1.0,0], np.inf)
# xdata for gamma dist fit (time series from crawford et al)

df = pd.read_csv("../output/data_rtm_CD4_arm.csv")
data = "crawford"
cell = "CD4"
inf = "arm"
rep_available = True

run_fit(df, data, cell, inf, rep_available)






