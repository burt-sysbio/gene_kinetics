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

alpha = 0.8
s = 60
g = sns.jointplot(data = df_all, x = "alpha1", y = "alpha_2_10", height = 10, 
                  hue = "study", kind = "scatter", alpha = alpha, s = s)
g.ax_joint.set(yscale="log", xscale = "log")
g.ax_joint.plot([0,100], [0,100], '--k') 


# zoom in
lower = 1e-4
upper = 1e2
g = sns.jointplot(data = df_all, x = "alpha1", y = "alpha_2_10", height = 10, 
                  hue = "study", kind = "scatter", alpha = alpha, s = s,
                  xlim = (lower,upper), ylim = (lower, upper))
g.ax_joint.set(yscale="log", xscale = "log")
g.ax_joint.plot([0,100], [0,100], '--k') 


# kick out genes with no annotation
df3 = df_all.loc[df_all.category != "other"]
lower = 1e-5
upper = 1e1
g = sns.jointplot(data = df3, x = "alpha1", y = "alpha_2_10", height = 10, 
                  hue = "category", kind = "scatter",
                  xlim = (lower,upper), ylim = (lower, upper))
g.ax_joint.set(yscale="log", xscale = "log")
g.ax_joint.plot([0,100], [0,100], '--k') 
#g.set(yscale="log", xscale = "log", )

# kick out genes with no annotation
df4 = df_craw.loc[df_craw.category != "other"]
lower = 1e-5
upper = 1e1
g = sns.jointplot(data = df3, x = "alpha1", y = "alpha_2_10", height = 10, 
                  hue = "category", kind = "scatter",
                  xlim = (lower,upper), ylim = (lower, upper))
g.ax_joint.set(yscale="log", xscale = "log")
g.ax_joint.plot([0,100], [0,100], '--k') 
#g.set(yscale="log", xscale = "log", )



# kick out genes with no annotation
df5 = df_peine.loc[df_peine.category != "other"]
lower = 1e-5
upper = 1e1
g = sns.jointplot(data = df3, x = "alpha1", y = "alpha_2_10", height = 10, 
                  hue = "category", kind = "scatter",
                  xlim = (lower,upper), ylim = (lower, upper))
g.ax_joint.set(yscale="log", xscale = "log")
g.ax_joint.plot([0,100], [0,100], '--k') 
#g.set(yscale="log", xscale = "log", )




# compare cl13 and arm
alpha = 0.8
s = 60
g = sns.relplot(data = df_craw, x = "alpha1", y = "alpha_2_10", 
                col = "infection", hue = "cell_type", kind = "scatter", alpha = alpha, s = s)
g.set(yscale="log", xscale = "log", xlim = (1e-6, 1e2), ylim = (1e-6, 1e2))



# old stuff

#df4 = df3[["cell_type", "gene_name", "err_alpha", "best_alpha", "infection", "study", "category"]]

#g = sns.catplot(data = df3, x = "best_alpha", y = "err_alpha", hue = "study", kind = "violin")

# get number of alphas
#best_alpha = df4.groupby(["best_alpha", "cell_type", "study", "category"])["gene_name"].count()
#best_alpha = best_alpha.to_frame()

#best_alpha = best_alpha.reset_index()

#g = sns.catplot(data = best_alpha, x = "best_alpha", hue = "cell_type", y = "gene_name", kind = "bar", col = "category")


# get genes for which
df_exp = df_craw[df_craw.best_alpha == 10]
df_exp = df_exp.sort_values("diff_err", ascending = False)
df_exp = df_exp.iloc[1:5,:]
df_exp2 = df_exp.drop(["alpha1", "alpha2", "alpha10", "best_alpha", "err_alpha", "alpha_2_10", "diff_err"], axis = 1)
id_vars = ["cell_type", "infection", "study", "category", "gene_name", "beta_1", "beta_2", "beta_10"]
df_exp3 = df_exp2.melt(id_vars = id_vars, var_name = "time", value_name = "exp_data")
df_exp3 = df_exp3.dropna()
df_exp3 = df_exp3.melt(id_vars = ["cell_type", "infection", "study", "category", "gene_name", "time", "exp_data"], 
                       var_name = "beta", value_name = "beta_fit")

df_exp3["time"] = pd.to_numeric(df_exp3["time"])
df_exp3["alpha"] = 1
df_exp3["alpha"][df_exp3.beta == "beta_2"] = 2
df_exp3["alpha"][df_exp3.beta == "beta_10"] = 10

def get_ydata(row):
    beta = row["beta_fit"]
    alpha = row["alpha"]
    t = row["time"]
    print(type(t), alpha, beta)
    model_fit = gamma_cdf(t, alpha, beta)
    return(model_fit)

df_exp3["model_fit"] = df_exp3.apply(get_ydata, axis = 1)