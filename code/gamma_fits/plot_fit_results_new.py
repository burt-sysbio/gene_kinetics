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

def tidy_df(df, timepoints, alpha, n_genes = 3):
    # keep cols where desired alpha is best fit and of those only good fits   
    df = df[df.best_alpha == alpha]
    thres = 0.1
    df = df[df.err_alpha < thres]
    
    # sort by bad fits for other alpha
    sort_by = "alpha10" if alpha == 1 else "alpha1"
    df = df.sort_values(sort_by, ascending = False)
    df = df.iloc[:n_genes,:]
    
    print(df[["gene_name", "alpha1", "alpha10"]])    
    # convert to long format for timepoints

    id_vars = ["gene_name", "beta_1", "beta_10"]
    df = df[id_vars + timepoints]
    df = df.melt(id_vars = id_vars, var_name = "time", value_name = "exp_data")
    df = df.dropna()
    df["time"] = pd.to_numeric(df["time"])
    
    # convert to long format for beta fit
    #df = df.melt(id_vars = ["gene_name", "time", "exp_data"], 
    #                       var_name = "beta", value_name = "beta_fit")

    # add alpha to each corresonding beta 
    #df["alpha"] = 1
    #df["alpha"][df.beta == "beta_10"] = 10    
    # now beta column is not needed 
    #df = df.drop(["beta"], axis = 1)
    
    # for each gene at each time compute gamma cdf for alpha and beta fit 
    #df["model_fit"] = df.apply(get_ydata, axis = 1)
    return(df)
    

def get_ydata(row):
    beta = row["beta_fit"]
    alpha = row["alpha"]
    t = row["time"]
    model_fit = gamma_cdf(t, alpha, beta)
    return(model_fit)


# focus on specific cell type and infection
my_df = df_craw[(df_craw.infection == "arm") & (df_craw.cell_type == "CD4")]
timepoints = [0,6,8,15,30]
timepoints = [str(val) for val in timepoints]

df_rtm = tidy_df(df_craw, timepoints, alpha = 10)
#g = sns.relplot(data = df_rtm, x = "time", y = "exp_data", col = "gene_name")

g = sns.relplot(data = df_rtm, x = "time", y = "exp_data", col = "gene_name",
                kind = "scatter")

# get unique genes
genes = df_rtm.gene_name.drop_duplicates()
genes = genes.values
for (ax, gene) in zip(g.axes.flat, genes):
    p = df_rtm.loc[df_rtm.gene_name == gene]
    print(p)
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
    