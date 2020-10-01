# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:35:29 2020

@author: Philipp
"""


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


sns.set(context = "poster", style = "ticks")
data = pd.read_csv("../output/crawford_annotation.csv")

genes = len(np.unique(data.gene_name))
print("%s uniquely annotated genes found" % genes)

# compute mean expr
data["mean_expr"] = data.groupby(['infection', 'cell_type', 'gene_name', 'time'])["value2"].transform("mean")

# get median expr across summed time points and filter
vals = data.value2.values

fig, ax = plt.subplots()
ax.hist(vals)
ax.set_xlabel("expression")
ax.set_ylabel("n genes")
plt.show()

# keep only genes that are greater than median for one condition
perc = 30
crit_percentile = np.percentile(vals, perc)
data2 = data[data.value2 > crit_percentile]
genes = np.unique(data2.gene_name)
n_genes = len(genes)

data = data[data["gene_name"].isin(genes)]
print("%s genes after expr threshold filter" % n_genes)


# get mean expr at time 0
d0 = data[data.time == 0]
d0 = d0[["infection", "cell_type", "gene_name", "mean_expr"]]
d0 = d0.rename(columns = {"mean_expr" : "mean_expr_d0"})
d0 = d0.drop_duplicates()
data2 = pd.merge(data, d0, on = ["infection", "cell_type", "gene_name"])

# compute fold changes
data2["fc_expr"] = data2.value2 / data2.mean_expr_d0
data2["fc_mean"] = data2.mean_expr / data2.mean_expr_d0

data2["log2FC"] = np.log2(data2.fc_expr)
data2["log2FC_mean"] = np.log2(data2.fc_mean)

# apply fold change criterion on averages
data3 = data2[["time", "infection", "cell_type", "gene_name", "log2FC_mean"]]
data3 = data3.drop_duplicates()
data3["log2FC_mean"] = np.abs(data3.log2FC_mean)

# get number of time points
n_timepoints = len(data3.time.drop_duplicates())


# convert to wide format
data3 = pd.pivot_table(data = data3, 
                       index = ["infection", "cell_type", "gene_name"], 
                       columns = "time")
data3 = data3.reset_index()

# get actual fc values
df_vals = data3.iloc[:,-n_timepoints:]
arr = df_vals.values

# apply fc criterion on log2fc absolute values
fc_crit = 0.25

arr2 = arr > fc_crit

arr3 = arr2[:,1:] & arr2[:,:-1]
# check if there is at least one pair of adjacent true vals
arr3 = np.sum(arr3, axis = 1)
adj_pairs = 1
arr3 = arr3 > adj_pairs

# subset df for genes that change
data4 = data3[arr3]
data4 = data4.reset_index(drop = True)

genes_kinetic = np.unique(data4.gene_name)
n_genes = len(genes_kinetic)
print("%s genes after fc crit 2subseq. days filter" % n_genes)


df_kinetic = data2[data2["gene_name"].isin(genes_kinetic)]


# compute sd
df_kinetic["log2fc_abs"] = np.abs(df_kinetic.log2FC)
df_kinetic["log2fc_mean_abs"] = np.abs(df_kinetic.log2FC_mean)
df_kinetic["max"] = df_kinetic.groupby(['infection', 'cell_type', 'gene_name'])["log2fc_mean_abs"].transform("max")
df_kinetic["sd"] = df_kinetic.groupby(['infection', 'cell_type', 'gene_name', 'time'])["log2fc_abs"].transform("std")
df_kinetic["expr_norm"] = df_kinetic["log2fc_mean_abs"] / df_kinetic["max"]


# add info if gene is kinetic in specific condition
df5 = data4.melt(id_vars = ["infection", "cell_type", "gene_name"])
df5 = df5[["infection", "cell_type", "gene_name", "time"]]
df5["kinetic"] = True
df6 = pd.merge(df_kinetic, df5, how = "left", on = ["infection", "cell_type", "gene_name", "time"])
df6.kinetic = df6.kinetic.fillna(False)

df6.to_csv("../output/kinetics_tidy.csv", index = False)


# plot fold change distr.
vals = data2.log2FC_mean.values
fig,ax = plt.subplots()
ax.hist(vals)
ax.set_ylim(0, 100000)
ax.set_xlim(-1,1)
ax.set_xlabel("expr fold-change")
ax.set_ylabel("n genes")