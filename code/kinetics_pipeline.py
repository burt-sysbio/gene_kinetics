# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:35:29 2020

@author: Philipp
"""


import pandas as pd
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.special import gammainc
import matplotlib.pyplot as plt

sns.set(context = "poster", style = "ticks")
data = pd.read_csv("../output/crawford_annotation.csv")


# compute mean expr
data["mean_expr"] = data.groupby(['infection', 'cell_type', 'gene_name', 'time'])["value2"].transform("mean")

# get median expr across summed time points and filter
data2 = data[["infection", "cell_type", "gene_name", "time", "mean_expr"]]
data2 = data2.drop_duplicates()
data2 = data2.reset_index(drop = True)

sum_times = data2.groupby(["gene_name", "infection", "cell_type"])["mean_expr"].agg("sum")
median = np.median(sum_times.values)

# keep only genes that are greater than median for one condition
sum_times = sum_times[sum_times > median]

genes = sum_times.reset_index()
genes = genes["gene_name"]
genes = genes.drop_duplicates()

fig, ax = plt.subplots()
ax.hist(sum_times.values)
ax.set_xlabel("median expression")
ax.set_ylabel("n genes")
plt.tight_layout()

# subset original data by genes whose summed score is greater than median
data = data[data["gene_name"].isin(genes.values)]

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


df_vals = data3.iloc[:,-n_timepoints:]
arr = df_vals.values

# apply fc criterion on log2fc absolute values
fc_crit = 0.5

arr2 = arr > fc_crit

arr3 = arr2[:,1:] & arr2[:,:-1]
# check if there is at least one pair of adjacent true vals
arr3 = np.sum(arr3, axis = 1)
adj_pairs = 1
arr3 = arr3 > adj_pairs

# subset df for genes that change
data3 = data3[arr3]
data3 = data3.reset_index(drop = True)
