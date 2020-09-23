#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:48:39 2020

@author: burt
"""
import pandas as pd
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks")
df = pd.read_csv("../output/avg_expression_norm.csv")

# only keep expression values
df_vals = df.iloc[:,-5:]
arr = df_vals.values

# normalization doncition: day 0
norm = arr[:,0]

# define thresold crit how much vals should be greater than norm condition
crit = 1.1
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

df4 = df2.copy()
df4.iloc[:,-5:] = df3

df4 = pd.melt(df4, id_vars = ["infection", "cell_type", "Gene symbol"],
              var_name = "time", value_name = "expr")

df4["time"] = pd.to_numeric(df4["time"])

g = sns.relplot(data = df4, x = "time", y = "expr", 
                hue = "Gene symbol", 
                col = "infection",
                row = "cell_type",
                kind = "line",
                facet_kws = {"margin_titles":True})

g.set_titles(row_template="{row_name}", col_template = "{col_name}")

