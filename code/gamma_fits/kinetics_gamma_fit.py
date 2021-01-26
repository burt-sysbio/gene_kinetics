#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:24:52 2020

@author: burt
fit gamma dist for output generated from R should work on all kinetic data together
"""

import pandas as pd
from utils import run_f_test

study = "powrie"
groups = ["innate"]
activations = ["colitis"]
month = "jan2021"

readdir = "../../data/data_rtm/"
savedir = "../../output/gamma_fits/" + month + "/"

for g in groups:
    for act in activations:

        file = readdir+study + "_rtm_" + g + "_" + act + ".csv"
        df = pd.read_csv(file)

        # if rna seq data, kick out nans generated through rlog trafo
        if study == "proserpio":
            genes_na = df.gene[df.SD.isna()]
            df = df[~df.gene.isin(genes_na)]

        # kick out NAs, at least in properio there are some nas for some reason
        # I need to double check this
        a,b,c = run_f_test(df)

        dummy = study + "_" + g + "_" + act + ".csv"
        file_a = savedir + "fit_expo_" + dummy
        file_b = savedir + "fit_gamma_" + dummy
        file_c = savedir + "ftest_" + dummy

        a.to_csv(file_a)
        b.to_csv(file_b)
        c.to_csv(file_c)


