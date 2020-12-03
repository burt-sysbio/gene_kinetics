#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:24:52 2020

@author: burt
fit gamma dist for output generated from R should work on all kinetic data together
"""

import pandas as pd
import numpy as np

from utils import run_fit
print("hi")
bounds=([1.0,0], np.inf)
#xdata for gamma dist fit (time series from crawford et al)
#
data = ["peine", "crawford"]
for d in data:
    if d == "peine":
        rep_available = False
        inf = "invitro"
        cell = ["Th0", "Th1", "Th2", "ThMix"]
        for c in cell:
            file = "output/rtm_data/"+d+"_rtm_"+c+"_"+inf+".csv"

            df = pd.read_csv(file)
            df2 = run_fit(df, data = d, cell = c, inf = inf, rep_available = rep_available)

    else:
        rep_available = True
        inf = ["arm", "cl13"]
        cell = ["CD4", "CD8"]

        for c in cell:
            for i in inf:
                file = "output/rtm_data/"+d+"_rtm_"+c+"_"+i+".csv"

                df = pd.read_csv(file)
                df2 = run_fit(df, data = d, cell = c, inf = i, rep_available = rep_available)