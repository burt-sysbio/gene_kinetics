#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:24:52 2020

@author: burt
fit gamma dist for output generated from R should work on all kinetic data together
"""

import pandas as pd
from utils import run_f_test
df = pd.read_csv("../../output/rtm_data/peine_rtm_Th0_invitro.csv")
fit_res = run_f_test(df)
fit_res.to_csv("../../output/gamma_fits/dec2020/gamma_fit_peine_Th0.csv")

