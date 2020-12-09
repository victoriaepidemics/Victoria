#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 09:27:13 2020

@author: jac

SEID model analysis of the Bubonic plague outbreak in Eyam, England

Data from the second black plague outbreak in the village of Eyam,
England, from 19 June to 1 November 1666 (Stojkova and D. A. Campbell 2017).
"""

import os
from datetime import date
from victoriaepi.seid import SEID
import matplotlib.pyplot as plt

workdir = "./"

# Residence rates 1/day, names and Erlang series
R_rates = \
        {'E':[1/5       , r'\sigma_1',  4],
         'I':[1/4      , r'\sigma_2',  3]}

trim=-30
eyam = SEID( Region="Eyam, England, Bubonic plague outbreak 1,666",\
                N=261, data_fnam="Eyam.csv", out_fnam="Eyam",\
                init_index=0, init=date( 1666, 6, 1), R_rates=R_rates, trim=trim, workdir= workdir)

os.makedirs(workdir+"output", exist_ok=True)

T=2000
pred=-trim+7*3
if T>0:
    eyam.RunMCMC( T=T, burnin=1000, pred=pred, plot_fit=False)

if os.environ["DISPLAY"]:
    eyam.PlotEvolution(pred=pred, cumm=True, right_axis=False)
    plt.show()
else:
    print ("No display found")