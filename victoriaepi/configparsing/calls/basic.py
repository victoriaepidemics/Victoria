#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic model
"""
import sys
import os
from datetime import date
import matplotlib.pyplot as plt
from $moduleName import $className

# Residence rates 1/day, names and Erlang series
R_rates=$R_rates
trim=$trim
workdir="$workdir"
$instanceName = $className(Region="$Region",
    N=$N, data_fnam="$data_fnam", out_fnam="$out_fnam",
    init_index=$init_index, init=$init_date, R_rates=R_rates, trim=trim, workdir=workdir)
T=$T
pred = $pred#-trim+7*3


os.makedirs(workdir, exist_ok=True)
os.makedirs(workdir+"output", exist_ok=True)
os.makedirs(workdir+"figs", exist_ok=True)
os.makedirs(workdir+"data", exist_ok=True)
os.makedirs(workdir+"csv", exist_ok=True)

if T>0:
    $instanceName.RunMCMC( T=T, burnin=$burnin, pred=pred, plot_fit=$plot_fit)


if os.environ["DISPLAY"]:
    $instanceName.PlotEvolution(pred=pred, cumm=True, right_axis=False)
    $instanceName.PlotEvolution(pred=pred, cumm=False, right_axis=False)
    plt.show()
else:
    print ("No display found")

