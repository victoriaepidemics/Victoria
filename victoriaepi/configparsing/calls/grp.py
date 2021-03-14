#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Marcos Capistran, Antonio capella y Andres Christen,
Marzo-Octubre 2020.

"""

import sys
from datetime import date, timedelta
import numpy as np
from numpy import array, ones
from scipy.stats import beta
from matplotlib import use
from matplotlib.pyplot import subplots, close, rcParams
from pandas import datetime, read_csv
import os
from victoriaepi.plotfrozen import PlotFrozenDist


from $moduleName import $className

# Number of groups
ngrp = $ngrp  # 4

# From Ferguson:
#F_g = array([0.004     , 0.0348    , 0.12333333, 0.2396    ])
#F_h = array([0.05      , 0.0552    , 0.17266667, 0.5112    ])

# Exit probs.,   f same x grp,   g,    h,   i same x grp
exit_probs = $exit_probs #[array([0.05]*ngrp), F_g,  F_h, array([0.5]*ngrp)]
"""
[array([0.05, 0.05, 0.05, 0.05]),
 array([0.004     , 0.0348    , 0.12333333, 0.2396    ]),
 array([0.05      , 0.0552    , 0.17266667, 0.5112    ]),
 array([0.5, 0.5, 0.5, 0.5])]
"""

"""
# ngrp X ngrp array with interaction matrix
"""

# Contact matrix:

Int_M = $Int_M # ones((ngrp,ngrp))/ngrp #Uniform



Pobs_I = $Pobs_I
Pobs_D = $Pobs_D

R_rates_V2 = $R_rates


# 'U^1':[1/5 era 1/7 ... 'H^2':[1/6 era 1/8
# R_rates_V2=\
#         { 'E'  :[1/5      , r'\sigma_1',  4],\
#           'I^S':[1/4      , r'\sigma_2',  3],\
#           'H^1':[1/1      , r'\sigma_3',  1],\
#           'U^1':[1/5      , r'\sigma_4',  3],\
#           'U^2':[1/1      , r'\sigma_5',  1],\
#           'I^A':[1/7      , r'\gamma_1',  3],\
#           'I^B':[1/7      , r'\gamma_2',  3],\
#           'H^2':[1/6      , r'\gamma_3',  3],\
#           'H^3':[1/3.3    , r'\gamma_4',  1]}




def AnalyzeZM( clave, T, burnin=5000, init=date(2020, 4, 1), trim=0, pred=100, plot_fit=True, workdir="./../"):

    # Translacion y uniformización de fecha de referencia para graficacion
    init=init
    if (ZMs[clave][3]+ timedelta(days=$daysdelay)) <= init:
        # Including the 4-day shift (symptoms onset)
        init_index=(init - (ZMs[clave][3] + timedelta(days=$daysdelay))).days
    else:
        init=ZMs[clave][3] + timedelta(days=$daysdelay)
        init_index=0

    if trim == 0:
        out_fnam = "$out_fnam"  # clave+"_grp"
    else:
        out_fnam = "$out_fnam" + "_trim%d" % (trim,)  # clave+"_trim%d" % (trim,)
    data_fnam= "$data_fnam"  # clave+'.csv'

    intervention_day = [(d-ZMs[clave][3]).days-4 for d in ZMs[clave][4:-(int(ZMs[clave][1])+1)]]

    exit_probs_copy = exit_probs.copy()
    if int(ZMs[clave][1]) > 0:
        relax_day = [(d - ZMs[clave][3]).days-4 for d in ZMs[clave][(-(int(ZMs[clave][1])+1)):-1]]
        exit_probs_copy[0] = array([0.4]*ngrp) #=f
        print("ama-2")
    else:
        relax_day = []

    N = ZMs[clave][2]
    age_prop = array($age_prop)
    # if clave == "9-01":
    #     age_prop = array([0.45213552, 0.36863423, 0.11800457, 0.06122568]) # Piramide para zm_vmx 2010
    # else:
    #     age_prop = array([0.48573906, 0.34280597, 0.10773953, 0.06371544]) # Piramide general para zm's del pais 2010

    zm =  $className( Region=ZMs[clave][0], data_fnam=data_fnam,\
        N=N, out_fnam=out_fnam, init_index=init_index, init=init,\
        intervention_day=intervention_day, relax_day=relax_day, trim=trim,\
        Pobs_I=Pobs_I, Pobs_D=Pobs_D,\
        R_rates=R_rates_V2, exit_probs=exit_probs_copy, workdir=workdir,\
        ngrp=ngrp, Int_M=Int_M, age_prop=age_prop)
    zm.age_groups = $age_groups  # [0, 25, 50, 65, 100]


    if T > 0:
        zm.RunMCMC(T=T, burnin=burnin, pred=pred, plot_fit=plot_fit)
    return zm




# reading: workdir/data/hospitales_datos/

dateparse = lambda x: datetime.strptime(x, '%Y-%m-%d') #data Reg IRAG


def PlotFigsZMs( zm, pred=99, q=[10,25,50,75,90], blue=True, workdir='./../'):

    close('all')


    try:
        zm_vmx_Hosp_RI = read_csv(workdir + "data/hosp/%s_DinHosp.csv" % (zm.clave,), parse_dates=['fecha'], date_parser=dateparse)
        hosp = True
    except:
        hosp = False

    out_fnam=zm.out_fnam

    fig, ax = subplots( num=1, figsize=(8,6))
    zm.PlotEvolution( pred=4, cumm=False, log=False, ty=0, ax=ax,\
            label='Median', q=q, blue=blue,\
            csv_fnam=workdir + "csv/%sI_short.csv" % (out_fnam,))
    fig.tight_layout()
    fig.savefig("%s%s_I_short.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots( num=2, figsize=(8,6))
    zm.PlotEvolution( pred=int(pred/2), cumm=False, log=False, ty=1, ax=ax,\
            label='Median', q=q, blue=blue, right_axis=False,\
            csv_fnam=workdir + "csv/%sD_short.csv" % (out_fnam,))
    fig.tight_layout()
    fig.savefig("%s%s_D_short.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots( num=3, figsize=(8,6))
    zm.PlotEvolution( pred=pred, cumm=False, log=False, ty=0, ax=ax, q=q, blue=blue, add_MRE=True,\
        label=r'Median', csv_fnam=workdir + "csv/%sI_long.csv" % (out_fnam,))
    fig.tight_layout()
    fig.savefig("%s%s_I_long.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots( num=4, figsize=(8,6))
    zm.PlotEvolution( pred=int(pred/2), cumm=True, log=False, ty=1, ax=ax, q=q, blue=blue,\
        label=r'Median', right_axis=False, csv_fnam=workdir + "csv/%sD_long.csv" % (out_fnam,))
    fig.tight_layout()
    fig.savefig("%s%s_D_long.png" % (workdir + 'figs/',out_fnam))


    colors=[ "black", "dimgrey", "grey", "darkgrey", "lightgrey"]
    fig, ax = subplots(figsize=(8,6))
    for grp in range(1,zm.ngrp):
        zm.PlotStateVar( 'H^2 H^3', grps=grp, pred=pred, every=1, ax=ax,\
            median_only=True, q=q, right_axis=False,\
            label=r'%d - %d' % (zm.age_groups[grp],zm.age_groups[grp+1]), color=colors[zm.ngrp-grp])
    zm.PlotStateVar( 'H^2 H^3', grps='all', pred=pred, every=1, ax=ax, q=q, blue=blue,\
        label=r'Sum', color='red',\
        csv_fnam=workdir + "csv/%sHs.csv" % (out_fnam,))
    if hosp:
        ax.bar(  array(zm_vmx_Hosp_RI['fecha']), array(zm_vmx_Hosp_RI['camasOcupadasHG']), color='red', width=0.5, alpha=0.5)
        ax.plot( zm_vmx_Hosp_RI['fecha'], zm_vmx_Hosp_RI['camasOcupadasHG'], 'ro', markersize=2)
    ax.set_xlabel(' ')
    fig.tight_layout()
    fig.savefig("%s%s_Hs.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots(figsize=(8,6))
    for grp in range(1,zm.ngrp):
        zm.PlotStateVar( 'U^1', grps=grp, pred=pred, every=1, ax=ax,\
            median_only=True, q=q, right_axis=False,\
            label=r'%d - %d' % (zm.age_groups[grp],zm.age_groups[grp+1]), color=colors[zm.ngrp-grp])
    zm.PlotStateVar( 'U^1', grps='all', pred=pred, every=1, ax=ax, q=q, blue=blue,\
        label=r'Sum', color='red',right_axis=False,\
        csv_fnam=workdir + "csv/%sU1.csv" % (out_fnam,))
        #ax.set_ylim((0,2100))
    if hosp:
        ax.bar( array(zm_vmx_Hosp_RI['fecha']), array(zm_vmx_Hosp_RI['camasOcupadasUCI_NoUCIVent']), color='red', width=0.5, alpha=0.5)
        ax.plot( zm_vmx_Hosp_RI['fecha'], zm_vmx_Hosp_RI['camasOcupadasUCI_NoUCIVent'], 'ro', markersize=1)
    ax.set_title("%s ICU hospital demand:\nNumber of beds required per day" % (zm.Region,))
    ax.set_xlabel(' ')
    fig.tight_layout()
    fig.savefig("%s%s_U1.png" % (workdir + 'figs/',out_fnam))

    # Estas dos figuras son nuevas
    if zm.relax_day != []:
        fig, ax = subplots( num=7, figsize=(8,6))
        zm.PlotStateVar( 'R', pred=pred, ax=ax, q=q, blue=True,\
            label=r'Median', color='red', right_axis=True,\
            csv_fnam=workdir + "csv/%sR.csv" % (out_fnam,))
        ax.set_title("%s: Infected population,\nf=%4.2f, N=%.2e" % (zm.Region,zm.f[0],zm.N_org))
        fig.tight_layout()
        fig.savefig("%s%s_R.png" % (workdir + 'figs/',out_fnam))

        fig, ax = subplots( num=8, figsize=(8,6))
        zm.PlotOmega(ax=ax)
        PlotFrozenDist(beta( 1+1/6, 1+1/3), ax=ax)
        ax.set_title("%s: Attack rate, f=%4.2f" % (zm.Region,zm.f[0]))
        fig.tight_layout()
        fig.savefig("%s%s_Omega.png" % (workdir + 'figs/',out_fnam))

    return ax



# csv files with data need to be in:
# workdir/data/id.csv
# where id is the identifiacion of the metro zone to be analyzed.
# The first column is daily deaths and the second daily confirmed cases
# The rest of the information is contained in a doctionary:

# General information for the metro zone or region to be analyzed:
# The actual instance of the covid_mcmc object is store in the last item of the list, once it has been instantiated by AnalyzeZM
#     id          Name   num_relax_days  Population   init date           intervention and relax dates
ZMs = $Zone
# The (optional) hospital occupancy data are stored in:
#  workdir/data/hosp/


rcParams.update({'font.size': 14})



if __name__=='__main__':

    q = $plottingQuantiles  # [10, 25, 50, 75, 90]

    workdir = "$workdir"
    os.makedirs(workdir, exist_ok=True)
    os.makedirs(workdir+"output", exist_ok=True)
    os.makedirs(workdir+"figs", exist_ok=True)
    os.makedirs(workdir+"data", exist_ok=True)
    #os.makedirs(workdir+"data/hospitales_datos", exist_ok=True)
    os.makedirs(workdir+"csv", exist_ok=True)
    #os.makedirs(workdir+"logs", exist_ok=True)
    #os.makedirs(workdir+"latex", exist_ok=True)
    #os.makedirs(workdir+"latex/images", exist_ok=True)


    # This a command line option (e.g. to run in batch):
    #                                T       id  pred
    # python AnalysisMZ.py --batch 200000  9-01  150
    if len(sys.argv) > 1:
        use('Agg') # Use ths backend, for ploting in batch
        if sys.argv[1] == '--batch':
            T=int(sys.argv[2])
            clave = sys.argv[3]
            pred = int(sys.argv[4])
            ZMs[clave][-1] = AnalyzeZM( clave, T=T, pred=pred, plot_fit=False)
    else:

        # This is to run interactively
        T= $T # 20000, Usual number of iterations used for MCMC
        burnin = $burnin
        clave = $Zone_id
        ZMs = $Zone
        trim = $trim
        init = $init_date  #init=date(2020, 4, 1)

        print( "\n%s, %s, T=%d, trim=%d" % ( clave, ZMs[clave][0], T, trim))
        ZMs[clave][-1] = AnalyzeZM( clave, T=T, trim=trim, init=init, pred=$pred, plot_fit=False, burnin=burnin, workdir=workdir)
        ZMs[clave][-1].clave = clave
        zm = ZMs[clave][-1]
        # Custom plots for a metro zone
        PlotFigsZMs(ZMs[clave][-1], pred=$plotpred, blue=True, workdir=workdir)
