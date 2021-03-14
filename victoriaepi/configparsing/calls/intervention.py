#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Marcos Capistran, Antonio capella y Andres Christen,
Marzo-Octubre 2020.

"""

import sys
import os
from datetime import date, timedelta

from numpy import array, random
from matplotlib import use
from matplotlib.pyplot import subplots, rcParams, close

from pandas import datetime, read_csv



from $moduleName import $className


"""Configuration parameters."""

# Exit probabilities
exit_probs = $exit_probs


out_f = exit_probs[0]

Pobs_I = $Pobs_I
Pobs_D = $Pobs_D

R_rates_V2 = $R_rates

def AnalyzeZM( clave, T, trim=0, pred=100, init=date(2020, 4, 1), plot_fit=True, workdir="./../", burnin=500):

    # Transfer and standardization of reference date for graphing
    init = init
    if (ZMs[clave][3]+ timedelta(days=$daysdelay)) <= init:
        # Including the 4-day shift (symptoms onset)
        init_index=(init - (ZMs[clave][3] + timedelta(days=$daysdelay))).days
    else:
        #init = ZMs[clave][3]
        init=ZMs[clave][3] + timedelta(days=$daysdelay)
        init_index = 0

    if trim == 0:
        out_fnam = "$out_fnam"
    else:
        out_fnam = "$out_fnam" + "_trim%d" % (trim,)  # clave+"_trim%d" % (trim,)
    data_fnam = "$data_fnam"

    intervention_day = [(d-ZMs[clave][3]).days for d in ZMs[clave][4:-(int(ZMs[clave][1])+1)]]

    exit_probs_copy = exit_probs.copy()
    if int(ZMs[clave][1]) > 0:
        relax_day = [(d - ZMs[clave][3]).days for d in ZMs[clave][(-(int(ZMs[clave][1])+1)):-1]]
        exit_probs_copy[0] = out_f#0.4 #=f
        print("ama-2")
    else:
        relax_day = []

    N = ZMs[clave][2]

    zm = $className( Region=ZMs[clave][0], data_fnam=data_fnam,\
        N=N, out_fnam=out_fnam, init_index=init_index, init=init,\
        intervention_day=intervention_day, relax_day=relax_day, trim=trim, Pobs_I=Pobs_I, Pobs_D=Pobs_D,\
        exit_probs=exit_probs_copy, R_rates=R_rates_V2, workdir=workdir)

    if T > 0:
        zm.RunMCMC(T=T, pred=pred, plot_fit=plot_fit, burnin=burnin)
    return zm





# Do a series of costomized plots using PlotEvolution, PlotStateVar and PlotOmega
def PlotFigsZMs( zm, pred=99, q=[10,25,50,75,90], blue=True, workdir='./../'):

    close('all')
    dateparse = lambda x: datetime.strptime(x, '%Y-%m-%d') #To read the IRAG csv's

    try:
        zm_vmx_Hosp_RI = read_csv(workdir + "data/hosp/%s_DinHospitales.csv" % (zm.clave,), parse_dates=['fecha'], date_parser=dateparse)
        hosp = True
    except:
        hosp = False

    out_fnam=zm.out_fnam

    fig, ax = subplots(figsize=(8,6))
    zm.PlotEvolution( pred=4, cumm=False, log=False, ty=0, ax=ax,\
            label='Median', q=q, blue=blue, right_axis=False,\
            csv_fnam=workdir + "csv/%s_I_short.csv" % (out_fnam,))
    ax.set_title('%s' % (zm.Region,))
    ax.set_xlabel("Date (day.month)")
    ax.set_ylabel("Confirmed cases")
    fig.tight_layout()
    fig.savefig("%s%s_I_short.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots(figsize=(8,6))
    zm.PlotEvolution( pred=int(pred/2), cumm=False, log=False, ty=1, ax=ax,\
            label='Median', q=q, blue=blue, right_axis=False,\
            csv_fnam=workdir + "csv/%s_D_short.csv" % (out_fnam,))
    ax.set_title('%s' % (zm.Region,))
    ax.set_xlabel("Date (day.month)")
    ax.set_ylabel("Deaths")
    fig.tight_layout()
    fig.savefig("%s%s_D_short.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots(figsize=(8,6))
    zm.PlotEvolution( pred=pred, cumm=False, log=False, ty=0, ax=ax, q=q, blue=blue, right_axis=False,\
        label=r'Median', csv_fnam=workdir + "csv/%s_I_long.csv" % (out_fnam,))
    ax.set_title('%s' % (zm.Region,))
    ax.set_xlabel("Date (day.month)")
    ax.set_ylabel("Confirmed cases")
    fig.tight_layout()
    fig.savefig("%s%s_I_short.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots(figsize=(8,6))
    zm.PlotEvolution( pred=int(pred/2), cumm=True, log=False, ty=1, ax=ax, q=q, blue=blue,\
        label=r'Median', right_axis=False, csv_fnam=workdir+"csv/%s_D_long.csv" % (out_fnam,))
    ax.set_title('%s' % (zm.Region,))
    ax.set_xlabel("Date (day.month)")
    ax.set_ylabel("Accumulated deaths")
    fig.tight_layout()
    fig.savefig("%s%s_D_long.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots(figsize=(8,6))
    zm.PlotStateVar( 'H^1 H^2 H^3', pred=pred, ax=ax, q=q, blue=blue,\
        label=r'Median', color='red', right_axis=False,\
        csv_fnam=workdir + "csv/%s_Hs.csv" % (out_fnam,))
    if hosp:
        ax.bar(  array(zm_vmx_Hosp_RI['fecha']), array(zm_vmx_Hosp_RI['camasOcupadasHG']), color='red', width=0.5, alpha=0.5)
        ax.plot( zm_vmx_Hosp_RI['fecha'], zm_vmx_Hosp_RI['camasOcupadasHG'], 'ro', markersize=2)
    ax.set_title('%s' % (zm.Region,))
    ax.set_xlabel("Date (day.month)")
    ax.set_ylabel("Occupancy no ICU beds")
    fig.tight_layout()
    fig.savefig("%s%s_Hs.png" % (workdir + 'figs/',out_fnam))

    fig, ax = subplots(figsize=(8,6))
    zm.PlotStateVar( 'U^1', pred=pred, ax=ax, q=q, blue=blue,\
        label=r'Median', color='red',right_axis=False,\
        csv_fnam=workdir + "csv/%s_U1.csv" % (out_fnam,))
    if hosp:
        ax.bar( array(zm_vmx_Hosp_RI['fecha']), array(zm_vmx_Hosp_RI['camasOcupadasUCI_NoUCIVent']), color='red', width=0.5, alpha=0.5)
        ax.plot( zm_vmx_Hosp_RI['fecha'], zm_vmx_Hosp_RI['camasOcupadasUCI_NoUCIVent'], 'ro', markersize=1)
    ax.set_title('%s' % (zm.Region,))
    ax.set_xlabel("Date (day.month)")
    ax.set_ylabel("Occupancy ICU beds")
    fig.tight_layout()
    fig.savefig("%s%s_U1.png" % (workdir + 'figs/',out_fnam))

    # Optional, recovered (see paper):
    fig, ax = subplots( num=7, figsize=(8,6))
    zm.PlotStateVar( 'R', pred=pred, ax=ax, q=q, blue=blue,\
        label=r'Median', color='black', right_axis=False)
    ax.set_title('')
    ax.set_xlabel("Date (day.month)")
    ax.set_ylabel("Accumulated cases")
    ax.legend().set_visible(False)
    fig.tight_layout()
    fig.savefig("%s%s_R.png" % (workdir + 'figs/',out_fnam))

    if zm.relax_day != []:
        fig, ax = subplots( num=8, figsize=(8,6))
        zm.PlotOmega(ax=ax)
        ax.set_title('')
        ax.set_ylabel("Density")
        fig.tight_layout()
        fig.savefig("%s%s_Omega.png" % (workdir + 'figs/',out_fnam))

        fig, ax = subplots( num=9, figsize=(8,6))
        zm.PlotOmega( ax=ax, sel=(0,1), mul_f=True)
        ax.set_title('')
        ax.set_ylabel("Density")
        fig.tight_layout()
        fig.savefig("%s%s_Omega_f.png" % (workdir + 'figs/',out_fnam))

    return fig, ax



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

    q = $plottingQuantiles  # [10,25,50,75,90] # Define the quantiles for plotting

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
    if len(sys.argv) > 5:
        use('Agg') # Use ths backend, for ploting in batch
        if sys.argv[1] == '--batch':
            T=int(sys.argv[2])
            clave = sys.argv[3]
            pred = int(sys.argv[4])
            ZMs[clave][-1] = AnalyzeZM( clave, T=T, pred=pred, plot_fit=False)
    else:


        T= $T

        burnin = $burnin

        clave = $Zone_id


        ZMs = $Zone

        print(ZMs[clave], flush=True)
        print( "\n" + clave + " T=", T, " burnin=", burnin, flush=True)

        trim = $trim
        init = $init_date
        #pred = $pred
        print( "\n%s, %s, T=%d, trim=%d" % ( clave, ZMs[clave][0], T, trim))
        ZMs[clave][-1] = AnalyzeZM( clave, T=T, trim=trim, init=init, pred=$pred, plot_fit=False, burnin=burnin, workdir=workdir)
        ZMs[clave][-1].clave = clave
        zm = ZMs[clave][-1]
        # Custom plots for zone
        PlotFigsZMs(ZMs[clave][-1], pred=$plotpred, blue=False, workdir=workdir)
