#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 13:22:09 2020



ama2 model:

Marcos Capistran, Antonio Capella y Andres Christen,
Marzo-Octubre 2020.


"""

import sys

import numpy as np
import scipy.stats as ss
from scipy import integrate
import matplotlib.pyplot as plt

from victoriaepi import victoria
from victoriaepi import pytwalk
import pickle



def odeint( rhs, X0, t_quad, args):
    return($odeint)

$defineModelMatrix


class $className(victoria.mcmc):
    def __init__( self, Region, N, intervention_day, data_fnam, out_fnam,\
                 Pobs_I, Pobs_D,\
                 init_index, init, exit_probs, R_rates, trim, relax_day=[], workdir="./../"):
        """
            ODE model, and parameter inference using
            Bayesian inference (MCMC with the twalk).

        N           # population size
        data_fnam   # Data file name, workdir/data/data_fnam
                    # Make sure to process data into a vertical text array
        out_fnam    # MCMC output file name, without .txt, workdir/output/out_fnam + '.txt'
        init_index  # day number from data start, where to start the plot
        init        # date of init_index
        exit_probs  # Exit probabilities of graph
        Pobs_I, Pobs_D # Prob of observation

        R_rates     # Residence rates
        trim        # how many data to trim

        intervention_day # None, single value or list in ascending order
        relax_day # LIST (always) of relaxation days, default []
        """

        super().__init__(Region=Region, N=N, data_fnam=data_fnam, out_fnam=out_fnam,\
                 init_index=init_index, init=init, trim=trim, workdir=workdir)

        # intervention_day: None, single value or list in ascending order
        if intervention_day == None:
            intervention_day = np.array([10000])
        elif isinstance(intervention_day, list):
            if relax_day == []:
                intervention_day = np.array(intervention_day + [10000])
            else: # Add a new beta at the intervetion day
                if np.max(intervention_day) > min(relax_day):
                    print("relax_day %d most be larger that intervention days %d",\
                          (relax_day, np.max(intervention_day)))
                    raise
                intervention_day = np.array(intervention_day + relax_day + [10000])
        else: # Single value
            if relax_day == []:
                intervention_day = np.array([intervention_day] + [10000])
            else: # Add a new beta at the intervetion day
                if np.max(intervention_day) > min(relax_day):
                    print("relax_day %d most be larger that intervention days %d",\
                          (relax_day, np.max(intervention_day)))
                    raise
                intervention_day = np.array([intervention_day] + relax_day + [10000])
        # Add a distant day in the future, so after last intervention day
        # np.where(t < self.intervention_day)[0][0] still works

        self.Init_fm_matrix( m = self.data.shape[0],\
                          intervention_day = intervention_day,\
                          exit_probs = exit_probs, R_rates = R_rates)

        if isinstance( Pobs_I, float):
            self.Pobs_I = Pobs_I           # probability of recording an infection
            self.Pobs_D = Pobs_D           # probability of recording an deaths
        else: #Nowcasting function, forming diferent obs probabilities, weighted average for nowcasting recent obs subcountns, ie. Pobs>1
            self.Pobs_I = Pobs_I(np.arange(self.m+1,1,-1))
            self.Pobs_D = Pobs_D(np.arange(self.m+1,1,-1))

        self.num_betas = self.num_pars - 3

        self.relax_day = relax_day
        if self.relax_day == []: #ama1
            self.solve_plain = self.solve_plain1 # Old without relaxation day
        else: #ama2
            self.solve_plain = self.solve_plain2 # with relaxation day
            self.num_pars += 1 + len(relax_day) # plus \omega_0 and \omega_i, last parameters
            self.num_omegas = 1 + len(relax_day)
            # Iterator for omegas:
            self.omega_ite = range(-self.num_omegas,0,1)
            self.N_org = self.N
            self.relax_day += [10000] #Dummy distant day




    def Init_fm_matrix( self, m, intervention_day, exit_probs, R_rates):
        """Init the forward map:
            m is the number of observed days, ie. sample size.
           intervention_day, number from day 0 in data
           f, g, h, i = exit_probs
           R_rates, dictionary of residence rates,
             e.g. `R_rates={ 'E'  :[1/1.5, r'\sigma_1'], 'I^S':[1/2  , r'\sigma_2']}`
        """
        self.num_pars = $num_pars + len(intervention_day) # Number of parameters to be inferred

        self.intervention_day = intervention_day
        self.factor_foi = $factor_foi  # 0.1*20
        self.R_rates = R_rates
        m_list=[self.R_rates[v][2] for v in self.R_rates.keys()]
        e_list = list(self.R_rates.keys())
        # Known parameters, set exit probs
        $exitProbsDefinitions

        # Define the graph matrix describing the model (see above)
        self.T = Model_$className(m_list, e_list, $inlineargsequals prn = False)
        # The rhs will be ('@' matrix multiplication of arrays in Python 3):
        #Graph matrix     State vars     Erlang mask      par (see below)
        #rhs(x)= M      @       x      *       E        @  par
        #       qxq            qx1            qxn          nx1
        # Before this, par is filled with non-linear terms etc. like the forze of infection
        # n original number of state variables
        # q final number of state variables after Erlang series


        self.n = self.T.n
        self.q = self.T.q # Total number of state variables

        # Call the base class Init_fm_matrix
        # p = 2 size of the return list for solve, incidence and deths
        # quad_k = 10 number of subdivions in day for quadrature
        super().Init_fm_matrix( m, num_state_vars = self.q, p = 2, quad_k = 10)

        # "S E I^A I^S I^B H^1 H^2 H^3 U^1 U^2 D R"
        self.par = np.zeros(self.n) # Par mask

        self.R_rates = R_rates
        # Known parameters, set residence rates:
        for v in R_rates.keys():
            self.par[self.T.ConvertVar(v)] = R_rates[v][0]

        # Auxiliars to be used in solve
        $AuxiliarDefinitionsForSolve

        # The masks to select variables from list of state variables
        $SelfMaskCode
        #self.psi = 10
        #IA_k = np.where(self.mask_IA == 1)[0][0]
        #self.T.M[0,IA_k] = -self.psi/(1+self.psi)

        self.X0 = np.zeros((self.q,))

    $solvingMethods

    $plottingMethods

    $additionalMethods