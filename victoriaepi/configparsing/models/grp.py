#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 13:06:40 2020

@author: jac
"""


import sys
import pickle
import datetime as dt

import numpy as np
import scipy.stats as ss
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


from victoriaepi import pytwalk

from victoriaepi import victoria




def odeint( rhs, X0, t_quad, args):
    return($odeint)

$defineModelMatrix


class $className(victoria.mcmc):
    def __init__( self, Region, N, intervention_day, data_fnam, out_fnam,\
                 Pobs_I, Pobs_D,\
                 init_index, init, trim, R_rates,\
                 ngrp, exit_probs, Int_M, age_prop,\
                 relax_day=[], workdir="./../"):
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
        trim        # how many data to trim
        Pobs_I, Pobs_D # Prob of observation

        R_rates     # Residence rates
        ngrp number of groups
        f, g, h, i = exit_probs, each a list of length ngrp
        Int_M ngrp X ngrp array with the interaction
        age_prop list of age proportions of length ngrp

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
                    print("relax_day %d must be larger that intervention days %d",\
                          (relax_day, np.max(intervention_day)))
                    raise
                intervention_day = np.array([intervention_day] + relax_day + [10000])
        # Add a distant day in the future, so after last intervention day
        # np.where(t < self.intervention_day)[0][0] still works

        self.Init_fm_matrix( m=self.data.shape[0],\
                          intervention_day=intervention_day,\
                          R_rates=R_rates,\
                          ngrp=ngrp, exit_probs=exit_probs, Int_M=Int_M, age_prop=age_prop)

        if isinstance( Pobs_I, float):
            self.Pobs_I = Pobs_I           # probability of recording an infection
            self.Pobs_D = Pobs_D           # probability of recording an deaths
        else: #Nowcasting function, forming diferent obs probabilities, weighted average for nowcasting recent obs subcountns, ie. Pobs>1
            self.Pobs_I = Pobs_I(np.arange(self.m+1,1,-1))
            self.Pobs_D = Pobs_D(np.arange(self.m+1,1,-1))

        self.num_betas = self.num_pars - 3

        self.relax_day = relax_day
        if self.relax_day == []: #ama1
            self.solve_plain = self.solve_plain1 # ama1, without relaxation day
        else: #ama2
            self.solve_plain = self.solve_plain2 # with relaxation day
            self.num_pars += 1 + len(relax_day) # plus \omega_0 and \omega_i, last parameters
            self.num_omegas = 1 + len(relax_day)
            # Iterator for omegas:
            self.omega_ite = range(-self.num_omegas,0,1)
            self.N_org = self.N
            self.relax_day += [10000] #Dummy distant day
        #


    def Init_fm_matrix( self, m, intervention_day, R_rates,\
                  ngrp, exit_probs, Int_M, age_prop):
        """m is the number of observed days, ie. sample size
           intervention_day, number from day 0 in data
           R_rates, dictionary of residence rates,
             e.g. `R_rates={ 'E'  :[1/1.5, r'\sigma_1'], 'I^S':[1/2  , r'\sigma_2']}`
           ngrp number of groups
           f, g, h, i = exit_probs, each a list of length ngrp
           Int_M ngrp X ngrp array with the interaction
           age_prop list of age proportions of length ngrp
        """

        self.num_pars = $num_pars + len(intervention_day) # Number of parameters to be inferred

        self.intervention_day = intervention_day

        self.factor_foi = $factor_foi  # 0.1*20

        self.intervention_day = intervention_day

        # Known parameters, set exit probs
        $exitProbsDefinitions
        #self.f, self.g, self.h, self.i = exit_probs
        #f           # fraction of severe infections
        #g   # fraction of severe infections that require hospitalization
        #h           # fraction of hospitations that require ICU
        #i           # fraction of ICU patients who die

        self.R_rates = R_rates
        m_list=[self.R_rates[v][2] for v in self.R_rates.keys()]
        e_list=list(self.R_rates.keys())

        self.ngrp = ngrp # Number of groups
        # Define the graph matrix describing the model
        # is the same graph, only some parameters may change (eg. f)
        self.Ts = [Model_$className(m_list, e_list, $inlineargsequals prn=False) for grp in range(self.ngrp)]
        # The rhs will be ('@' matrix multiplication of arrays in Python 3):
        #Graph matrix     State vars     Erlang mask      par (see below)
        #for grp in range(self.ngrp):
        # rt[(grp*self.q):((grp+1)*self.q)] = M      @       x      *       E        @  par
        #       qXq            qX1            qXn          nX1
        # Before this, par[grp] is filled with non-linear terms in its group
        # n original number of state variables in each
        # q final number of state variables after Erlang series
        # q*self.groups total number of state variables

        self.T = self.Ts[0] # To acces the shared graph etc.
        self.n = self.T.n
        self.q = self.T.q

        # Call the base class Init_fm_matrix
        # p=2 size of the return list for solve, incidence and deths
        # quad_k=10 number of subdivions in day for quadrature
        super().Init_fm_matrix( m, num_state_vars=self.q*self.ngrp, p=2, quad_k=10)

        # To hold the return of the rhs
        self.rt = np.zeros(self.q*self.ngrp)

        # List with each group indices:
        self.sel_ = [np.arange( grp*self.q, (grp+1)*self.q) for grp in range(self.ngrp)]

        # self.ngrp X self.ngrp interaction matrix
        self.Int_M = Int_M
        self.age_prop = age_prop

        # FOR THE MOMENT residence parameters are the same in each group
        # This can be easily changed (see rhs)
        # "S E I^A I^S I^B H^1 H^2 H^3 U^1 U^2 D R"
        self.par = np.zeros(self.n) # Par mask

        # Known parameters, set residence rates:
        for v in R_rates.keys():
            self.par[self.T.ConvertVar(v)] = R_rates[v][0]

        # Auxiliars
        $AuxiliarDefinitionsForSolve
        # self.sigma_1 = self.par[self.T.ConvertVar('E')] # Used in solve
        # self.gamma_1 = self.par[self.T.ConvertVar('I^A')]
        # self.sigma_2 = self.par[self.T.ConvertVar('I^S')]

        # The masks are the same for all groups since the graph is the same
        $SelfMaskCode
        # self.mask_S   = self.T.SelectMask('S')
        # self.mask_E   = self.T.SelectMask('E')
        # self.mask_IA  = self.T.SelectMask('I^A')
        # self.mask_IAs_flat = self.T.SelectMask('I^A', E_range='all')
        # self.mask_ISs = self.T.SelectMask('I^S', E_range='all', as_col_vec=True)
        # self.mask_ISs_flat = self.T.SelectMask('I^S', E_range='all')
        # self.mask_IS  = self.T.SelectMask('I^S')
        # self.mask_D   = self.T.SelectMask('D')

        self.X0 = np.zeros((self.q*self.ngrp,)) # Total of state variables



    def RunMCMC( self, T, burnin=1000, pred=100, plot_fit=True):
        """Run twalk MCMC, T = number of iterations.
           burnin, thining = IAT."""
        # This is a different RunMCMC, since grps are repetitions of a single model

        self.SetTime( 0 )
        self.twalk = pytwalk.pytwalk(n = self.num_pars, U = self.energy, Supp = self.support)
        self.twalk.Run( T=T, x0 = self.sim_init(), xp0 = self.sim_init())

        self.mcmc_samples = True

        self.iat = int(self.twalk.IAT(start=burnin)[0,0])
        self.burnin = burnin
        print("\nEffective sample size: %d" % ((T-burnin)/self.iat,))
        self.samples = self.twalk.Output[burnin::(self.iat),:] # Burn in and thining
        self.essize = self.samples.shape[0]

        self.SetTime( pred )
        # solutions I and D
        self.solns = [np.zeros(( self.essize, self.m + pred)) for i in range(self.p)]
        self.solns_plain = np.zeros(( self.essize, self.m + pred, self.q*self.ngrp))
        print("Sampling %d model solutions." % ( self.essize,))
        for index,m in enumerate(self.samples):
            tmp = list(self.solve(m[:-1]))
            self.solns_plain[index,:,:] = self.soln[10::10,:]
            for i, sl in enumerate(tmp):
                self.solns[i][index,:] = np.cumsum(sl)
            if ((index+1) % 100) == 0:
                print( index+1, end=' ')
        print("\nSaving files in ", self.workdir + 'output/' + self.out_fnam + '_*.pkl')
        pickle.dump( self.samples, open(self.workdir + 'output/' + self.out_fnam + '_samples.pkl', 'wb'))
        pickle.dump( self.solns, open(self.workdir + 'output/' + self.out_fnam + '_solns.pkl', 'wb'))
        pickle.dump( self.solns_plain, open(self.workdir + 'output/' + self.out_fnam + '_solns_plain.pkl', 'wb'))

        if plot_fit:
            self.PlotFit()


    $solvingMethods

    $plottingMethods

    $additionalMethods








