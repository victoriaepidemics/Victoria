#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 20:32:32 2020

@author: jac
"""


import sys

import numpy as np
import scipy.stats as ss
from scipy import integrate
import matplotlib.pyplot as plt


from victoriaepi import victoria



def odeint( rhs, X0, t_quad, args):
    return($odeint)



$defineModelMatrix



class $className(victoria.mcmc):
    def __init__( self, Region, N, data_fnam, out_fnam,\
                 init_index, init, R_rates, trim, workdir="./../"):

        """
            Simple SEIR model, and parameter inference using
            Bayesian inference (MCMC with the twalk).

        Region      # Region
        N           # population size
        data_fnam   # Data file name, workdir/data/data_fnam
                    # Make sure to process data into a vertical text array
        out_fnam    # MCMC output file name, without .txt, workdir/output/out_fnam + '.txt'
        init_index  # day number from data start, where to start the plot
        init        # date of init_index

        R_rates     # Residence rates
        trim        # how many data points to trim
        """

        super().__init__(Region=Region, N=N, data_fnam=data_fnam, out_fnam=out_fnam,\
                 init_index=init_index, init=init, trim=trim, workdir=workdir)


        self.Init_fm_matrix( m=self.data.shape[0], R_rates=R_rates)

    def Init_fm_matrix( self, m, R_rates):
        """Init the forward map:
            m is the number of observed days, ie. sample size.
            R_rates, dictionary of residence rates,
             e.g. `R_rates={ 'E'  :[1/1.5, r'\sigma_1'], 'I^S':[1/2  , r'\sigma_2']}`
        """
        self.num_pars = $num_pars # Number of parameters to be inferred: S(0) and contact rate

        self.R_rates = R_rates
        m_list=[self.R_rates[v][2] for v in self.R_rates.keys()]
        e_list=list(self.R_rates.keys())

        # Define the graph matrix describing the model (see above)
        self.T = Model_$className( m_list, e_list, prn=True)
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
        # p=1 size of the return list for solve, daily deaths
        # quad_k=1 number of subdivions in day, no quadrature is needed
        super().Init_fm_matrix( m, num_state_vars=self.q, p=1, quad_k=1)

        self.result_D = np.zeros(self.nn-1) # To hold result of quadrature

        # ""S E I D""
        self.par = np.zeros(self.n) # Par mask

        self.R_rates = R_rates
        # Known parameters, set residence rates:
        for v in R_rates.keys():
            self.par[self.T.ConvertVar(v)] = R_rates[v][0]

        # The masks to select variables from list of state variables
        $SelfMaskCode

        self.X0 = np.zeros((self.q,))

    def GetMask( self, v, E_range='all', as_col_vec=False):
        """Returns a mask to select variable v from grp_list.

           E_range = [0] (default), first in the list, or original variable if no Erlang list.
           E_range = 'all' use the whole Erlang list for variable
           or provide E_range list manually.
        """
        return self.T.SelectMask( v, E_range=E_range, as_col_vec=as_col_vec)

    $solvingMethods

    $plottingMethods

    $additionalMethods

