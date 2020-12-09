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

from . import victoria



def odeint( rhs, X0, t_quad, args):
    """
    See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

    Parameters
    ----------
    rhs : callable(y, t, …)
        Computes the derivative of y at t.
    X0 : array
        Initial condition on y (can be a vector).
    t_quad : array
        A sequence of time points for which to solve for y. The initial value point should be the first element of this sequence. This sequence must be monotonically increasing or monotonically decreasing; repeated values are allowed.
    args : tuple
        Extra arguments to pass to function.


    Returns::

        scipy.integrate.odeint( rhs, X0, t_quad, args)



    """
    return integrate.odeint( rhs, X0, t_quad, args)
    #return rungekuttaodeint( rhs, X0, t_quad, h=1/15, args=args, Method=3)[0]

def Model_SEID( m_list, e_list, prn=True):
    """
    Define the graph matrix describing the model.

    .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: Model_SEID
            :lines: 1,10-

    """


    T = victoria.AuxMatrix(names="S E I D", prn=prn)

    T.BaseVar('S')

    T.Exit( 'S', 'E')
    T.Exit( 'E', 'I')
    T.Exit( 'I', 'D')
    T.NoExit('D')

    T.End()

    # Split in Erlang series of length m
    T.SplitErlang( e_list, m_list)
    return T


class SEID(victoria.mcmc):
    """
    ODE model of a SEID model, and parameter inference using
    Bayesian inference (MCMC with the twalk).

    Args:
        Region:      Region
        N:           population size
        data_fnam:   Data file name, workdir/data/data_fnam
                    Make sure to process data into a vertical text array
        out_fnam:    MCMC output file name, without .txt, workdir/output/out_fnam + '.txt'
        init_index:  day number from data start, where to start the plot
        init:        date of init_index

        R_rates:     Residence rates
        trim:        how many data to trim
    """
    def __init__( self, Region, N, data_fnam, out_fnam,\
                 init_index, init, R_rates, trim, workdir="./../"):


        super().__init__(Region=Region, N=N, data_fnam=data_fnam, out_fnam=out_fnam,\
                 init_index=init_index, init=init, trim=trim, workdir=workdir)


        self.Init_fm_matrix( m=self.data.shape[0], R_rates=R_rates)

    def Init_fm_matrix( self, m, R_rates):
        """
        Init the forward map.

        Args:
            m (int): is the number of observed days, ie. sample size.
            R_rates (dict): dictionary of residence rates, `R_rates={ 'E'  :[1/1.5, r'\sigma_1'], 'I^S':[1/2  , r'\sigma_2']}`
        """
        self.num_pars = 2 # Number of parameters to be inferred: S(0) and contact rate

        self.R_rates = R_rates
        m_list=[self.R_rates[v][2] for v in self.R_rates.keys()]
        e_list=list(self.R_rates.keys())

        # Define the graph matrix describing the model (see above)
        self.T = Model_SEID( m_list, e_list, prn=True)
        # The rhs will be ('@' matrix multiplication of arrays in Python 3):
        #Graph matrix     State vars     Erlang mask      par (see below)
        #rhs(x)= M      @       x      *       E        @  par
        #       qxq            qx1            qxn          nx1
        # Before this, par is filled with non-linear terms etc. like the force of infection
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
        self.mask_S   = self.T.SelectMask('S')
        """Masks to select variables from list of state variables."""
        self.mask_E   = self.T.SelectMask('E')
        """Masks to select variables from list of state variables."""
        self.mask_I  = self.T.SelectMask('I')
        """Masks to select variables from list of state variables."""
        self.mask_D   = self.T.SelectMask('D')
        """Masks to select variables from list of state variables."""

        self.X0 = np.zeros((self.q,))

    def GetMask( self, v, E_range='all', as_col_vec=False):
        """Returns a mask to select variable `v` from `grp_list`.

           E_range = [0] (default), first in the list, or original variable if no Erlang list.
           E_range = 'all' use the whole Erlang list for variable
           or provide E_range list manually.
        """
        return self.T.SelectMask( v, E_range=E_range, as_col_vec=as_col_vec)

    def rhs( self, x, t, p):
        """ See :meth:`victoriaepi.seid.odeint` to check usage.

        .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: SEID.rhs
            :lines: 1,9-
        """

        beta = p[1] #I(0) is p[0]

        I = np.sum(x * self.mask_I) # total number of infectious

         #force of infection beta1*I^A/N + some factor of  beta1*I^S/N
        foi = I/self.N * beta

        self.par[self.T.ConvertVar('S')] = foi

        return self.T.M @ (x * (self.T.par_mask @ self.par))

    def solve_plain(self, p, quad=True):
        """
        Solve the initial value problem.

        .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: SEID.solve_plain
            :lines: 1,9-
        """

        self.X0 *= 0
        self.X0 += p[0]*self.mask_I #Initial infected
        self.X0 += (self.N-p[0])*self.mask_S #suceptible
        if quad:
            return odeint(self.rhs, self.X0, self.t_quad, args=(p,))
        else:
            return odeint(self.rhs, self.X0, self.time, args=(p,))

    def solve( self, p):
        """
        Solve the initial value problem.
        Integral of incidence between observation times

        .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: SEID.solve
            :lines: 1,10-
        """
        # Use the solver:

        self.soln = self.solve_plain( p, quad=False )

        return [np.diff(self.soln[::self.quad_k,:] @ self.mask_D)] # list of size self.p=1

    def llikelihood( self, p):
        """
        Log likelihood.

        .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: SEID.llikelihood
            :lines: 1,9-
        """
        #if support(p): # Not necessary, the twalk checks it already before acllin energy support(p):
        # negative binomial likelihood
        mu_D = self.solve(p)[0]
        mu_D +=3
        # negative binomial likelihood for deaths
        omega = 2.0
        theta = 0.5   #antonio 0.5
        r = mu_D/(omega-1.0+theta*mu_D)
        q = 1.0/(omega+theta*mu_D)
        log_likelihood = np.sum(ss.nbinom.logpmf( self.data+3,r,q))

        return log_likelihood

    def lprior( self, p):
        """
        Log prior.

        .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: SEID.lprior
            :lines: 1,9-
        """
        # Log priors:
        log_prior = 0.0
        # gamma prior distribution parameters for I(0)
        log_prior += ss.gamma.logpdf(p[0],1.0,scale=10.0)

        # log-normal prior distribution parameters for beta
        log_prior += np.sum(ss.lognorm.logpdf(p[1], 1.0, scale=1.0)) #scale=np.exp(0.0)

        return log_prior

    def support( self, p):
        """
        Support.

        .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: SEID.support
            :lines: 1,10-
        """

        rt = True
        rt &= (0.0 < p[0] < 10.0**2)
        # beta in [0,20]
        rt &= all((0.0 < p[1]) * (p[3:] < 20.0))

        return rt

    def sim_init(self):
        """Simulate initial values for mcmc.

        .. literalinclude:: ../victoriaepi/seid.py
            :pyobject: SEID.sim_init
            :lines: 1,8-
        """
        p = np.zeros(self.num_pars)
        p[0] = np.random.uniform(low = 0.01, high = 10.0)
        p[1] = np.random.uniform(low = 0.01, high = 5.0)

        return p

    def PlotEvolution( self, pred, cumm=False, log=False, ax=None,\
                       csv_fnam=None, q=[ 10, 25, 50, 75, 90], blue=True, add_MRE=False,\
                       color='red', color_q='black', label='Mediana', right_axis=True, label_cases=True):
        """ Plot Evolution.

        Args:

            pred: number of days to predict
            ty: 0 = Infected,1 = deaths
            cumm: True if cumulative, default False
            log: True if y log scale, default False
            ax: axis where to print the plot (optional)
            csv_fnam: name of file to save the data for the plot (optional)

        """
        if ax == None:
            fig = plt.figure(figsize=(12,10))
            ax = fig.gca()
        else:
            fig = None

        data = self.data # Deaths REPORTED
        data_trimed = self.data_trimed
        title = 'Deaths'

        # cumulative or prevalanece, prepapr solns
        if cumm:
            prevalence = np.cumsum(data) # aggregate observed data
            self.future = prevalence[-1] + np.cumsum(data_trimed)
            solns = self.solns[0]
            ylabel = 'Accumulated cases'
            title = 'Accumulated ' + title
        else:
            prevalence = data # aggregate observed data
            self.future = data_trimed
            solns = np.diff( np.append( np.zeros((self.solns[0].shape[0],1)), self.solns[0], axis=1), axis=1)
            ylabel = 'Num. cases'
            title = 'Incidence of ' + title

        self.PlotEvolution_fm( solns=solns, prevalence=prevalence, pred=pred, ylabel=ylabel, log=log, ax=ax,\
                       csv_fnam=csv_fnam, q=q, blue=blue, add_MRE=add_MRE,\
                       color=color, color_q=color_q, label=label, right_axis=right_axis, label_cases=label_cases)

        ax.set_title(self.Region + '. ' + title)

