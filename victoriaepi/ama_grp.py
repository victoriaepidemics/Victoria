#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import sys
import pickle
import datetime as dt

import numpy as np
import scipy.stats as ss
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


from . import pytwalk

from . import victoria

from .ama import Model_ama



def odeint( rhs, X0, t_quad, args):
    """
    See:

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

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




class ama2_grp(victoria.mcmc):
    """
        ODE model of COVID19 pandemics, and parameter inference using
        Bayesian inference (MCMC with the twalk).

        Args:

            N           : population size
            data_fnam   :Data file name, workdir/data/data_fnam
                        Make sure to process data into a vertical text array
            out_fnam    :MCMC output file name, without .txt, workdir/output/out_fnam + '.txt'
            init_index  :day number from data start, where to start the plot
            init        :date of init_index
            exit_probs  :Exit probabilities of graph
            trim        :how many data to trim
            Pobs_I :Prob of observation
            Pobs_D :Prob of observation
            intervention_day: None, or single value or list in ascending order

            R_rates     :Residence rates
            ngrp: number of groups
            exit_probs: f, g, h, i each a list of length ngrp
            Int_M: ngrp X ngrp array with the interaction
            age_prop: list of age proportions of length ngrp
            relax_day (date): list of dates

    """
    def __init__( self, Region, N, intervention_day, data_fnam, out_fnam,\
                 Pobs_I, Pobs_D,\
                 init_index, init, trim, R_rates,\
                 ngrp, exit_probs, Int_M, age_prop,\
                 relax_day=[], workdir="./../"):


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
        """``self.num_betas = self.num_pars - 3``"""
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

    def PrintConfigParameters( self, f=sys.stdout):
        """Prints all configuration parameters to file f (default sys.stdout)."""
        keys=['Pobs_I', 'Pobs_D',\
                 'ngrp', 'f', 'g', 'h', 'i', 'Int_M', 'age_prop']
        print("\nR_rates:", file=f)
        for key in self.R_rates.keys():
            print("%4s, %7s = %6.4f (%4.1f d), E_list= %2d" %\
                  (key,self.R_rates[key][1],self.R_rates[key][0],1/self.R_rates[key][0],self.R_rates[key][2]),\
                  file=f)
        print("", file=f)
        for key in keys:
            print("%s :" % (key,), file=f)
            print(self.__dict__[key], file=f)
            print("", file=f)

    def Init_fm_matrix( self, m, intervention_day, R_rates,\
                  ngrp, exit_probs, Int_M, age_prop):
        """
        Args:
            m: is the number of observed days, ie. sample size
           intervention_day:  number from day 0 in data
           R_rates:  dictionary of residence rates, e.g. `R_rates={ 'E'  :[1/1.5, r'\sigma_1'], 'I^S':[1/2  , r'\sigma_2']}`
           ngrp: number of groups
           exit_probs: f, g, h, i  each a list of length ngrp
           Int_M: ngrp X ngrp array with the interaction
           age_prop: list of age proportions of length ngrp

        """

        self.num_pars = 3 + len(intervention_day) # Number of parameters to be inferred

        self.intervention_day = intervention_day

        self.factor_foi = 0.1*20

        self.intervention_day = intervention_day

        # Known parameters, set exit probs
        self.f, self.g, self.h, self.i = exit_probs
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
        self.Ts = [Model_ama( m_list, e_list, f=self.f[grp], g=self.g[grp], h=self.h[grp], i=self.i[grp], prn=False) for grp in range(self.ngrp)]
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
        self.sigma_1 = self.par[self.T.ConvertVar('E')] # Used in solve
        self.gamma_1 = self.par[self.T.ConvertVar('I^A')]
        self.sigma_2 = self.par[self.T.ConvertVar('I^S')]

        # The masks are the same for all groups since the graph is the same
        self.mask_S   = self.T.SelectMask('S')
        self.mask_E   = self.T.SelectMask('E')
        self.mask_IA  = self.T.SelectMask('I^A')
        self.mask_IAs = self.T.SelectMask('I^A', E_range='all')
        """``self.mask_IAs = self.T.SelectMask('I^A', E_range='all')``"""
        self.mask_ISs = self.T.SelectMask('I^S', E_range='all', as_col_vec=True)
        """``self.mask_ISs = self.T.SelectMask('I^S', E_range='all', as_col_vec=True)``"""
        self.mask_ISs_flat = self.T.SelectMask('I^S', E_range='all')
        """``self.mask_ISs_flat = self.T.SelectMask('I^S', E_range='all')``"""
        self.mask_IS  = self.T.SelectMask('I^S')
        self.mask_D   = self.T.SelectMask('D')

        self.X0 = np.zeros((self.q*self.ngrp,)) # Total of state variables

    def SetTime( self, shift=0):
        """Extends :meth:`victoriaepi.victoria.mcmc.SetTime`

        .. literalinclude:: ../victoriaepi/ama_grp.py
                :pyobject: ama2_grp.SetTime
                :lines: 1,8-
        """
        super().SetTime(shift=shift)
        self.result_I = np.zeros(self.nn-1) # To hold result of quadrature
        self.result_D = np.zeros(self.nn-1) # To hold result of quadrature
        self.result_H1 = np.zeros(self.nn-1) # To hold result of quadrature
        self.result_H1_aux = np.zeros(self.nn-1) # aux to claculate result of quadrature

    def GetMask( self, v, grps='all', E_range='all', as_col_vec=False):
        """Returns a mask to select variable v from grp_list.
           If grps = int or list of groups, if 'all', return the mask for all groups.
           E_range = [0] (default), first in the list, or original variable if no Erlang list.
           E_range = 'all' use the whole Erlang list for variable
           or provide E_range list manually.
        """
        if isinstance( grps, int):
            grps = [grps]
        if isinstance( grps, str): #'all' groups
            grps = range(self.ngrp)

        mask = self.T.SelectMask( v, E_range=E_range, as_col_vec=False)
        rt = np.zeros((self.q*self.ngrp,1))
        for grp in grps:
            rt[self.sel_[grp], 0] = mask #Repeat the mask for each group
        if as_col_vec:
            return rt
        else:
            return rt.flatten()

    def rhs( self, x, t, p):
        """ See :meth:`victoriaepi.ama_grp.odeint` to check usage.

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.rhs
            :lines: 1,9-
        """

        beta1 = p[3 + np.where(t < self.intervention_day)[0][0]]

        # total number of asymptomatic infections in each group
        I_A = np.array([np.sum(x[self.sel_[grp]] * self.mask_IAs) for grp in range(self.ngrp)])
        I_S = np.array([np.sum(x[self.sel_[grp]] * self.mask_ISs_flat) for grp in range(self.ngrp)])
        # array of length self.ngrps

        #force of infection in each group
        foi = beta1/self.N * (self.Int_M @ I_A + self.factor_foi*self.Int_M @ I_S)

        for grp in range(self.ngrp):
            self.par[self.Ts[grp].ConvertVar('S')] = foi[grp]

        for grp in range(self.ngrp):
            self.rt[self.sel_[grp]] = self.Ts[grp].M @ (x[self.sel_[grp]] * (self.Ts[grp].par_mask @ self.par))

        return self.rt

    def solve_plain1( self, p, quad=True):

        """
        Solve the initial value problem.

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.solve_plain1
            :lines: 1,9-
        """
        #beta_1 = p[0] # asymptomatic infection rate
        #beta_2 = p[1] # hospitalized infection rate

        for grp in range(self.ngrp):
            self.X0[self.sel_[grp]] *= 0
            self.X0[self.sel_[grp]] += p[0]*self.age_prop[grp]*self.mask_E
            self.X0[self.sel_[grp]] += p[1]*self.age_prop[grp]*self.mask_IA
            self.X0[self.sel_[grp]] += p[2]*self.mask_IS
            self.X0[self.sel_[grp]] += (self.N-p[0]-p[1]-p[2])*self.age_prop[grp]*self.mask_S

        if quad:
            return odeint(self.rhs, self.X0, self.t_quad, args=(p,))
        else:
            return odeint(self.rhs, self.X0, self.time, args=(p,))

    def solve_plain2(self, p, quad=True):
        """
        Solve the initial value problem with relaxation day

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.solve_plain2
            :lines: 1,9-
        """
        self.N = p[-self.num_omegas]*self.N_org #\omega_0*N_org, first N_eff

        for grp in range(self.ngrp):
            self.X0[self.sel_[grp]] *= 0
            self.X0[self.sel_[grp]] += p[0]*self.age_prop[grp]*self.mask_E
            self.X0[self.sel_[grp]] += p[1]*self.age_prop[grp]*self.mask_IA
            self.X0[self.sel_[grp]] += p[2]*self.mask_IS
            self.X0[self.sel_[grp]] += (self.N-p[0]-p[1]-p[2])*self.age_prop[grp]*self.mask_S

        if quad:
            rt = odeint(self.rhs, self.X0,\
                          self.t_quad[self.t_quad < self.relax_day[0]], args=(p,))
        else:
            rt = odeint(self.rhs, self.X0,\
                          self.time[self.time < self.relax_day[0]], args=(p,))

        for i in self.omega_ite[1:]:
            N = self.N
            self.N = p[i]*self.N_org # \omega_i*N_org = self.N_new > N = \omega_{i-1}*N_org
            X0 = rt[-1,:] # last values for state vars, initial conds after next relaxation day
            X0[0] += self.N - N #New amount of suceptibles: The remaining + the added bunch
            if quad:
                rt2 = odeint(self.rhs, X0,\
                          self.t_quad[(self.relax_day[i-1] <= self.t_quad) *  (self.t_quad < self.relax_day[i])], args=(p,))
            else:
                rt2 = odeint(self.rhs, X0,\
                          self.time[(self.relax_day[i-1] <= self.time) * (self.time < self.relax_day[i])], args=(p,))
            rt = np.append( rt, rt2, axis=0)

        return rt

    def solve( self, p):
        """
        Solve the initial value problem:
        Integral of incidence between observation times

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.solve
            :lines: 1,10-
        """
        # Use the solver:
        self.soln = self.solve_plain(p)

        self.result_H1 *= 0
        self.result_D *= 0
        # Integrate to produce accumulated number of infected
        for grp in range(self.ngrp):
            for k in range( 0, self.nn-1):
                #x_e = self.soln[ (10*k):(10*(k+1)+1),self.sel_[grp]] @ self.mask_E
                #incidence = self.f[grp]*self.sigma_1*x_e #To hospital
                #self.result_I_aux[k] = self.dt*np.dot(self.weigths,incidence)
                x_s = np.sum( self.soln[ (self.quad_k*k):(self.quad_k*(k+1)+1), self.sel_[grp]] @ self.mask_ISs, axis=1)
                incidence = self.g[grp]*self.sigma_2*x_s #To hospital
                self.result_H1_aux[k] = self.dt*np.dot(self.weigths,incidence)
            self.result_D += np.diff(self.soln[::self.quad_k, self.sel_[grp]] @ self.mask_D)
            self.result_H1 += self.result_H1_aux
        return self.result_H1, self.result_D

    def llikelihood( self, p):
        """
        Log likelihood.

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.llikelihood
            :lines: 1,9-
        """
        #if support(p): # Not necessary, the twalk checks it already before acllin energy support(p):
        # negative binomial likelihood
        mu_I, mu_D = self.solve(p)
        mu_I +=3
        mu_D +=3
        mu_I *= self.Pobs_I
        mu_D *= self.Pobs_D
        # negative binomial likelihood for infectious individuals
        omega = 2.0 #1
        theta = 1.0 #0.01
        r = mu_I/(omega-1.0+theta*mu_I)
        q = 1.0/(omega+theta*mu_I)
        log_likelihood_I = np.sum(ss.nbinom.logpmf( self.data[:,1]+3, r, q))
        # negative binomial likelihood for deaths
        omega = 2.0
        theta = 0.5   #antonio 0.5
        r = mu_D/(omega-1.0+theta*mu_D)
        q = 1.0/(omega+theta*mu_D)
        log_likelihood_D = np.sum(ss.nbinom.logpmf( self.data[:,0]+3,r,q))
        # joint likelihood
        log_likelihood = log_likelihood_D + log_likelihood_I

        return log_likelihood

    def lprior( self, p):
        """
        Log prior.

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.lprior
            :lines: 1,9-
        """
        # Log priors:
        log_prior = 0.0
        # gamma prior distribution parameters for E(0)
        log_prior += ss.gamma.logpdf(p[0],1.0,scale=10.0)
        # gamma prior distribution parameters for IA(0)
        log_prior += ss.gamma.logpdf(p[1],1.0,scale=10.0)
        # gamma prior distribution parameters for IS(0)
        log_prior += ss.gamma.logpdf(p[2],1.0,scale=10.0)



        # \beta_i and \omega_i's run in pairs
        # log-normal prior distribution parameters for beta_i
        log_prior += np.sum(ss.lognorm.logpdf(p[3], 1.0, scale=1.0)) #scale=np.exp(0.0)
        for i in range(1, self.num_betas):
            log_prior += np.sum(ss.lognorm.logpdf(p[3+i], 1.0, scale=p[3+i-1])) #scale=np.exp(log(p[i-1]))


        if self.relax_day != []:
            for i in self.omega_ite:
                log_prior += ss.beta.logpdf( p[i], 1+1/6, 1+1/3)

        return log_prior

    def support( self, p):
        """
        Support.

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.support
            :lines: 1,10-
        """
        rt = True
        rt &= (0.0 < p[0] < 10.0**3)
        rt &= (0.0 < p[1] < 10.0**3)
        rt &= (0.0 < p[2] < 10.0**3)
        # All betas (and omegas, to make it simple!) in [0,20]
        rt &= all((0.0 < p[3:]) * (p[3:] < 20.0))
        if not(rt):
            return False # Fast track, avoid checking the omegas, if any

        if self.relax_day != []:
            rt &= (0.0 < p[-self.num_omegas] < 1.0) #\omega0
            rt &= ( p[-self.num_omegas]*self.N_org-p[0]-p[1]-p[2] > 0.0)
            for i in self.omega_ite[:-1]:
                rt &= (p[i] < p[i+1] < 1.0) #\omega_i < \omega_{i+1}

        return rt

    def sim_init(self):
        """Simulate initial values for mcmc.

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.sim_init
            :lines: 1,8-
        """
        p = np.zeros(self.num_pars)
        p[0] = np.random.uniform(low = 0.01, high = 10.0)
        p[1] = np.random.uniform(low = 0.01, high = 10.0)
        p[2] = np.random.uniform(low = 0.01, high = 10.0)

        p[3:] = np.random.uniform(low = 0.01, high = 5.0, size=self.num_pars-3)


        if self.relax_day != []:
            p[-self.num_omegas] = ss.beta.rvs( 2.0, 4.0)
            for i in self.omega_ite[:-1]:
                p[i+1] = p[i] + (1-p[i])*ss.uniform.rvs()

        return p

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

    def PlotEvolution( self, pred, ty=0, cumm=False, log=False, ax=None,\
                       csv_fnam=None, q=[ 10, 25, 50, 75, 90], blue=True, add_MRE=False,\
                       color='red', color_q='black', label='Mediana', right_axis=True, label_cases=True):
        """Plot Evolution

            Args:

                pred: number of days to predict
                ty: 0 = Infected, 1 = deaths
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

        # Prepare data for comparisons
        if ty == 0:
            data = self.data[:,1] # Infected REPORTED
            data_trimed = self.data_trimed[:,1]
            title = 'Daily confirmed cases'
        else:
            data = self.data[:,0] # Deaths REPORTED
            data_trimed = self.data_trimed[:,0]
            title = 'Decesos'

        if isinstance(self.Pobs_I,float):
            Pobs = (self.Pobs_I, self.Pobs_D)
        else: #nowcasting, for plotting use only the limit nowcasting proportion, ie ignore the nowcasting
            Pobs = (self.Pobs_I[0], self.Pobs_D[0])

        # cumulative or prevalanece, prepapr solns
        if cumm:
            prevalence = np.cumsum(data) # aggregate observed data
            self.future = prevalence[-1] + np.cumsum(data_trimed)
            solns = Pobs[ty]*self.solns[ty]
            ylabel = 'Amount of accumulated cases'
            title = 'Acumulados de ' + title
        else:
            prevalence = data # aggregate observed data
            self.future = data_trimed
            solns = np.diff( np.append( np.zeros((self.solns[ty].shape[0],1)), Pobs[ty]*self.solns[ty], axis=1), axis=1)
            ylabel = 'Number of cases'
            title = 'incidencia de ' + title

        self.PlotEvolution_fm( solns=solns, prevalence=prevalence, pred=pred, ylabel=ylabel, log=log, ax=ax,\
                       csv_fnam=csv_fnam, q=q, blue=blue, add_MRE=add_MRE,\
                       color=color, color_q=color_q, label=label, right_axis=right_axis, label_cases=label_cases)

        # Put a marker at intervention days, if within plotting area
        for d in self.intervention_day[:-1]:
            if self.init_index <= d < self.init_index+self.shift:
                ax.axvline( self.days[d-self.init_index], color='black')
        # Put a red marker at relax date, if within plotting area
        for d in self.relax_day:
            if self.init_index <= d < self.init_index+self.shift:
                ax.axvline( self.days[d-self.init_index], color='red')

        ax.set_title(self.Region + ' ' + title)

    def PlotPanelEvolution( self, pred1, pred2, ty):
        """ Plot Panel Evolution.

        .. literalinclude:: ../victoriaepi/ama_grp.py
            :pyobject: ama2_grp.PlotPanelEvolution
            :lines: 1,9-
        """

        fig, axs = plt.subplots( nrows=2, ncols=2, sharex='col', figsize=(12,10), tight_layout=True)
        self.PlotEvolution( ax=axs[0,0], pred=pred1, ty=ty, cumm=False, log=False, every=1)
        self.PlotEvolution( ax=axs[1,0], pred=pred1, ty=ty, cumm=True, log=False, every=1)
        axs[1,0].set_title('')
        axs[1,0].get_legend().remove()

        self.PlotEvolution( ax=axs[0,1], pred=pred2, ty=ty, cumm=False, log=False, every=1)
        self.PlotEvolution( ax=axs[1,1], pred=pred2, ty=ty, cumm=True, log=False, every=1)
        axs[1,1].set_title('')
        axs[1,1].get_legend().remove()
        #fig.tight_layout()

    def PlotFit( self ):
        """Plots basic results of the MCMC."""
        if not(self.mcmc_samples):
            print("MCMC samples not available.")
            return

        fig, axs = plt.subplots( nrows=2, ncols=2, figsize=(9,6), tight_layout=True)
        self.PlotEvolution( ax=axs[0,0], pred=3, ty=0, cumm=False, log=False, every=1)
        self.PlotEvolution( ax=axs[1,0], pred=3, ty=1, cumm=False, log=False, every=1)
        axs[0,0].set_title(self.Region)
        axs[0,0].get_legend().remove()
        axs[0,0].set_ylabel('Confirmed cases')
        axs[1,0].set_title('')
        axs[1,0].get_legend().remove()
        axs[1,0].set_ylabel('Decesos')

        axs[0,1].plot(-self.twalk.Output[:,-1], 'k-')
        axs[0,1].plot(-self.twalk.Output[:self.burnin,-1], 'b-')
        axs[0,1].set_ylabel('LogPost')
        axs[1,1].plot(-self.samples[:,-1], 'k-')
        axs[1,1].set_ylabel('LogPost')
        fig.tight_layout()

    def PlotStateVar( self, v, pred, grps='all', E_range = 'all', log=False, every=1, title=None,\
                     ax=None, median_only=False, q=[ 10, 25, 50, 75, 90], blue=True,\
                     color='red', color_q='black', label=' ', right_axis=True, csv_fnam=None):
        """ Plot evolution of state var(s) in grps.
            If v are many state vars then sum their values.

            eg.::

                PlotStateVar('H^2', pred=20)
                PlotStateVar('H^2 H^3', pred=60)

            Args:


                grps: may be a list or a number.
                    If grps is a list, add the values of all state vars in all grps.
                    If grps = 'all' take all groups

                pred: number of days to predict

                E_range: What Erlang list to include in variable, default 'all'.
                log: True if y log scale, default False
                every: prediction boxes given 'every' days, default 1
                ax(optional): axis where to print the plot
                csv_fnam (optional): name of file to save the data for the plot

        """
        # This is different from single grp:

        if pred > self.solns[0].shape[1] - self.m:
            print("Maximum prediction days pred=%d" % (self.solns_plain.shape[1] - self.m))
            return

        mask = np.zeros((self.q*self.ngrp,1))
        for var in v.split():
            mask += self.GetMask( var, grps=grps, E_range=E_range, as_col_vec=True)

        if isinstance( grps, int):
            grps = [grps]
        if isinstance( grps, str): #'all' groups
            grps = range(self.ngrp)

        if ax == None:
            fig = plt.figure(figsize=(12,10))
            ax = fig.gca()
        else:
            fig = None

        length = self.m - self.init_index
        shift = length + pred

        # Prepare solns with the solves from the sum of selected variables
        solns = np.zeros(self.solns[0].shape) # Same shape as solns
        for index in range(self.solns[0].shape[0]):#enumerate(self.samples):
            tmp = self.solns_plain[ index, :, :]
            solns[index,:] = (tmp @ mask).flatten()


        days = mdates.drange( self.init, self.init+dt.timedelta(shift), dt.timedelta(days=1)) # how often do de plot
        #days_pred = mdates.drange( self.init+dt.timedelta(shift), self.init+dt.timedelta(shift_pred), dt.timedelta(days=7)) # how often do de plot
        self.days = days

        # To save all the data for the plot, len(mexico.days) rows with days
        # columns: year, month, day, datum, datum_pred, map, q_05, q_25, q_50, q_75, q_95
        #             0      1   2     3            4    5     6    7     8     9     10
        sv = -np.ones(( len(days), 11))
        for i,day in enumerate(days):
            d = dt.date.fromordinal(int(day))
            sv[ i, 0] = d.year
            sv[ i, 1] = d.month
            sv[ i, 2] = d.day

        if title==None:
            title = r'$' + v + r'$'
        ylabel = 'Number of cases'
        if log:
            ylabel += ' (log)'


        # Add the map
        #sv[ :shift, 5] = my_soln[self.init_index:(self.init_index+shift)].flatten()

        # plot uncertainty with boxplots
        #yl=ax.get_ylim()
        #xl=ax.get_xlim()
        #ax.boxplot( solns[ :, self.init_index:(self.init_index+length)], positions=days[:length], widths=0.25, sym=' ', whis=[q[0],q[-1]])#, manage_xticks=True)
        #ax.boxplot( solns[ :, (self.init_index+length):(self.init_index+shift):every], positions=days[length:shift:every], widths=0.25, sym=' ', whis=[q[0],q[-1]])#, manage_xticks=True)
        #ax.set_xlim(xl)
        #ax.set_ylim(yl)

        for i in range(length):
            sv[ i, 6:11] = np.quantile( solns[ :, self.init_index+i], q=np.array(q)/100)
        for i in range( length, shift, every):
            sv[ i, 6:11] = np.quantile( solns[ :, self.init_index+i], q=np.array(q)/100)

        self.acme = []
        for d in np.argmax(solns[ :, self.init_index:], axis=1):
            self.acme += [self.init + dt.timedelta(int(d))]

       # plot model prediction at median, with shaped areas for quantiles
        if log:
            ax.semilogy( days[:shift], sv[:shift, 8], '-', linewidth=2, color=color, label=label)
            if blue: #Blue shaowed quantiles
                ax.fill_between( days[:shift], sv[:shift, 6], sv[:shift, 10], color='blue', alpha=0.25)
                ax.fill_between( days[:shift], sv[:shift, 7], sv[:shift, 9], color='blue', alpha=0.25)
            else:
                ax.semilogy( days[:shift], sv[:shift, 6], '--', color=color_q, linewidth=1)
                ax.semilogy( days[:shift], sv[:shift, 10], '--', color=color_q, linewidth=1)
        else:
            ax.plot( days[:shift], sv[:shift, 8], '-', linewidth=2, color=color, label=label)
            if blue: #Blue shaowed quantiles
                ax.fill_between( days[:shift], sv[:shift, 6], sv[:shift, 10], color='blue', alpha=0.25)
                ax.fill_between( days[:shift], sv[:shift, 7], sv[:shift, 9], color='blue', alpha=0.25)
            else:
                ax.plot( days[:shift], sv[:shift, 6], '--', color=color_q, linewidth=1)
                ax.plot( days[:shift], sv[:shift, 10], '--', color=color_q, linewidth=1)

        # format plot
        ax.set_title(self.Region + ' ' + title)
        ax.legend(loc=0, shadow = True)
        # x-axis
        ax.set_xlabel("Date (day.month)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d.%m'))
        """
        day_jumps = np.array([ 1, 7, 7*2, 7*4, 7*8])
        w = np.where((shift / day_jumps) < 9)[0][0]
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=day_jumps[w]))
        if w > 1:
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            #ax.xaxis.set_minor_formatter(mdates.DateFormatter('%d.%m'))
            #ax.tick_params( which='minor', axis='x', labelsize=8)
        """
        #ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        if shift < 190:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        else:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        ax.tick_params( which='major', axis='x', labelsize=12)#, labelrotation=40)
        plt.setp( ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        # y-axis
        ax.set_ylabel(ylabel)
        ax.set_ylim((0, 1.1*np.max(sv[:,-1]))) #Max of upper quantiles
        # Grid
        ax.grid(color='grey', which='both', linestyle='--', linewidth=0.5)

        if right_axis and not(log):
            ax_p = ax.twinx()
            y1, y2 = ax.get_ylim()
            ax_p.set_ylim( y1*1e5/self.N, y2*1e5/self.N)
            ax_p.set_ylabel('per 100,000')

        if csv_fnam != None:
            q_str = ', '.join(["q_%02d" % (qunt,) for qunt in q])
            np.savetxt( csv_fnam, sv, delimiter=', ', fmt='%.1f', header="año, mes, dia, datum, datum_pred, map, " + q_str, comments='')
        return ax

    def StateSlice( self, v, pred, grps='all', E_range = 'all'):
        """:meta private:"""

        mask = np.zeros((self.q*self.ngrp,1))
        for var in v.split():
            mask += self.GetMask( var, grps=grps, E_range=E_range, as_col_vec=True)

        if pred > self.solns[0].shape[1] - self.m:
            print("Maximum prediction days pred=%d" % (self.solns[0].shape[1] - self.m))
            return

        if isinstance( grps, int):
            grps = [grps]
        if isinstance( grps, str): #'all' groups
            grps = range(self.ngrp)

        length = self.m - self.init_index
        shift = length + pred

        if shift < 0:
            print("Negative (in the past)  pred=%d before initial date." % (pred,))
            return

        self.SetTime( pred )
        solns = np.zeros(self.samples.shape[0])
        for index,m in enumerate(self.samples):
            sv = self.solns_plain[index,shift,:]
            solns[index] = (sv @ mask).flatten()
        return solns

    def PlotStateSlice( self, v, pred, grps='all', E_range = 'all', title=None, bins='auto',\
                     ax=None, color='blue', alpha=0.5, prop_axis=False, csv_fnam=None):
        """ Plot posterior of state variable(s) at ``pred``.
            Can be used to state long term evolution, ie when ``pred`` is big.
            If ``v`` represents many state vars then sum their values.

            Args:

                v: state var(s).
                pred: number of days.
                E_range: What Erlang list to include in variable, default ``'all'``.
                ax (optional): axis where to print the plot
                csv_fnam (optional): name of file to save the data for the plot

        """

        if ax == None:
            fig = plt.figure(figsize=(12,10))
            ax = fig.gca()
        else:
            fig = None

        solns = self.StateSlice( v, pred, grps=grps, E_range = E_range)
        length = self.m - self.init_index
        shift = length + pred

        ax.hist( solns, density=True, bins=bins, color=color, alpha=alpha)
        ax.set_title(self.Region)
        ax.set_xlabel(r"$%s$, %s" % (v,(self.init+dt.timedelta(shift)).__str__()))
        ax.set_ylabel("Densidad")
        if prop_axis:
                ax.set_xticklabels(["%4.2f" % (tick/self.N,) for tick in ax.get_xticks()])

    def PlotOmega(self, ax=None, sel=None, mul_f=False,\
                  colors=['red', 'orange', 'gold', 'lawngreen', 'lightseagreen', 'royalblue', 'blueviolet']):
        """Plot the estimated :math:`\omega` population multipliers,
           for ``relax_day != []`` .

           Args:

               sel: selection of :math:`\omega`s, `(c1,c2)`, if None, select all
               mul_f: if ``True``, multiply by ``self.f``

        """
        if self.relax_day == []:
            print("No relax_day, no \omega estimation.")
            return

        if ax == None:
            fig = plt.figure(figsize=(8,6))
            ax = fig.gca()
        else:
            fig = None

        if sel != None:
            ite = self.omega_ite[sel[0]:sel[1]]
        else:
            ite = self.omega_ite

        if mul_f:
            mul = self.f
        else:
            mul = 1

        for c,i in enumerate(ite):
            ax.hist(mul*self.samples[:,i-1], density=True, color='grey', alpha=0.2) #After relax_day
            ax.hist(mul*self.samples[:,i-1], density=True, color=colors[c], histtype='step', linewidth=2) #After relax_day
        if mul_f:
            ax.set_xlabel(r"$f \ \omega$")
        else:
            ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(r"Densidad")
        return ax








