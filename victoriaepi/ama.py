#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import sys

import numpy as np
import scipy.stats as ss
from scipy import integrate
import matplotlib.pyplot as plt

from . import victoria
from . import pytwalk
import pickle



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




def Model_ama( m_list, e_list, f=5.0/100.0, g=0.03, h=0.4, i=0.5, prn=False):
    """
    Define the graph matrix describing the model.

    .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: Model_ama
            :lines: 1,10-

    """

    T = victoria.AuxMatrix(names="S E I^A I^S I^B H^1 H^2 H^3 U^1 U^2 D R", prn=prn)

    T.BaseVar('S')

    T.Exit( 'S', 'E')
    T.SplitExit( 'E', 'I^A', 'I^S', 1-f, prob_symb=['1-f','f'])
    T.Exit( 'I^A', 'R')
    #T.Exit( 'E', 'sIS', f)

    T.SplitExit( 'I^S', 'H^1', 'I^B', g, prob_symb=['g', '1-g'])
    T.Exit( 'I^B', 'R')

    T.SplitExit( 'H^1', 'U^1', 'H^2', h, prob_symb=['h', '1-h'])
    T.Exit( 'H^2', 'R')
    T.Exit( 'U^1', 'U^2')
    T.SplitExit( 'U^2', 'D', 'H^3', i, prob_symb=['i', '1-i'])
    T.Exit( 'H^3', 'R')

    T.NoExit('R')
    T.NoExit('D')
    #T.NoExit('SIS')

    T.End()

    # Split in Erlang series of length m
    T.SplitErlang( e_list, m_list)
    return T


class ama2(victoria.mcmc):
    """
        ODE model of COVID19 pandemics, and parameter inference using
        Bayesian inference (MCMC with the twalk).

        Args:
            N (int): population size
            data_fnam (str): Data file name, workdir/data/data_fnam Make sure to process data into a vertical text array
            where the first column is daily deaths and the second daily confirmed cases.
            out_fnam: MCMC output file name, without .txt, workdir/output/out_fnam + '.txt'
            init_index: day number from data start, where to start the plot
            init: date of init_index
            exit_probs: Exit probabilities of graph
            Pobs_I, Pobs_D: Prob of observation probability of recording an infection/death

            R_rates: Residence rates
            trim: how many data to trim

            intervention_day: None, single value or list in ascending order
            relax_day (list): LIST (always) of relaxation days, default []
    """
    def __init__( self, Region, N, intervention_day, data_fnam, out_fnam,\
                 Pobs_I, Pobs_D,\
                 init_index, init, exit_probs, R_rates, trim, relax_day=[], workdir="./../"):


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
                          exit_probs=exit_probs, R_rates=R_rates)

        if isinstance( Pobs_I, float):
            self.Pobs_I = Pobs_I           # probability of recording an infection
            """ Probability of recording an infection."""
            self.Pobs_D = Pobs_D           # probability of recording a death
            """ Probability of recording a death"""
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
        #

    def PrintConfigParameters( self, f=sys.stdout):
        """Prints all configuration parameters to file f (default sys.stdout)."""
        keys=['Pobs_I', 'Pobs_D', 'f', 'g', 'h', 'i']
        print("\nR_rates:", file=f)
        for key in self.R_rates.keys():
            print("%4s, %7s = %6.4f (%4.1f d), E_list= %2d" %\
                  (key,self.R_rates[key][1],self.R_rates[key][0],1/self.R_rates[key][0],self.R_rates[key][2]), file=f)
        print("", file=f)
        for key in keys:
            print("%s :" % (key,), file=f)
            print(self._dict_[key], file=f)
            print("", file=f)

    def Init_fm_matrix( self, m, intervention_day, exit_probs, R_rates):
        """Init the forward map.

           Args:
               m : is the number of observed days, ie. sample size.
               intervention_day : number from day 0 in data
               exit_probs : f, g, h, i
               R_rates : dictionary of residence rates, eg::

                   R_rates={ 'E'  :[1/1.5, r'\sigma_1'], 'I^S':[1/2  , r'\sigma_2']}

        """
        self.num_pars = 3 + len(intervention_day) # Number of parameters to be inferred

        self.intervention_day = intervention_day

        # Known parameters, set exit probs
        f, g, h, i = exit_probs
        self.f = f           # fraction of severe infections
        """fraction of severe infections """
        self.g = g           # fraction of severe infections that require hospitalization
        """fraction of severe infections that require hospitalization """
        self.h = h           # fraction of hospitations that require ICU
        """fraction of hospitations that require ICU """
        self.i = i           # fraction of ICU patients who die
        """fraction of ICU patients who die """

        self.factor_foi = 0.1*20


        self.R_rates = R_rates
        m_list=[self.R_rates[v][2] for v in self.R_rates.keys()]
        e_list=list(self.R_rates.keys())

        # Define the graph matrix describing the model (see above)
        self.T = Model_ama( m_list, e_list, f=f, g=g, h=h, i=i, prn=False)
        # The rhs will be ('@' matrix multiplication of arrays in Python 3):
        #Graph matrix     State vars     Erlang mask      par (see below)
        #rhs(x)= M      @       x      *       E        @  par
        #       qxq            qx1            qxn          nx1
        # Before this, par is filled with non-linear terms etc. like the force of infection
        # n original number of state variables
        # q final number of state variables after Erlang series


        self.n = self.T.n
        self.q = self.T.q # Total number of state variables
        """Total number of state variables"""
        # Call the base class Init_fm_matrix
        # p=2 size of the return list for solve, incidence and deths
        # quad_k=10 number of subdivions in day for quadrature
        super().Init_fm_matrix( m, num_state_vars=self.q, p=2, quad_k=10)

        # "S E I^A I^S I^B H^1 H^2 H^3 U^1 U^2 D R"
        self.par = np.zeros(self.n) # Par mask

        self.R_rates = R_rates
        # Known parameters, set residence rates:
        for v in R_rates.keys():
            self.par[self.T.ConvertVar(v)] = R_rates[v][0]

        # Auxiliars
        self.sigma_1 = self.par[self.T.ConvertVar('E')] # Used in solve
        self.gamma_1 = self.par[self.T.ConvertVar('I^A')]
        self.sigma_2 = self.par[self.T.ConvertVar('I^S')]

        # The masks to select variables from list of state variables
        self.mask_S   = self.T.SelectMask('S')
        self.mask_E   = self.T.SelectMask('E')
        self.mask_IA  = self.T.SelectMask('I^A')
        self.mask_IAs = self.T.SelectMask('I^A', E_range='all')
        self.mask_ISs = self.T.SelectMask('I^S', E_range='all', as_col_vec=True)
        self.mask_IBs = self.T.SelectMask('I^B', E_range='all', as_col_vec=True)
        self.mask_ISs_flat = self.T.SelectMask('I^S', E_range='all')
        self.mask_IS  = self.T.SelectMask('I^S')
        self.mask_D   = self.T.SelectMask('D')

        #self.psi = 10
        #IA_k = np.where(self.mask_IA == 1)[0][0]
        #self.T.M[0,IA_k] = -self.psi/(1+self.psi)

        self.X0 = np.zeros((self.q,))

    def SetTime( self, shift=0):
        """Extends :meth:`victoriaepi.victoria.mcmc.SetTime`

        .. literalinclude:: ../victoriaepi/ama.py
                :pyobject: ama2.SetTime
                :lines: 1,8-
        """
        super().SetTime(shift=shift)
        self.result_I = np.zeros(self.nn-1) # To hold result of quadrature
        self.result_D = np.zeros(self.nn-1) # To hold result of quadrature
        self.result_H1 = np.zeros(self.nn-1) # To hold result of quadrature

    def GetMask( self, v, E_range='all', as_col_vec=False):
        """Returns a mask to select variable v from grp_list.

            Args:
                v: variable
                E_range: = [0] (default), first in the list, or original variable if no Erlang list.
                E_range: = 'all' use the whole Erlang list for variable or provide E_range list manually.
        """
        return self.T.SelectMask( v, E_range=E_range, as_col_vec=as_col_vec)

    def rhs( self, x, t, p):
        """ See :meth:`victoriaepi.ama.odeint` to check usage.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.rhs
            :lines: 1,9-
        """



        beta1 = p[3 + np.where(t < self.intervention_day)[0][0]]

        I_A = np.sum(x * self.mask_IAs) # total number of asymptomatic infections
        I_S = np.sum(x * self.mask_ISs_flat) # total number of asymptomatic infections
        #S = np.sum(x * self.mask_S) # foi *= (2*S/self.N)

        #
         #force of infection beta1*I^A/N + some factor of  beta1*I^S/N
        foi = (I_A + self.factor_foi*I_S)/self.N * beta1

        self.par[self.T.ConvertVar('S')] = foi

        return self.T.M @ (x * (self.T.par_mask @ self.par))

    def solve_plain1(self, p, quad=True):
        """
        Solve the initial value problem.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.solve_plain1
            :lines: 1,9-
        """
        #beta_1 = p[0] # asymptomatic infection rate
        #beta_2 = p[1] # hospitalized infection rate

        self.X0 *= 0
        self.X0 += p[0]*self.mask_E
        self.X0 += p[1]*self.mask_IA
        self.X0 += p[2]*self.mask_IS
        self.X0 += (self.N-p[0]-p[1]-p[2])*self.mask_S
        if quad:
            return odeint(self.rhs, self.X0, self.t_quad, args=(p,))
        else:
            return odeint(self.rhs, self.X0, self.time, args=(p,))

    def solve_plain2(self, p, quad=True):
        """
        Solve the initial value problem with relaxation day.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.solve_plain2
            :lines: 1,9-
        """
        self.N = p[-self.num_omegas]*self.N_org #\omega_0*N_org, first N_eff

        self.X0 *= 0
        self.X0 += p[0]*self.mask_E
        self.X0 += p[1]*self.mask_IA
        self.X0 += p[2]*self.mask_IS
        self.X0 += (self.N-p[0]-p[1]-p[2])*self.mask_S
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

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.solve
            :lines: 1,10-
        """
        # Use the solver:

        self.soln = self.solve_plain( p )

        # Integrate to produce accumulated number of infected
        for k in range( 0, self.nn-1):
            #x_e = self.soln[ (10*k):(10*(k+1)+1),:] @ self.mask_E
            #incidence = self.f*self.sigma_1*x_e #To hospital
            #self.result_I[k] = self.dt*np.dot(self.weigths,incidence)
            x_s = np.sum( self.soln[ (10*k):(10*(k+1)+1),:] @ self.mask_ISs, axis=1)
            incidence = self.g*self.sigma_2*x_s #To hospital
            self.result_H1[k] = self.dt*np.dot(self.weigths,incidence)
        return self.result_H1, np.diff(self.soln[::10,:] @ self.mask_D)#, self.result_H1

    def llikelihood( self, p):
        """
        Log likelihood.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.llikelihood
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
        theta = 4.0#1.0 #0.01
        r = mu_I/(omega-1.0+theta*mu_I)
        q = 1.0/(omega+theta*mu_I)
        log_likelihood_I = np.sum(ss.nbinom.logpmf( self.data[:,1]+3, r, q))
        # negative binomial likelihood for deaths
        omega = 2.0
        theta = 4.0 #0.5   #antonio 0.5
        r = mu_D/(omega-1.0+theta*mu_D)
        q = 1.0/(omega+theta*mu_D)
        log_likelihood_D = np.sum(ss.nbinom.logpmf( self.data[:,0]+3,r,q))
        # joint likelihood
        log_likelihood = log_likelihood_D + log_likelihood_I

        return log_likelihood

    def lprior( self, p):
        """
        Log prior.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.lprior
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

        """#"""
        if self.relax_day != []:
            for i in self.omega_ite:
                log_prior += ss.beta.logpdf( p[i], 1+1/6, 1+1/3)

        return log_prior

    def support( self, p):
        """
        Support.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.support
            :lines: 1,10-
        """
        rt = True
        rt &= (0.0 < p[0] < 10.0**3)
        rt &= (0.0 < p[1] < 10.0**3)
        rt &= (0.0 < p[2] < 10.0**3)
        # All betas (and omegas, to make it simple!) in [0,20]
        rt &= all((0.0 < p[3:]) * (p[3:] < 20.0))
        if not(rt):
            return False  # Fast track, avoid checking the omegas, if any


        if self.relax_day != []:
            rt &= (0.0 < p[-self.num_omegas] < 1.0) #\omega0
            rt &= ( p[-self.num_omegas]*self.N_org-p[0]-p[1]-p[2] > 0.0)
            for i in self.omega_ite[:-1]:
                rt &= (p[i] < p[i+1] < 1.0) #\omega_i < \omega_{i+1}

        return rt

    def sim_init(self):
        """Simulate initial values for mcmc.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.sim_init
            :lines: 1,8-
        """
        p = np.zeros(self.num_pars)
        p[0] = np.random.uniform(low = 0.01, high = 10.0)
        p[1] = np.random.uniform(low = 0.01, high = 10.0)
        p[2] = np.random.uniform(low = 0.01, high = 10.0)

        p[3:] = np.random.uniform(low = 0.01, high = 5.0, size=self.num_pars-3)

        """#"""
        if self.relax_day != []:
            p[-self.num_omegas] = ss.beta.rvs( 2.0, 4.0)
            for i in self.omega_ite[:-1]:
                p[i+1] = p[i] + (1-p[i])*ss.uniform.rvs()

        return p


    def PlotEvolution( self, pred, ty=0, cumm=False, log=False, ax=None,\
                       csv_fnam=None, q=[ 10, 25, 50, 75, 90], blue=True, add_MRE=False,\
                       color='red', color_q='black', label='Mediana', right_axis=True, label_cases=True):
        """Plot Evolution.

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

    def PlotPanelEvolution( self, pred1, pred2, ty):
        """ Plot Panel Evolution.

        .. literalinclude:: ../victoriaepi/ama.py
            :pyobject: ama2.PlotPanelEvolution
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

