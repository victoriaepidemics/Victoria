# -*- coding: utf-8 -*-




import os
import pickle
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from . import pytwalk










AuxMatrixLatexHead = r"""
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Latex file automatically created by covid_fm.AuxMatrix
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\documentclass[11pt]{article}

%% tickz stuff %%
\usepackage{tikz}
\usetikzlibrary{fit,positioning}
\usetikzlibrary{matrix}
\usetikzlibrary{calc}
\usetikzlibrary{arrows}
\tikzstyle{int}=[draw, fill=white!8, minimum size=2em]
\tikzstyle{init} = [pin edge={to-,thin,black}]
\newcommand\encircle[1]{%
  \tikz[baseline=(X.base)]
    \node (X) [draw, shape=circle, inner sep=0.5] {\strut #1};}

%%%<
\usepackage{verbatim}
\usepackage[active,tightpage]{preview}
\PreviewEnvironment{tikzpicture}
\setlength\PreviewBorder{5pt}%
%%%>

\begin{document}
\pagestyle{empty}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{tikzpicture}[node distance=2.5cm,auto,>=latex']

"""

AuxMatrixLatexTail = r"""
\end{tikzpicture}

\end{document}
"""

#import os

class AuxMatrix:
    """Auxiliary class to define the ODE transition matrix to build models.


    Args:
        names: defaults to None
        num_state_vars:  defaults to None
        prn: to print or not the matrix after each operation, defaults to False
        tex_fnam: defaults to None
        pdflatex: defaults to "/Library/TeX/texbin/pdflatex"


    """
    def __init__( self, names=None, num_state_vars=None, prn=False, tex_fnam=None, pdflatex="/Library/TeX/texbin/pdflatex"):

        if names == None:
            self.names = [r"V^{%d}" % (i,) for i in range(num_state_vars)]
        else:
            self.names = names.split()
        self.names_org = self.names.copy()

        if num_state_vars == None:
            self.q = len(self.names)
        else:
            self.q = num_state_vars
        self.n = self.q
        self.Erlang_len = [1]*self.n # Length of Erlang list for each variable, we start with 1 (no Erland list in fact)

        self.M = np.diagflat( -np.ones(self.q), k=0 ) # Initial matrix, every var has an exit (-1 in the diagonal)
        self.par_mask = np.diagflat( np.ones(self.q), k=0 )
        self.par_vec = np.zeros((self.q,1))
        self.prn = prn # Print info commments

        # To produce the latex diagram of the graph
        if tex_fnam == None:
            self.tex_f = open( 'AuxMatrixGraph.tex', 'w')
            self.tex_fnam = 'AuxMatrixGraph.tex'
        else:
            self.tex_f = open( tex_fnam, 'w')
            self.tex_fnam = tex_fnam
        self.pdflatex = pdflatex

    def ConvertVar( self, v):
        """If ``v`` is an `int` return it as is, if it is a string, return its index
           in the list of variable names.
        """
        if isinstance( v, int):
            return v
        else:
            return self.names_org.index(v)

    def PrintMatrix( self, j1=None, j2=None):
        """Prints the matrix and sets the default to ``prn=True``."""
        if j1 == None:
            j1 = 0
        if j2 == None:
            j2 = self.q

        print("\n%9s" % (" ",), end=' ')
        for j in range( j1, j2):
            print("%9s" % (self.names[j],), end=' ')
        print("")
        for i in range( j1, j2):
            print("%9s" % (self.names[i],), end=' ')
            for j in range( j1, j2):
                print("%9.2f" % (self.M[i,j],), end=' ')
            print("")

    def BaseVar( self, v):
        """Define a base of source var with no entry."""
        # Noting in the ode in the matrix, this is part of the ODE initial conditions
        # Only plot the node
        v = self.ConvertVar(v)
        print(AuxMatrixLatexHead, file=self.tex_f)
        print(r"\node [int] (%s) {$%s$};" % (self.names[v],self.names[v]), file=self.tex_f)
        print("", file=self.tex_f)

    def SplitExit( self, v, v1, v2, prob, prob_symb=None):
        """Split the exit of the variable ``v`` (row ``v+1``) to variables ``v1`` and ``v2``,
           with the proportions ``prob`` and ``1-prob``, respectively.
           """
        v = self.ConvertVar(v)
        v1 = self.ConvertVar(v1)
        v2 = self.ConvertVar(v2)
        self.M[v+1, v] = 0 # Delete its exit
        self.M[ v1, v] = prob # Split
        self.M[ v2, v] = 1-prob

        if  prob_symb == None:
            p = "%4.2f" % (prob,)
            p1 = "%4.2f" % (1-prob,)
        else:
            p = prob_symb[0]
            p1 = prob_symb[1]
        print(r"\node [int] (%s) [above right of = %s] {$%s$};" % (self.names[v1],self.names[v],self.names[v1]), file=self.tex_f)
        print(r"\node [int] (%s) [below right of = %s] {$%s$};" % (self.names[v2],self.names[v],self.names[v2]), file=self.tex_f)
        print(r"\path[->] (%s) edge node  {$%s$} (%s);" % (self.names[v],p,self.names[v1]), file=self.tex_f)
        print(r"\path[->] (%s) edge node  {$%s$} (%s);" % (self.names[v],p1,self.names[v2]), file=self.tex_f)
        print("", file=self.tex_f)

        if self.prn:
            print("%s --> split --> %s (%4.2f) and --> %s (%4.2f)" %\
                  ( self.names[v], self.names[v1], prob, self.names[v2], 1-prob))

    def Exit( self, v, v1, w=1):
        """Put an exit to variable ``v1`` from variable ``v``, with weight ``w`` (positive)."""
        v = self.ConvertVar(v)
        v1 = self.ConvertVar(v1)
        self.M[ v1, v] = w
        print(r"\node [int] (%s) [right of = %s] {$%s$};" % (self.names[v1],self.names[v],self.names[v1]), file=self.tex_f)
        print(r"\path[->] (%s) edge node {} (%s);" % (self.names[v],self.names[v1]), file=self.tex_f)
        print("", file=self.tex_f)

        if self.prn:
            print("%s --> %s" % (self.names[v],self.names[v1]))

    def NoExit( self, v):
        """Variable ``v`` has no exit."""
        v = self.ConvertVar(v)
        self.M[ v, v]  = 0

    def End( self):
        """Does nothing to the matrix and flushes the latex output."""
        print(AuxMatrixLatexTail, file=self.tex_f)
        print("", file=self.tex_f)
        self.tex_f.close()
        """
        try:
            os.system(self.pdflatex + " " + self.tex_fnam + "> AuxMatrixGraph_latex_log")
        except:
            print("Could not create graph png plot, perhaps /Library/TeX/texbin/pdflatex not available.")
        """

    def SplitErlang( self, v_list, m_list):
        """Split the exit of the variable list ``v_list`` into a chain of m variables, where
           m is taken from ``m_list``. If ``m_list = m``,and ``int``, use the same ``m``
           for all variables.

           Args:

               v_list: variable list
               m_list: m list. If m_list is ``int``, use the same value for all variables


        """
        if isinstance( m_list, int):
            m_list = [m_list]*len(v_list) # Same m for all
        # Sort in ascending order the variables
        v_list = [self.ConvertVar(v) for v in v_list]
        tmp_L = [ (v_list[i],i) for i in range(len(v_list)) ]
        tmp_L.sort()
        v_list, permutation = zip(*tmp_L) # Sorted list of vars and their corresponding
        m_list = [m_list[permutation[i]] for i in range(len(m_list))] # Erlang list lengths
        # Add the Erlang lists
        E_sum = 0
        for k,v_org in enumerate(v_list): # v_org  original variable number
            m = m_list[k]
            if m == 1:
                continue #Nothing to do
            elif m < 1:
                raise NameError('covid_fm: Negative Erlang list')
            v = v_org + E_sum #k*(m-1)
            nm = self.names[v]
            l = 0
            self.names[v] += r"_{%d}" % (l,)
            self.par_mask[v, v_org] = m
            for i in range( v+1, v+m):
                self.M = np.insert( self.M, i, 0, axis=1)
                self.M = np.insert( self.M, i, 0, axis=0)
                self.q += 1
                self.M[ i:,  i] = self.M[ i:, i-1] # Move the column below, to move all exists
                self.M[ i:,i-1] *= 0 # Change all previous exists to zero
                self.M[ i, i-1] = 1 # Add an entry to new variable
                self.M[ i,   i] = -1 # Add its exit, to the previous exit of variable v

                l += 1
                self.names.insert( i, nm + r"_{%d}" % (l,))
                self.par_mask = np.insert( self.par_mask, i, 0, axis=0)
                self.par_mask[ i, v_org] = m
            self.Erlang_len[v_org] = m
            E_sum += m-1 # Add m-1 state variables
        return True

    def GetE_range( self, v):
        """ Returns self.Erlang_len[v]."""
        v = self.ConvertVar(v)
        return self.Erlang_len[v]

    def SelectMask( self, v, E_range=[0], as_col_vec=False):
        """Create a selection mask with 1's in variable ``v`` and, maybe, within
           its list of the 'Erlang series' of variables ``E_range``.
           ``E_range = [0]`` first in the list as default, or with the original variable if there is no Erlang list.
           ``E_range = 'all'`` use the whole Erlang list
           or provide ``E_range`` list manually.

           Args:
               v: state variable
               E_range: list of 'Erlang series'


        """
        v = self.ConvertVar(v)

        if isinstance( E_range, str):
            if E_range == 'all':
                E_range = range(self.Erlang_len[v])
            else:
                raise NameError('covid_fm: Unknown option E_range=%s' % (E_range,))
        if max(E_range) > self.Erlang_len[v]:
            raise NameError('covid_fm: var in Erlang list beyond list length %d:' % (self.Erlang_len[v],), E_range)

        # Make a vector with a 1 at row v in a column vector of zeros with self.n rows
        tmp = np.zeros((self.n,1))
        tmp[v,0] = 1
        # with the mask we obtain all variables in the Erlang list (possibly length 1)
        tmp2 = self.par_mask @ tmp
        # Where is the first variable in the list:
        i0 = np.where(tmp2 != 0)[0][0]
        # Slect the variables in the Erlang list
        tmp2 *= 0
        for i in E_range:
            tmp2[i0+i] = 1
        if as_col_vec:
            return tmp2
        else:
            return tmp2.flatten()



class fm:
    """
        A template class, in which one needs to fill the methods: ``rhs``, ``solve_plain``, ``solve``, in order to set the dynamics of the boxes in the forward map.

        Parameters
        ----------
        m : int
             is the number of observed days, ie. sample size

        num_state_vars : int
            is the number of state variables
        p : int
            length of solve list (see solve below)
        quad_k : int
            quad_k=10 time of subdivisions for quadrature



        """
    def __init__( self, m, num_state_vars, p, quad_k):
        self.m = m
        """the number of observed days, ie. sample size"""
        self.num_state_vars = num_state_vars

        #Auxiliary arrays

        self.X0 = np.zeros(self.num_state_vars)  # Initial conditions
        """Initial conditions."""

        # For quadrature:
        self.weigths = np.ones(11)   # quadrature weights
        """quadrature weights."""
        self.weigths[0] = 0.5
        self.weigths[-1] = 0.5
        self.dt = 1.0/(10.0) # dt for quadrature

        self.quad_k = quad_k # subdivisions in between days, for quadrature
        self.p = p
        self.SetTime() # Set the time ranges for quadrature and the fm

    def SetTime( self, shift=0):
        """Set the time range from day 0 to m+shift days inclusive."""
        self.time = np.linspace( -1.0, self.m-1+shift, num=(self.m+shift)+1, endpoint=True) # observation time
        self.nn = len(self.time)
        self.t_quad = np.linspace( -1.0, self.m-1+shift, num=self.quad_k*(self.m+shift)+1, endpoint=True) # Grid for quadrature
        self.n_quad = len(self.t_quad)

    def rhs( self, x, t, p):
        """Template rhs of the ODE, typically using a matrix created with AuxMatrix, with parameters p.

           The rhs will be ('@' matrix multiplication of arrays in Python 3)

           .. code-block:: text

                Graph matrix     State vars     Erlang mask      par (see below)
                rhs(x)= M      @       x      *       E        @  par
                        qxq            qx1            qxn          nx1


        Before this, ``par`` is filled with non-linear terms like the force of infection,
        ``n``, the original number of the state variables ,and
        ``q``, the final number of the state variables after Erlang series.

        """
        pass

    def solve_plain(self, p, quad=True):
        """Template. Solve the system for parameters ``p``."""
        pass

    def solve( self, p):
        """Template. Solve the forward map with the obs functional ready for the likelihood.
        """
        # template class
        # Returns a *list* of solves of size self.k
        pass

    """ %%%%% end of the former fm_matrix class %%%%%%%%% """







class mcmc(fm):
    """
         Victoria epidemics ODE model, where Bayesian inference is implemented to get the initial conditions (MCMC with the twalk).

        Args:

            Region (str): Region name
            N (int): population size
            data_fnam (Str): Data file name, workdir + 'data/' + data_fnam make sure to process data into a vertical text array
            out_fnam  (Str): MCMC output file name, without .txt, e.g. workdir/output/out_fnam + '.txt'
            init_index  (int) : day number from data start, where to start the plot
            init (date): date of init_index
            trim (int): how many data to trim to the past, NEGATIVE NUMBER

    """
    def __init__( self, Region, N, data_fnam, out_fnam,\
                 init_index, init, trim, workdir="./../"):

        self.Region = Region
        self.N = N        # population size

        # Load number of new cases per day below.

        #data_fnam  # Data file name, workdir/data/data_fnam
        data = np.loadtxt(workdir + 'data/' + data_fnam)
        self.init_index = init_index
        self.init = init # # Initial position for plotting
        self.trim = trim
        if self.trim < 0:
            if len(data.shape) == 1: #one data column read only
                self.data = data[:self.trim]
                self.data_trimed = data[self.trim:]
            else: # several data columns read
                self.data = data[:self.trim,:]
                self.data_trimed = data[self.trim:,:]
        else:
            self.data = data
            if len(data.shape) == 1: #one data column read only
                self.data_trimed = 20+np.zeros(2) # Dummy
            else:
                self.data_trimed = 20+np.zeros((2,data.shape[1])) # Dummy

        self.out_fnam = out_fnam
        if os.path.isfile(workdir + 'output/' + self.out_fnam + '_solns.pkl'): # samples file exists
            print("File with mcmc samples exists, loading samples ...", end=' ')
            self.samples = pickle.load(open(workdir + 'output/' + self.out_fnam + '_samples.pkl', 'rb'))
            self.solns = pickle.load(open(workdir + 'output/' + self.out_fnam + '_solns.pkl', 'rb'))
            self.solns_plain = pickle.load(open(workdir + 'output/' + self.out_fnam + '_solns_plain.pkl', 'rb'))
            self.essize = self.samples.shape[0]
            print("effective sample size: %d" % (self.essize,))
        else:
            print("File with mcmc samples does not exist, run RunMCMC first.")
        self.workdir = workdir
        self.mcmc_samples = False

        # Typically Init_fm_matrix should be called here as well
        #self.Init_fm_matrix( m=self.data.shape[0], num_state_vars=self.q)



    """
    This used to be the class fm_matrix
                covid_fm_matrix.py

    Marcos Capistran, Antonio capella y Andres Christen,
    Marzo-Octubre 2020.


    # Derived forward map class from template class fm
    class fm_matrix(fm):
        Init_fm_matrix used to be __init__
    """
    def Init_fm_matrix( self, m, num_state_vars, p, quad_k):
        """
        Init the forward map.
        Template class, here the forward map (model) is created and initialized.

        Args:

            m (int): is the number of observed days, ie. sample size
            num_state_vars (int): is the number of state variables

            p: =1 Need to initialize the size of the fm return list, e.g.: 1 if only infectious are observed, 2 if I and D are observed, etc

            quad_k: number of subdivisions of day for possible quadrature in solve

        """
        super().__init__( m=m, num_state_vars=num_state_vars, p=p, quad_k=quad_k)



    def llikelihood( self, p):
        """Template. Calculates the log-likelihood of p."""
        # template class
        return 1.0

    def lprior( self, p):
        """Template. Calculates the log-prior of p."""
        # template class
        return 1.0

    def energy( self, p):
        """
        Returns minus log-post of p::

            return -1*(self.llikelihood(p) + self.lprior(p))

        """
        return -1*(self.llikelihood(p) + self.lprior(p))

    def support( self, p):
        """Template. Returns True if p is in the parameter support, False otherwise."""
        # template class
        return True

    def sim_init(self):
        """Template. Simulate initial values for mcmc."""
        pass

    def RunMCMC( self, T, burnin=1000, pred=100, plot_fit=True):
        """Run twalk MCMC.

        Args:

            T: number of iterations.
            burnin: burn in iterations

        """

        self.SetTime( 0 )
        self.twalk = pytwalk.pytwalk(n = self.num_pars, U=self.energy, Supp =self.support)
        self.twalk.Run( T=T, x0 = self.sim_init(), xp0 = self.sim_init())

        self.mcmc_samples = True

        self.iat = int(self.twalk.IAT(start=burnin)[0,0])
        self.burnin = burnin
        print("\nEffective sample size: %d" % ((T-burnin)/self.iat,))
        self.samples = self.twalk.Output[burnin::(self.iat),:] # Burn in and thining
        self.essize = self.samples.shape[0]

        self.SetTime( pred )

        self.solns = [np.zeros(( self.essize, self.m + pred)) for i in range(self.p)]
        self.solns_plain = np.zeros(( self.essize, self.m + pred, self.q))
        print("Sampling %d model solutions." % ( self.essize,))
        for index,m in enumerate(self.samples):
            tmp = list(self.solve(m[:-1]))
            self.solns_plain[index,:,:] = self.soln[self.quad_k::self.quad_k,:]
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

    def CalculateMAP(self):
        """
        :meta private:
        """
        nmap = np.int(np.where(self.samples[:,-1]==self.samples[:,-1].min())[0][0])
        self.pmap = self.samples[nmap,:]


    def PlotEvolution_fm( self, solns, prevalence, pred, ylabel, log=False, ax=None,\
                       csv_fnam=None, q=[ 10, 25, 50, 75, 90], blue=True, add_MRE=False,\
                       color='red', color_q='black', label='Median', right_axis=True, label_cases=True):
        """Plot Evolution.

        Args:

            pred: number of days to predict
            ty: 0 = Infected,1 = deaths
            cumm: True if cumulative, default False
            log: True if y log scale, default False
            ax: axis where to print the plot (optional)
                csv_fnam: name of file to save the data for the plot (optional)

        """
        # Old, depreciated
        every=1

        if ax == None:
            fig = plt.figure(figsize=(12,10))
            ax = fig.gca()
        else:
            fig = None

        if pred > solns.shape[1] - self.m:
            print("Maximum prediction days pred=%d" % (self.solns[0].shape[1] - self.m))
            return
        pred = max( pred, -self.trim) #Pred should cover the trimmed data

        # length and shift for plotting
        length = self.m - self.init_index
        shift = length + pred
        self.SetTime( pred )

        # The time frame
        days = mdates.drange( self.init, self.init+dt.timedelta(shift), dt.timedelta(days=1)) # how often do de plot
        #days_pred = mdates.drange( self.init+dt.timedelta(shift), self.init+dt.timedelta(shift_pred), dt.timedelta(days=7)) # how often do de plot
        #If some method other requires this

        #Solution to a reported bug with mdates.range
        j=0
        for x in range( self.init.toordinal(),  (self.init+dt.timedelta(shift)).toordinal(), 1):
            days[j] = x
            j = j+1
        self.days = days
        self.shift = shift

        # To save all the data for the plot, len(mexico.days) rows with days
        # columns: year, month, day, datum, datum_pred, map, q_05, q_25, q_50, q_75, q_95
        #             0      1   2     3            4    5     6    7     8     9     10
        sv = -np.ones(( len(days), 11))
        for i,day in enumerate(days):
            d = dt.date.fromordinal(int(day))
            sv[ i, 0] = d.year
            sv[ i, 1] = d.month
            sv[ i, 2] = d.day
        # Save data and predicted data
        sv[:length, 3] = prevalence[self.init_index:]
        if self.trim < 0:
            sv[ length:(length-self.trim), 4] = self.future

        self.MRE = np.zeros(length)
        # Calculate quantiles
        for i in range(length):
            sv[ i, 6:11] = np.quantile( solns[ :, self.init_index+i], q=np.array(q)/100)
            self.MRE[i] = np.mean(np.abs(solns[ :, self.init_index+i]-prevalence[self.init_index+i])/(1+prevalence[self.init_index+i]))
        for i in range( length, shift, every):
            sv[ i, 6:11] = np.quantile( solns[ :, self.init_index+i], q=np.array(q)/100)
        self.PE_solns = solns

        if add_MRE:
            MRE_q = np.quantile( self.MRE, q=[ 0.05, 0.5, 0.95])
            ax.annotate( "MRE: [%5.2f-%5.2f-%5.2f]-%5.2f*" % ( MRE_q[0],MRE_q[1],MRE_q[2], np.max(self.MRE)),\
                    (0,1), (10, -10), xycoords='axes fraction',\
                    textcoords='offset points', va='top', fontsize=10)

        self.acme = []
        for d in np.argmax(solns[ :, self.init_index:], axis=1):
            self.acme += [self.init + dt.timedelta(int(d))]
        self.acme_qs = []
        for a_q in np.quantile(np.argmax(solns[ :, self.init_index:], axis=1), q=np.array(q)/100):
            self.acme_qs += [self.init + dt.timedelta(int(a_q))]

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

        # Plot data and prediction
        if log:
            # #plot data
            if label_cases:
                ax.semilogy( days[:length], prevalence[self.init_index:], 'k*', lw=1, label='Cases')
            else:
                ax.semilogy( days[:length], prevalence[self.init_index:], 'k*', lw=1)
            if self.trim < 0:
                if label_cases:
                    ax.semilogy( days[length:(length-self.trim)], self.future, 'r*', lw=1, label='Future cases')
                else:
                    ax.semilogy( days[length:(length-self.trim)], self.future, 'r*', lw=1)
            ylabel += ' (log)'
        else:
            # #plot data
            if label_cases:
                ax.bar( days[:length], prevalence[self.init_index:], color='blue', width=0.5, alpha=0.5)
                ax.plot( days[:length], prevalence[self.init_index:], 'bo', markersize=2, label='Cases')
            else:
                ax.bar( days[:length], prevalence[self.init_index:], color='blue', markersize=2, width=0.5, alpha=0.5)
                ax.plot( days[:length], prevalence[self.init_index:], 'bo')
            if self.trim < 0:
                if label_cases:
                    ax.bar( days[length:(length-self.trim)], self.future, color='grey', width=0.5, alpha=0.1)
                    ax.plot( days[length:(length-self.trim)], self.future, 'k*', markersize=2, label='Future cases')
                else:
                    ax.bar( days[length:(length-self.trim)], self.future, color='grey', width=0.5, alpha=0.1, lw=1)
                    ax.plot( days[length:(length-self.trim)], self.future, 'k*', markersize=2)

        # format plot
        ax.legend(loc=0, shadow = True)
        # x-axis
        ax.set_xlabel("Date (day.month)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d.%m'))
        #ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        if shift < 190:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        else:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
        ax.tick_params( which='major', axis='x', labelsize=12)#, labelrotation=40)
        plt.setp( ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        # y-axis
        ax.set_ylabel(ylabel)
        ax.set_ylim((0, max( 1.1*np.max(sv[:,-1]), 1.1*np.max(prevalence)) ) ) #Max of upper quantiles
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

    def PlotAcmeHist(self, ax=None, **kwargs):
        """ Plots a histogram of self.acme,
            last produced by :meth:`PlotEvolution_fm` or :meth:`PlotStateVar`.
            ``kwargs`` are passed to ``hist``.

        """
        if ax == None:
            fig = plt.figure(figsize=(5,5))
            ax = fig.gca()
        else:
            fig = None

        ax.hist( self.acme, density=True, **kwargs)
        ax.set_xlabel("Date (day.month)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d.%m'))
        if (max(self.acme) - min(self.acme)) < 15:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        else:
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        ax.tick_params( which='major', axis='x', labelsize=12)
        plt.setp( ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        # y-axis
        ax.set_ylabel("Density")
        return ax

    def PlotFit( self ):
        """Template. Plots basic results of the MCMC."""
        if not(self.mcmc_samples):
            print("MCMC samples not available.")
            return
        pass #Nothing template method


    def PlotStateVar( self, v, pred, E_range = 'all', log=False, title=None,\
                     ax=None, q=[ 10, 25, 50, 75, 90], blue=True,\
                     color='red', color_q='black', label=' ', right_axis=True, csv_fnam=None):
        """ Plot evolution of state var(s)

            Args:
                v: state vars. If it is more than one, then sum their values
                pred: number of days to predict
                E_range: Erlang list to include in variable, as default 'all'
                log: True if y log scale, as default False
                ax: axis where to print the plot (optional)
                csv_fnam: name of file to save the data for the plot (optional)

            Example::

                PlotStateVar('H^2', pred=20), PlotStateVar('H^2 H^3', pred=60)


        """

        # Depreciated
        every=1

        if pred > self.solns[0].shape[1] - self.m:
            print("Maximum prediction days pred=%d" % (self.solns[0].shape[1] - self.m))
            return

        mask = np.zeros((self.q,1))
        for var in v.split():
            mask += self.GetMask( var, E_range=E_range, as_col_vec=True)

        if ax == None:
            fig = plt.figure(figsize=(12,10))
            ax = fig.gca()
        else:
            fig = None


        # length and shift for plotting
        length = self.m - self.init_index
        shift = length + pred

        # Prepare solns with the solves from the sum of selected variables
        solns = np.zeros(self.solns[0].shape) # Same shape as solns
        for index in range(self.solns[0].shape[0]):#enumerate(self.samples):
            tmp = self.solns_plain[ index, :, :]
            solns[index,:] = (tmp @ mask).flatten()

        # the time frame for plotting
        days = mdates.drange( self.init, self.init+dt.timedelta(shift), dt.timedelta(days=1)) # how often do de plot
        #days_pred = mdates.drange( self.init+dt.timedelta(shift), self.init+dt.timedelta(shift_pred), dt.timedelta(days=7)) # how often do de plot

        #Solution to a reported bug with mdates.range
        j=0
        for x in range( self.init.toordinal(),  (self.init+dt.timedelta(shift)).toordinal(), 1):
            days[j] = x
            j = j+1
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

    def StateSlice( self, v, pred, E_range = 'all'):
        """:meta private:"""

        mask = np.zeros((self.q,1))
        for var in v.split():
            mask += self.GetMask( var, E_range=E_range, as_col_vec=True)

        if pred > self.solns[0].shape[1] - self.m:
            print("Maximum prediction days pred=%d" % (self.solns[0].shape[1] - self.m))
            return

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

    def PlotStateSlice( self, v, pred, E_range = 'all', title=None, bins='auto',\
                     ax=None, color='blue', alpha=0.5, prop_axis=False, csv_fnam=None):
        """ Plot posterior of state variable(s) at ``pred``.
            Can be used to state long term evolution, ie when ``pred`` is big.
            If ``v`` represents many state vars then sum their values.

            Args:

                v: state var(s)
                pred: number of days
                E_range: What Erlang list to include in variable, default ``'all'``
                ax (optional): axis where to print the plot
                csv_fnam (optional): name of file to save the data for the plot

        """

        if ax == None:
            fig = plt.figure(figsize=(12,10))
            ax = fig.gca()
        else:
            fig = None

        solns = self.StateSlice( v, pred, E_range = E_range)
        length = self.m - self.init_index
        shift = length + pred

        ax.hist( solns, density=True, bins=bins, color=color, alpha=alpha)
        ax.set_title(self.Region)
        ax.set_xlabel(r"$%s$, %s" % (v,(self.init+dt.timedelta(shift)).__str__()))
        ax.set_ylabel("Density")
        if prop_axis:
                ax.set_xticklabels(["%4.2f" % (tick/self.N,) for tick in ax.get_xticks()])
