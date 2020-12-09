"""
This is the Python implementation of t-walk by Andrés Christen.

The t-walk is a "A General Purpose Sampling Algorithm for Continuous Distributions" to sample from many objective functions (especially suited for posterior distributions using non-standard models that would make the use of common algorithms and software difficult); it is an MCMC that does not required tuning.

See:  http://www.cimat.mx/~jac/twalk/

Related paper is available online at: http://projecteuclid.org/euclid.ba/1340218339

See also: http://www.cimat.mx/~jac/twalk/pytwalktutorial.py
"""



from numpy.random import uniform, normal
from numpy import ones, zeros, cumsum, shape, mat, cov, mean, ceil, matrix, sqrt
from numpy import floor, exp, log, sum, pi, savetxt, loadtxt, array

from time import time, localtime, strftime

from numpy import random


try:
    from matplotlib.pyplot import plot, hist, xlabel, ylabel, title
except:
    print("pytwalk: WARNING: matplotlib.pyplot module not available, Ana, TS and Hist methods will fail.")

# Some auxiliar functions and constants

def SqrNorm(x):
    """Square of the norm.

    :param x: return sum(x*x)
    :return: return sum(x*x)

    """

    return sum(x*x)

log2pi = log(2*pi)
log3 = log(3.0)

def Remain( Tr, it, sec1, sec2):
    """
    Returns messages about the remaining time to finish.

    Parameters
    ----------
    Tr : int
        total iterations
    it : int
        current iteration
    sec1 : float
        start time
    sec2 : float
        current time, as returned by time() (floats)

    Returns
    -------
    float
        Remaining time

    """

    # how many seconds remaining
    ax = int( (Tr - it) *  ((sec2 - sec1)/it) )


    if (ax < 1):

        return " "

    if (ax < 60):

        return "Finish in approx. %d sec." % (ax,)

    if (ax <= 360):

        return "Finish in approx. %d min and %d sec." % ( ax // 60, ax % 60)

    if (ax > 360):

        ax += sec2  # current time plus seconds remaining=end time
        return "Finish by " + strftime("%a, %d %b %Y, %H:%M.", localtime(ax))






class pytwalk:
    """T-walk class. A generic self adjusting MCMC.

        Initiates defining the dimension n, Supp= defines the support,
        returning True if x is within the support and False otherwise,
        and U= -log of the objective function, eg::

            >>> Mytwalk = pytwalk( n=3, U=MyMinusLogf, Supp=MySupportFunction)

        or

        t positive, u= -log likelihood and w= -log prior::

            >>> Mytwalk = pytwalk( n=3, t=0.5, u=MyMinusLoglikelihood, w=MyMinusLogPrior, Supp=MySupportFunction)

        In this case the objective function is U= t*u + w and u, for x (not xp)
        is saved in self.Output_u.  This is a backward-compatible implementation for the penalized likelihood for the thermodynamic integral to estimate normalizing constants.

        Then use the :meth:`pytwalk.Run` method with default values as in the paper, which usually  do not need to be changed.

        Args:
            n: dimension
            ww: the prob. of choosing each kernel, aw, at, n1phi (see inside twalk.py)
            aw: For the walk move.
            at: For the Traverse move.


        """

    def __init__( self, n, U=(lambda x: sum(0.5*x**2)), Supp=(lambda x: True),
        t=-1, u=(lambda x: sum(0.5*x**2)), w=(lambda x: 0.0),
        ww=[0.0000, 0.4918, 0.4918, 0.0082, 0.0082], aw=1.5, at=6.0, n1phi=4.0,
        silent=False):
        # Careful the Hop move does not work!!
        self.n = n
        self.t = t
        if self.t >= 0: # Penalized likelihood
            self.LikelihoodEnergy = u
            self.PriorEnergy = w
            self.Output_u = array([0.0])
        else:  # Usual case
            self.PriorEnergy = (lambda x: 0.0)
            self.LikelihoodEnergy = U
            self.t = 1.0
        self.U = (lambda x: self.Energy(x))
        self.Supp = Supp
        self.Output = zeros((1, n+1)) # No data (MCMC output) yet
        self.Output_u = array([0.0]) # To save ll_e, the likelihood energy
        self.T = 1
        self.Acc = zeros(6)  # To save the acceptance rates of each kernel, and the global acc. rate

        # Kernel probabilities
        self.Fw = cumsum(ww)
        """Kernel probabilities"""

        # Parameters for the propolsals
        self.aw = aw  # For the walk move
        """For the walk move"""
        self.at = at # For the Traverse move
        """For the Traverse move"""

        #n1phi = 5 # expected value of parameters to move
        self.pphi = min( n, n1phi)/(1.0*n) # Prob. of choosing each par.

        self.WAIT = 30

        self.silent=silent

    def Energy( self, x):
        self.ll_e = self.LikelihoodEnergy(x)
        self.prior_e = self.PriorEnergy(x)
        return self.t*self.ll_e + self.prior_e

    def _SetUpInitialValues( self, x0, xp0):
        """Private method."""

        # Check x0 and xp0 in the support

        if any(abs(x0 -xp0) <= 0):
            print("pytwalk: ERROR, not all entries of initial values different.")
            return [ False, 0.0, 0.0]

        if not(self.Supp(x0)):
            print("pytwalk: ERROR, initial point x0 out of support.")
            return [ False, 0.0, 0.0]
        u = self.U(x0)

        if not(self.Supp(xp0)):
            print("pytwalk: ERROR, initial point xp0 out of support.")
            return [ False, u, 0.0]
        up = self.U(xp0)

        return [ True, u, up]



    def Run( self, T, x0, xp0, t=1):
        """Run the twalk.

           Args:

               T:  Number of iterations.
               xp0,x0: two initial points within the support, *each entry of x0 and xp0 must be different*.
        """

        self.t = t

        if not(self.silent):
            sec = time()
            print("pytwalk: Running the twalk with %d iterations"\
                    % (T,), end=' ')
            if self.t == 1:
                print(". ",  strftime("%a, %d %b %Y, %H:%M:%S.", localtime(sec)))
            else:
                print(" (%f). " % (self.t,), strftime("%a, %d %b %Y, %H:%M:%S.", localtime(sec)))

        # Check x0 and xp0 are in the support
        [ rt, u, up] = self._SetUpInitialValues( x0, xp0)

        if (not(rt)):
            return 0


        # send an estimation for the duration of the sampling if
        # evaluating the ob. func. twice (in self._SetUpInitialValues) takes more than one second

        if not(self.silent):
            sec2 = time() # last time we sent a message
            print("       " + Remain( T, 2, sec, sec2))

        x = x0     # Use x and xp by reference, so we can retrive the last values used
        xp = xp0

        # Set the array to place the iterations and the U's ... we donot save up's
        self.Output = zeros((T+1, self.n+1))
        self.Output_u = zeros(T+1)
        self.T = T+1
        self.Acc = zeros(6)
        kercall = zeros(6) # Times each kernel is called


        self.Output[ 0, 0:self.n] = x.copy()
        self.Output[ 0, self.n] = u
        self.Output_u[0] = self.ll_e

        j1=1
        j=0

        # Sampling
        for it in range(T):

            y, yp, ke, A, u_prop, up_prop = self.onemove( x, u, xp, up)

            kercall[ke] += 1
            kercall[5] += 1
            if (uniform() < A):
                x = y.copy()   # Accept the propolsal y
                u = u_prop
                xp = yp.copy()   # Accept the propolsal yp
                up = up_prop

                self.Acc[ke] += 1
                self.Acc[5] += 1


            # To retrive the current values
            self.x = x
            self.xp = xp
            self.u = u
            self.up = up

            self.Output[it+1,0:self.n] = x.copy()
            self.Output[it+1,self.n] = u
            self.Output_u[it+1] = self.ll_e

            # Estimate the remaing time, every 2**j1 iterations
            if not(self.silent):
                if ((it % (1 << j1)) == 0):

                    j1 += 1
                    j1 = min( j1, 10)  # check the time at least every 2^10=1024 iterations
                    ax = time()
                    if ((ax - sec2) > (1 << j)*self.WAIT): # Print an estimation every WAIT*2**j
                        print("pytwalk: %10d iterations so far. " % (it,) + Remain( T, it, sec, ax))
                        sec2 = ax
                        j += 1
                        j1 -= 1 # check the time as often

        if not(self.silent):
            if (self.Acc[5] == 0):
                print("pytwalk: WARNING,  all propolsals were rejected!")
                print(strftime("%a, %d %b %Y, %H:%M:%S.", localtime(time())))
                return 0
            else:
                print("pytwalk: finished, " + strftime("%a, %d %b %Y, %H:%M:%S.", localtime(time())))

        for i in range(6):
            if kercall[i] != 0:
                self.Acc[i] /= kercall[i]
        return 1


    def  onemove( self, x, u, xp, up):
        """One move of the twalk.  This is basically the raw twalk kernel.
           It is useful if the twalk is needed inside a more complex MCMC.

           Args:

               x, xp: two points WITHIN the support. Note: *each entry of x0 and xp0 must be different*.
               u, up: the value of the objective at x, and xp u=U(x), up=U(xp).

           Returns
           -------
           [y, yp, ke, A, u_prop, up_prop]: list
           y, yp:
               the proposed jump
           ke:
               The kernel used, 0=nothing, 1=Walk, 2=Traverse, 3=Blow, 4=Hop
           A:
               the M-H ratio
           u_prop, up_prop:
               The values for the objective func. at the proposed jumps

        """

        # Make local references for less writing
        n = self.n
        U = self.U
        Supp = self.Supp
        Fw = self.Fw

        ker = uniform() # To choose the kernel to be used
        ke = 1
        A = 0

        # Kernel nothing exchange x with xp, not used
        if ((0.0 <= ker) & (ker < Fw[0])):
            ke = 0
            y = xp.copy()
            up_prop = u
            yp = x.copy()
            u_prop = up
            # A is the MH acceptance ratio
            A = 1.0  #always accepted


        # The Walk move
        if ((Fw[0] <= ker) & (ker < Fw[1])):

            ke = 1

            dir = uniform()

            if ((0 <= dir) & (dir < 0.5)):  # x as pivot

                yp = self.SimWalk( xp, x)

                y = x.copy()
                u_prop = u

                if ((Supp(yp)) & (all(abs(yp - y) > 0))):
                    up_prop = U(yp)
                    A = exp(up - up_prop)
                else:
                    up_prop = None
                    A = 0; #out of support, not accepted

            else:  # xp as pivot

                y = self.SimWalk( x, xp)

                yp = xp.copy()
                up_prop = up

                if ((Supp(y)) & (all(abs(yp - y) > 0))):
                    u_prop = U(y)
                    A = exp(u - u_prop)
                else:
                    u_prop = None
                    A = 0; #out of support, not accepted


        # The Traverse move
        if ((Fw[1] <= ker) & (ker < Fw[2])):

            ke = 2
            dir = uniform()

            if ((0 <= dir) & (dir < 0.5)):  # x as pivot

                beta = self.Simbeta()
                yp = self.SimTraverse( xp, x, beta)

                y = x.copy()
                u_prop = u

                if Supp(yp):
                    up_prop = U(yp)
                    if (self.nphi == 0):
                        A = 1 #Nothing moved
                    else:
                        A = exp((up - up_prop) +  (self.nphi-2)*log(beta))
                else:
                    up_prop = None
                    A = 0 #out of support, not accepted
            else:            # xp as pivot

                beta = self.Simbeta()
                y = self.SimTraverse( x, xp, beta)

                yp = xp.copy()
                up_prop = up

                if Supp(y):
                    u_prop = U(y)
                    if (self.nphi == 0):
                        A = 1 #Nothing moved
                    else:
                        A = exp((u - u_prop) +  (self.nphi-2)*log(beta))
                else:
                    u_prop = None
                    A = 0 #out of support, not accepted

        # The Blow move
        if ((Fw[2] <= ker) & (ker < Fw[3])):

            ke = 3
            dir = uniform()

            if ((0 <= dir) & (dir < 0.5)):  # x as pivot
                yp = self.SimBlow( xp, x)

                y = x.copy()
                u_prop = u
                if ((Supp(yp)) & all(yp != x)):
                    up_prop = U(yp)
                    W1 = self.GBlowU( yp, xp,  x)
                    W2 = self.GBlowU( xp, yp,  x)
                    A = exp((up - up_prop) + (W1 - W2))
                else:
                    up_prop = None
                    A = 0 #out of support, not accepted
            else:  # xp as pivot
                y = self.SimBlow( x, xp)

                yp = xp.copy()
                up_prop = up
                if ((Supp(y)) & all(y != xp)):
                    u_prop = U(y)
                    W1 = self.GBlowU(  y,  x, xp)
                    W2 = self.GBlowU(  x,  y, xp)
                    A = exp((u - u_prop) + (W1 - W2))
                else:
                    u_prop = None
                    A = 0 #out of support, not accepted


        # The Hop move
        if ((Fw[3] <= ker) & (ker < Fw[4])):

            ke = 4
            dir = uniform()

            if ((0 <= dir) & (dir < 0.5)):  # x as pivot
                yp = self.SimHop( xp, x)

                y = x.copy()
                u_prop = u
                if ((Supp(yp)) & all(yp != x)):
                    up_prop = U(yp)
                    W1 = self.GHopU( yp, xp,  x)
                    W2 = self.GHopU( xp, yp,  x)
                    A = exp((up - up_prop) + (W1 - W2))
                else:
                    up_prop = None
                    A = 0 #out of support, not accepted
            else:  # xp as pivot
                y = self.SimHop( x, xp)

                yp = xp.copy()
                up_prop = up
                if ((Supp(y)) & all(y != xp)):
                    u_prop = U(y)
                    W1 = self.GHopU(  y,  x, xp)
                    W2 = self.GHopU(  x,  y, xp)
                    A = exp((u - u_prop) + (W1 - W2))
                else:
                    u_prop = None
                    A = 0 #out of support, not accepted

        return [y, yp, ke, A, u_prop, up_prop]



#
# Auxiliars for the kernels

    def SimWalk( self, x, xp):
        """Used by the Walk kernel"""
        aw = self.aw
        n = self.n

        phi = (uniform(size=n) < self.pphi) # parametrs to move
        self.nphi = sum(phi)
        z = zeros(n)

        for i in range(n):
            if phi[i]:
                u = uniform()
                z[i] = (aw/(1+aw))*(aw*u**2.0 + 2.0*u - 1.0)

        return x + (x - xp)*z


    def Simbeta(self):
        """Used by the Traverse kernel"""
        at = self.at
        if (uniform() < (at-1.0)/(2.0*at)):
            return exp(1.0/(at+1.0)*log(uniform()))
        else:
            return exp(1.0/(1.0-at)*log(uniform()))

    def SimTraverse( self,  x, xp, beta):
        n = self.n

        phi = (uniform(size=n) < self.pphi)
        self.nphi = sum(phi)

        rt = x.copy()
        for i in range(n):
            if (phi[i]):
                rt[i] = xp[i] + beta*(xp[i] - x[i])

        return rt



    def SimBlow( self, x, xp):
        """ Used by the Blow kernel"""
        n = self.n

        self.phi = (uniform(size=n) < self.pphi)
        self.nphi = sum(self.phi)

        self.sigma = max(self.phi*abs(xp - x))

        rt = x.copy()
        for i in range(n):
            if (self.phi[i]):
                rt[i] = xp[i] + self.sigma * normal()

        return rt


    def GBlowU( self, h, x, xp):
        nphi = self.nphi
        self.sigma = max(self.phi*abs(xp - x)) #recalculate sigma, but same phi
        if (nphi > 0):
            return (nphi/2.0)*log2pi + nphi*log(self.sigma) + 0.5*SqrNorm(h - xp)/(self.sigma**2)
        else:
            return 0



    def SimHop( self, x, xp):
        """ Used by the Hop kernel"""
        n = self.n

        self.phi = (uniform(size=n) < self.pphi)
        self.nphi = sum(self.phi)

        self.sigma = max(self.phi*abs(xp - x))/3.0

        rt = x.copy()
        for i in range(n):
            if (self.phi[i]):
                rt[i] = x[i] + self.sigma * normal()

        return rt


    def GHopU( self, h, x, xp): # It is actually equal to GBlowU!
        nphi = self.nphi
        self.sigma = max(self.phi*abs(xp - x))/3.0 #Recalculate sigma, but same phi

        if (nphi > 0): #Mistake until 20AUG2020, formely pivot was left to xp:
            return (nphi/2.0)*log2pi + nphi*log(self.sigma) + 0.5*SqrNorm(h - x)/(self.sigma**2) #
        else:
            return 0




#  Output analysis auxiliar methods

    def IAT( self, par=-1, start=0, end=0, maxlag=0):
        """Calculate the Integrated Autocorrelation Times of parameters par
           the default value par=-1 is for the IAT of the U's"""
        if (end == 0):
            end = self.T

        if (self.Acc[5] == 0):
            print("twalk: IAT: WARNING,  all propolsals were rejected!")
            print("twalk: IAT: Cannot calculate IAT, fixing it to the sample size.")
            return self.T

        iat = IAT( self.Output, cols=par, maxlag=maxlag, start=start, end=end)

        return iat


    def TS( self, par=-1, start=0, end=0):
        """Plot time series of parameter par (default = log f) etc."""
        if par == -1:
            par = self.n

        if (end == 0):
            end = self.T

        if (par == self.n):
            plot( list(range( start, end)), -1*self.Output[ start:end, par])
            ylabel("Log of Objective")
        else:
            plot( list(range( start, end)), self.Output[ start:end, par])
            ylabel("Parameter %d" % par)
        xlabel("Iteration")


    def Ana( self, par=-1, start=0, end=0):
        """Output Analysis, TS plots, acceptance rates, IAT, etc."""
        if par == -1:
            par = self.n

        if (end == 0):
            end = self.T

        print("Acceptance rates for the Walk, Traverse, Blow and Hop kernels:" + str(self.Acc[1:5]))
        print("Global acceptance rate: %7.5f" % self.Acc[5])

        iat = self.IAT( par=par, start=start, end=end)
        print("Integrated Autocorrelation Time: %7.1f, IAT/n: %7.1f" % (iat, iat/self.n))

        self.TS( par=par, start=start, end=end)

        return iat


    def Hist( self, par=-1, start=0, end=0, g=(lambda x: x[0]), xlab=None, bins=20, density=False):
        """Basic histograms and output analysis.  If par=-1, use g.
           The function g provides a transformation to be applied to the data,
           e.g. g=(lambda x: abs(x[0]-x[1]) would plot a histogram of the distance
           between parameters 0 and 1, etc."""

        if (end == 0):
            end = self.T

        if (par == -1):
            ser = zeros(end-start)
            for it in range(end-start):
                ser[it] = g(self.Output[ it+start, :-1])
            if (xlab == None):
                xlab = "g"
        else:
            ser = self.Output[ start:end, par]
            if (xlab == None):
                xlab = "parameter %d" % (par,)

        xlabel(xlab)
        print("Mean for %s= %f" % ( xlab, mean(ser)))
        return hist( ser, bins=bins, density=density)



    def Save( self, fnam, start=0, end=-1, thin=1):
        """Saves the Output as a text file, starting at start (burn in), with thinning (thin)."""

        print("Saving output, all pars. plus the U's in file", fnam)

        savetxt( fnam, self.Output[ start:end:thin,:])



    def Load( self, fnam, start=0, thin=1):
        """Loads the Output from a text file, typically written with the Save method.
           It will overwrite any other twalk output.  Updates the dimension n and the sample size T."""

        print("Loading output from file", fnam)

        self.Output = loadtxt(fnam)
        self.T, self.n = self.Output.shape
        self.n -= 1


# A simple Random Walk M-H
    def RunRWMH( self, T, x0, sigma):
        """Run a simple Random Walk M-H"""

        sec = time() # last time we sent a message
        print("pytwalk: This is the Random Walk M-H running with %d iterations." % T)
        # Local variables
        x = x0.copy()
        if not(self.Supp(x)):
            print("pytwalk: ERROR, initial point x0 out of support.")
            return 0
        self.T = T

        u = self.U(x)
        n = self.n

        sec2 = time() # last time we sent a message
        print("       " + Remain( T, 2, sec, sec2))

        # Set the array to place the iterations and the U's
        self.Output = zeros((T+1, n+1))
        self.Acc = zeros(6)

        # Make local references for less writing
        Output = self.Output
        U = self.U
        Supp = self.Supp
        Acc = self.Acc

        Output[ 0, 0:n] = x.copy()
        Output[ 0, n] = u

        j1=1
        j=0

        y = x.copy()
        for it in range(T):
            y = x + normal(size=n)*sigma # each entry with sigma[i] variance
            if Supp(y):        # If it is within the support of the objective
                uprop = U(y)   # Evaluate the objective
                if (uniform() < exp(u-uprop)):
                    x = y.copy()   # Accept the propolsal y
                    u = uprop
                    Acc[5] += 1

            # Estimate the remaing time, every 2**j1 iterations
            if ((it % (1 << j1)) == 0):

                j1 += 1
                j1 = min( j1, 10)  # check the time at least every 2^10=1024 iterations
                ax = time()
                if ((ax - sec2) > (1 << j)*self.WAIT): # Print an estimation every WAIT*2**j
                    print("pytwalk: %10d iterations so far. " % (it,) + Remain( T, it, sec, ax))
                    sec2 = ax
                    j += 1
                    j1 -= 1 # check the time as often

            Output[it+1,0:n] = x
            Output[it+1,n] = u

        if (Acc[5] == 0):
            print("pytwalk: WARNING,  all propolsals were rejected!")
            return 0

        Acc[5] /= T;
        return 1








#
# Auxiliary functions to calculate Integrated autocorrelation times of a time series



def AutoCov( Ser, c, la, T=0):
    """
    Calculates an autocovariance 2x2 matrix at lag la in column c of matrix Ser with T rows. The variances of each series are in the diagonal and the (auto)covariance in the off diag.

    :param Ser: matrix
    :param c: column
    :param la: lag
    :param T: rows, defaults to 0



    """
    if (T == 0):
        T = shape(Ser)[0]  # Number of rows in the matrix (sample size)

    return cov( Ser[0:(T-1-la), c], Ser[la:(T-1), c], bias=1)




def AutoCorr( Ser, cols=0, la=1):
    """
    Calculates the autocorrelation from lag 0 to lag la of columns cols (list)
    for matrix Ser.

    :param Ser: matrix
    :param cols: columns, defaults to 0
    :param la: lag, defaults to 1


    """
    T = shape(Ser)[0]  # Number of rows in the matrix (sample size)

    ncols = shape(mat(cols))[1] # Number of columns to analyse (parameters)

    #if ncols == 1:
    #    cols = [cols]

    # Matrix to hold output
    Out = matrix(ones((la+1)*ncols)).reshape( la+1, ncols)

    for c in range(ncols):
        for l in range( 1, la+1):
            Co = AutoCov( Ser, cols[c], l, T)
            Out[l,c] = Co[0,1]/(sqrt(Co[0,0]*Co[1,1]))

    return Out



def MakeSumMat(lag):
    """
    Makes an upper band matrix of ones to add the autocorrelation matrix.

    :code:`gamma = auto[2*m+1,c]+auto[2*m+2,c]` etc.

    :code:`MakeSumMat(lag) * AutoCorr( Ser, cols=c, la=lag)` to make the gamma matrix

    Args:
        lag (int): lag


    """
    rows = (lag)//2   # Integer division!
    Out = mat(zeros([rows,lag], dtype=int))

    for i in range(rows):
        Out[i,2*i] = 1
        Out[i,2*i+1] = 1

    return Out



def Cutts(Gamma):
    """Finds the cutting time when gammas become negative.
    """
    cols = shape(Gamma)[1]
    rows = shape(Gamma)[0]
    Out = mat(zeros([1,cols], dtype=int))
    Stop = mat(zeros([1,cols], dtype=bool))

    if (rows == 1):
        return Out

    i = 0
    #while (not(all(Stop)) & (i < (rows-1))):
    for i in range(rows-1):
        for j in range(cols):  # while Gamma stays positive and decreasing
            if (((Gamma[i+1,j] > 0.0) & (Gamma[i+1,j] < Gamma[i,j])) & (not Stop[0,j])):
                Out[0,j] = i+1 # the cutting time for colomn j is i+i
            else:
                Stop[0,j] = True
        i += 1


    return Out



def AutoMaxlag( Ser, c, rholimit=0.05, maxmaxlag=20000):
    """ Automatically find a maxlag for IAT calculations.

    :param Ser: matrix

    :param c: c

    :param rholimit: defaults to 0.05

    :param maxmaxlag: defaults to 20000



    """
    Co = AutoCov( Ser, c, la=1)
    rho = Co[0,1]/Co[0,0]  # lag one autocorrelation

    # if autocorrelation is like exp(- lag/lam) then, for lag = 1
    lam = -1.0/log(abs(rho))

    # Our initial guess for maxlag is 1.5 times lam (eg. three times the mean life)
    maxlag = int(floor(3.0*lam))+1

    # We take 1% of lam to jump forward and look for the
    # rholimit threshold
    jmp = int(ceil(0.01*lam)) + 1

    T = shape(Ser)[0]  # Number of rows in the matrix (sample size)

    while ((abs(rho) > rholimit) & (maxlag < min(T//2,maxmaxlag))):
        Co = AutoCov( Ser, c, la=maxlag)
        rho = Co[0,1]/Co[0,0]
        maxlag = maxlag + jmp
        #print("maxlag=", maxlag, "rho", abs(rho), "\n")

    maxlag = int(floor(1.3*maxlag));  #30% more

    if (maxlag >= min(T//2,maxmaxlag)): #not enough data
        fixmaxlag = min(min( T//2, maxlag), maxmaxlag)
        print("AutoMaxlag: Warning: maxlag= %d > min(T//2,maxmaxlag=%d), fixing it to %d" % (maxlag, maxmaxlag, fixmaxlag))
        return fixmaxlag

    if (maxlag <= 1):
        fixmaxlag = 10
        print("AutoMaxlag: Warning: maxlag= %d ?!, fixing it to %d" % (maxlag, fixmaxlag))
        return fixmaxlag

    print("AutoMaxlag: maxlag= %d." % maxlag)
    return maxlag



def IAT( Ser, cols=-1,  maxlag=0, start=0, end=0):
    """Find the IAT.


    :param cols: defaults to -1

    :param maxlag: defaults to 0

    :param start: defaults to 0

    :param end: defaults to 0

    """

    ncols = shape(mat(cols))[1] # Number of columns to analyse (parameters)
    if ncols == 1:
        if (cols == -1):
            cols = shape(Ser)[1]-1 # default = last column
        cols = [cols]

    if (end == 0):
        end = shape(Ser)[0]

    if (maxlag == 0):
        for c in cols:
            maxlag = max(maxlag, AutoMaxlag( Ser[start:end,:], c))

    #print("IAT: Maxlag=", maxlag)

    #Ga = MakeSumMat(maxlag) * AutoCorr( Ser[start:end,:], cols=cols, la=maxlag)

    Ga = mat(zeros((maxlag//2,ncols)))
    auto = AutoCorr( Ser[start:end,:], cols=cols, la=maxlag)

    # Instead of producing the maxlag/2 X maxlag MakeSumMat matrix, we calculate the gammas like this
    for c in range(ncols):
        for i in range(maxlag//2):
            Ga[i,c] = auto[2*i,c]+auto[2*i+1,c]

    cut = Cutts(Ga)
    nrows = shape(Ga)[0]

    ncols = shape(cut)[1]
    Out = -1.0*mat(ones( [1,ncols] ))

    if any((cut+1) == nrows):
        print("IAT: Warning: Not enough lag to calculate IAT")

    for c in range(ncols):
        for i in range(cut[0,c]+1):
            Out[0,c] += 2*Ga[i,c]

    return Out


