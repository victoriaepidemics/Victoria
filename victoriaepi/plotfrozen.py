# -*- coding: utf-8 -*-
"""
Plots the instanced (frozen) univariate probability distribution.

"""



from numpy import arange, linspace
from pylab import subplots

def PlotFrozenDist( dist, q=1e-6, N=200, linestyle='-', marker='o', color='k', ax=None, **kwargs):
    """
    Plots the instanced (frozen) univariate probability distribution `dist`, the range
    to plot is chosen from the `q` quantile to the `1-q` quantile.

    Automatically checks if `dist` is discrete or continuous.  It is assumed that all discrete distributions
    are supported in the integers.  For continuous `dist` the plot range is divided
    in `N` equal parts for plotting. q may be a size 2 iterable such that the plotting
    range is from the `q[0]` quantile to the `q[1]` quantile. The rest  of the arguments
    control plotting and are passed to the function plot.

    Args:
        dist : univariate probability distribution
        q : quantile, float or size 2 iterable
        N (int): For continuous `dist` the plot range is divided into `N` equal parts for plotting
        *args, **kwargs: passed to the function plot

    Example::

        from plotfrozen import PlotFrozenDist, ExamplePlotFrozenDist
        from pylab import figure, subplot
        from scipy.stats import norm, gamma, poisson, skellam
        from  matplotlib.pyplot import vlines

        # This same example is reproduced by calling ExamplePlotFrozenDist().

        figure()
        subplot(421)
        PlotFrozenDist( norm() )
        subplot(422)
        PlotFrozenDist( norm( loc=10, scale=0.5) )
        subplot(423)
        PlotFrozenDist( gamma(3) )
        subplot(424)
        ga = gamma( 3, loc=5)
        PlotFrozenDist( ga, q=1e-6, color='g')
        mn = float(ga.stats()[0]) # add a red vline at the mean
        vlines( mn, 0, ga.pdf(mn),  colors='r', linestyles='-')
        subplot(425)
        po = poisson(5)
        PlotFrozenDist(po)
        subplot(426)
        PlotFrozenDist( poisson(10), color='r', ms=3)
        subplot(427)
        PlotFrozenDist( skellam( 5, 10), marker='*')
        subplot(428)
        PlotFrozenDist( skellam( 100, 10))




    """

    if (isinstance(q, float)):
        q0 = q
        q1 = 1.0-q
    else:
        q0 = q[0]
        q1 = q[1]

    if ax==None:
        fig, ax = subplots( figsize=(5,5))

    # Check the type of distribution:
    continuous = True
    try:
        dist.pdf(0)
    except:
        continuous = False

    if (continuous):
        x = linspace( dist.ppf(q0), dist.ppf(q1), N)
        y = dist.pdf(x)
        ax.plot(x, y, linestyle=linestyle, color=color, **kwargs)
    else:
        x = arange( dist.ppf(q0), dist.ppf(q1), dtype=int)
        y = dist.pmf(x)
        ax.plot(x, y, marker, color=color, **kwargs)
        ax.vlines(x, 0, y,  colors=color, linestyles=linestyle)

    return ax


def ExamplePlotFrozenDist():
    """Example of PlotFrozenDist.
    See :meth:`plotfrozen.PlotFrozenDist`
    """

    from pylab import figure, subplot
    from scipy.stats import norm, gamma, poisson, skellam
    from  matplotlib.pyplot import vlines
    figure()
    subplot(421)
    PlotFrozenDist( norm() )
    subplot(422)
    PlotFrozenDist( norm( loc=10, scale=0.5) )
    subplot(423)
    PlotFrozenDist( gamma(3) )
    subplot(424)
    ga = gamma( 3, loc=5)
    PlotFrozenDist( ga, q=1e-6, color='g')
    mn = float(ga.stats()[0]) # add a red vline at the mean
    vlines( mn, 0, ga.pdf(mn),  colors='r', linestyles='-')
    subplot(425)
    PlotFrozenDist(poisson(5.3))
    subplot(426)
    PlotFrozenDist( poisson(10), color='r', ms=3) # Pass other plotting arguments
    subplot(427)
    PlotFrozenDist( skellam( 5.3, 10), marker='*')
    subplot(428)
    PlotFrozenDist( skellam( 100, 10))

