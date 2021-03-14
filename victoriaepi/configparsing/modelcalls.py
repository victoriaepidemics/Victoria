#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:08:51 2020

"""

from dataclasses import dataclass
import collections
import dacite
from .generic import codeconfig_getvars
from . import calls_dir
from datetime import date
import re
import typing
import warnings
import os
from tempfile import gettempdir

@dataclass
class rrate ():
    """Residence rate 1/day, names, and Erlang series."""

    residenceRate: str
    """A number or expression."""
    name: str
    """Name of the rate, eg. r'\sigma_1'."""
    erlang: int
    """Erlang series."""

    def __post_init__(self):
        assert isinstance(eval(self.residenceRate),float)|isinstance(eval(self.residenceRate),int),"Residence rate must be numeric after eval(\""+self.residenceRate+"\")"
        assert self.erlang>0
    def __repr__(self):
        return("["+str(self.residenceRate)+", r'"+self.name+"', "+str(self.erlang)+"]")

class Rrates (collections.OrderedDict):
    """Residence rates"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for k in self.keys():
            super().__setitem__(k,
                                dacite.from_dict(data_class=rrate, data=self[k]))
    def __repr__(self):
        """
        Creates default dictionary representation
        Residence rates 1/day, names and Erlang series:
        R_rates=\
        { 'E':[1/5, r'\sigma_1',  4],\
          'I':[1/4, r'\sigma_2',  3]}


        """
        res = "{"
        for k in self.keys():
            res+=" '"+str(k)+"':"+str(self[k])+","
        res=res[:-1]+" }"
        return res
    def AuxiliarDefinitionsForSolve(self):
        res=""
        for k in self.keys():
            nombrevar=re.sub(r"[^\w]","",self[k].name)
            res+="self."+nombrevar+" = self.par[self.T.ConvertVar('"+k+"')]\n"
        return res

class strevaldate (str):
    """ Simple str extension used to validate date strings."""
    def __new__(cls, val, *args, **kwargs):
        return super(strevaldate, cls).__new__(cls, val)
    def __init__(self,val):
        assert isinstance(eval(val),date),val + " is not valid date."


@dataclass
class Zone:
    """
    General information for the metro zone or region to be analyzed.

    The actual instance of the covid_mcmc object is stored in the last item of the list:
   
    * id 
    * Name   
    * num_relax_days  
    * Population
    * init date
    * intervention and relax dates

    """

    id: str
    name: str
    """ Name of the region. """
    num_relax_days: int
    population: float
    init_date: str
    """eg. date(2020, 2, 27)"""
    intervention_n_relax_dates: typing.List
    """See ``intervention_day`` in :class:``victoriaepi.ama.ama2``. However here is no need to end with ``None``
    eg ``[ "date(2020, 3, 22)", "date(2020, 4, 3)", "date(2020, 5, 10)", "date(2020, 6, 3)"]``.
    """

    def __post_init__(self):
        self.init_date = strevaldate(self.init_date)
        self.intervention_n_relax_dates = \
            [strevaldate(x) for x in self.intervention_n_relax_dates]

    def __repr__(self):
        """
        {"9-01": "Mexico city", 2, 21.942666e6, date(2020, 2, 27),  date(2020, 3, 22), date(2020, 4,  3), date(2020, 5, 10), date(2020, 6, 3), None]}
        :return: code representation
        :rtype: str
        """
        res = "{" + '"' + self.id + '" : ["' + self.name + '", ' + str(self.num_relax_days)\
            + ', ' + str(self.population) + ', ' + self.init_date + ', ' \
            + ', '.join(self.intervention_n_relax_dates) + ", None]}"
        return res

    @classmethod
    def from_dict(cls, odic):
        """Initialize instance from dictionary."""
        return dacite.from_dict(data_class=cls, data=odic)


@dataclass
class ModelCall:
    """Creates model call."""

    Region: str
    """
    Region name

    .. deprecated:: 0.0.2
        Use Zone instead

    """
    complexity: str
    """ See :any:`victoriaepi.configparsing.config.ModelConfig.complexity`"""
    className: str
    """ See :any:`victoriaepi.configparsing.config.ModelConfig.className`"""
    moduleName: str
    """ See :any:`victoriaepi.configparsing.config.ModelConfig.moduleName`"""
    T: int
    """ Number of iterations. Should be greater than burnin. See :any:`victoriaepi.victoria.mcmc.RunMCMC`."""
    data_fnam: str
    """Data file name, workdir/data/data_fnam. Make sure to process data into a vertical text array."""
    out_fnam: str
    """MCMC output file name, without .txt, workdir/output/out_fnam + '.txt'"""
    init_index: int
    """Day number from data start, where to start the plot."""
    burnin: int
    """ Burnin iterations. See :any:`victoriaepi.victoria.mcmc.RunMCMC`."""
    init_date: str  # TODO: ok to override Zone
    """
    Date at init_index. eg "date(2020, 4, 1)"

    .. deprecated:: 0.0.2
        Use Zone instead

    """
    trim: int
    """How many data points to trim."""
    pred: typing.Union[int, str]
    """Number of days to predict. Can be an expression."""

    workdir: str  # ='./../'
    """ string with ./data/ and ./output/ folders  and data files."""
    R_rates: collections.OrderedDict  # Later cast as Rrates
    """ :class: `victoriaepi.configparsing.modelcalls.Zone`"""
    exit_probs: str = ""
    """exit probabilities."""
    plot_fit: str = ""
    """Plots basic results of the MCMC."""
    daysdelay: int = None
    """Delay for symptoms onset"""
    instanceName: str = None
    """
    Gives a specific name to the instance.

    .. deprecated:: 0.0.2
        Not very useful.

    """
    plotpred: typing.Union[int, str] = None
    """Number of predicted days to plot. Can be an expression."""
    N: int = None
    """
    Population size

    .. deprecated:: 0.0.2
        Use Zone instead

    """
    Zone: collections.OrderedDict = None
    """See :class: `victoriaepi.configparsing.modelcalls.Zone`"""

    plottingQuantiles: typing.List[typing.Union[int, float]] = None
    """ Plotting quantiles eg [10, 25, 50, 75, 90]"""
    Pobs_D: float = None
    """probability of recording a death"""
    Pobs_I: float = None
    """probability of recording an infection"""
    ngrp: int = None
    """ number of age groups"""
    age_prop: typing.List[float] = None
    """list of age proportions of length ngrp"""
    age_groups: typing.List[int] = None
    """ Ages of the age groups eg. [0, 25, 50, 65, 100]"""
    Int_M: str = None
    """ Contact matrix. ngrp X ngrp array with the interaction. Can be an expression such as ``np.ones((ngrp,ngrp))/ngrp``."""
    age_prop: typing.List[float] = None
    """ Population pyramid"""
    def __post_init__(self):
        """Model Call."""
        with open(calls_dir + self.complexity + ".py", "r") as f:
            self.code_cfg = f.read()
        # Late castings
        self.R_rates = Rrates(self.R_rates)
        # Consistency checks
        # TODO: def basic checks
        assert self.N > 1, "N should be >1"  # TODO: could be in Zone
        assert self.trim <= 0, "trim should be negative or zero"
        self.init_date = strevaldate(self.init_date)
        assert os.path.exists(self.workdir), "Workdir " + self.workdir + "does not exist."
        assert os.path.isfile(self.workdir + 'data/' + self.data_fnam), \
            "Data_fnam not found: "+ self.workdir + 'data/' + self.data_fnam
        if self.workdir[-1] != "/":
            self.workdir = self.workdir + "/"
        # Additional definitions depending on model complexity
        if self.Zone is not None:
            self.Zone = Zone.from_dict(self.Zone)
            self.Zone_id = "'" + self.Zone.id + "'"
        compile(str(self.pred), "< pred = " + str(self.pred) + ">", 'exec')
        compile(str(self.plotpred), "< plotpred = " + str(self.plotpred) + ">", 'exec')
        # TODO: add remaining vars

    def gencode(self):
        """Create Model Call code."""

        fullcode = self.code_cfg
        variables = codeconfig_getvars(fullcode)
        if len(variables) > 0:
            variables.sort(key=lambda x: len(x), reverse=True)
            for va in variables:
                if eval("self." + va[1:]) is None or eval("self." + va[1:]) == "":
                    warnings.warn("Inserting None for " + va, UserWarning)
                fullcode = fullcode.replace(va, str(eval("self." + va[1:])))
        # testing for sintax errors
        #compile(fullcode, "<Test ModelCall Code>\n" + fullcode, 'exec')
        tmpfile = os.path.join(gettempdir(), "TestModelCallCode." + self.complexity + ".py")
        with open(tmpfile, "w") as f: f.write(fullcode)
        compile(fullcode, tmpfile, 'exec')
        return fullcode

    def unusedVars(self):
        """Find defined, yet unused variables."""
        fullcode = self.code_cfg
        variables = set([x[1:] for x in codeconfig_getvars(fullcode)])
        exceptions = set(['complexity', 'code_cfg'])
        clsvars = set(vars(self).keys())
        nones = set(filter(lambda x: self.__dict__[x] is None, clsvars))
        nones = nones.union(set(filter(lambda x: str(self.__dict__[x]) == "", clsvars)))
        unused = clsvars - variables - exceptions - nones
        return unused
