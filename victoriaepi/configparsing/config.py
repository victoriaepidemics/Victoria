#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 21:32:27 2020

"""
from datetime import date
# import pathlib
import dacite
from dataclasses import dataclass
import collections
# import typing
import warnings
import re
# import string
import hjson
# import inspect
import os
import textwrap
from tempfile import gettempdir

from .generic import codeconfig_getvars
from .modelclass import ModelMatrix, plottingMethods, solvingMethods, functionCollection
from .modelcalls import ModelCall
from . import models_dir


@dataclass
class ModelConfig:
    """
    Parse config files.

    :return: ModelConfig
    :rtype: ModelConfig

    """
    className: str
    """str: class name"""
    complexity: str
    """Represents the complexity of the model to be configured. Can be one of `basic`, `intervention` or `grp`."""
    moduleName: str #
    """Module name. This will be used for naming the output files."""
    odeint: str
    """str: Raw code to fill odeint method See :meth:`victoriaepi.seid.odeint` """
    baseModelVersion: str #= "0.0.0"
    """ Base Model Version for internal reference only."""
    ModelCall: collections.OrderedDict #
    """ Starts as a ``collections.OrderedDict`` but is later cast as a :class:`victoriaepi.configparsing.modelcalls.ModelCall` instance."""
    ModelMatrix: collections.OrderedDict #
    """ Starts as a ``collections.OrderedDict`` but is later cast as a
    :class:`victoriaepi.configparsing.modelclass.ModelMatrix` instance. """
    num_pars: int
    """int: number of parameters in model"""
    solvingMethods: collections.OrderedDict
    """ Starts as a ``collections.OrderedDict`` but is later cast as a
    :class:`victoriaepi.configparsing.modelclass.solvingMethods` instance."""
    plottingMethods: collections.OrderedDict
    """ Starts as a ``collections.OrderedDict`` but is later cast as a
    :class:`victoriaepi.configparsing.modelclass.plottingMethods` instance."""

    factor_foi: str = ""
    """force of infection"""
    additionalMethods: collections.OrderedDict = None
    """Optional. Starts as a ``collections.OrderedDict`` but is later cast as a
    :class:`victoriaepi.configparsing.modelclass.functionCollection` instance.
    This is also useful for overriding :any:`victoriaepi.victoria.mcmc`'s methods.
    """
    additionalCallMethods: collections.OrderedDict = None
    """Optional. Starts as a ``collections.OrderedDict`` but is later cast as a
    :class:`victoriaepi.configparsing.modelclass.functionCollection` instance.
    Can be used to add additional methods to the call file.
    """


    def __post_init__(self):

        with open(models_dir+self.complexity+".py","r") as f:
            self.code_cfg = f.read()
        #Preparing forward map config TODO: def
        self.ModelMatrix["className"]=self.className
        self.ModelMatrix["complexity"]=self.complexity
        self.ModelMatrix = dacite.from_dict(data_class=ModelMatrix, data=self.ModelMatrix)  # TODO: ModelMatrix.fromdict

        #Preparing Call config TODO: def
        self.ModelCall["className"]=self.className
        self.ModelCall["moduleName"]=self.moduleName
        self.ModelCall["complexity"]=self.complexity
        self.ModelCall = dacite.from_dict(data_class = ModelCall, data = self.ModelCall)  # TODO: ModelCall.fromdict


        if self.__isMarkerPresent("SelfMaskCode"):
            self.SelfMaskCode = self.ModelMatrix.genMaskcode(self.__getcodeindent("SelfMaskCode"))
            """str: masking code"""
        #Basic forward map configuration
        self.defineModelMatrix = textwrap.indent(self.ModelMatrix.gencode(),
                                               self.__getcodeindent("defineModelMatrix"))
         #Late assert R_rates are in names
        for k in self.ModelCall.R_rates.keys():
            assert k in self.ModelMatrix.splitnames, k + " used in rates but not in names"
        if self.__isMarkerPresent("AuxiliarDefinitionsForSolve"):
            self.AuxiliarDefinitionsForSolve = self.__autoindent(\
                self.ModelCall.R_rates.AuxiliarDefinitionsForSolve(), "AuxiliarDefinitionsForSolve")
        # Additional configuration based on complexity
        if self.ModelMatrix.ExitFractions is not None:
            self.exitProbsDefinitions = self.ModelMatrix.ExitFractions.\
                exitprobsDefinitions(self.__getcodeindent("exitProbsDefinitions"))
            self.inlineargsequals = self.ModelMatrix.ExitFractions.inlineargsequals()
            self.ModelCall.exit_probs = self.ModelMatrix.ExitFractions.ModelCall_exit_probs()

        #Solving methods
        self.solvingMethods = solvingMethods(self.solvingMethods,
                                                 baseIndent = self.__getcodeindent("solvingMethods"))
        #Plotting methods
        self.plottingMethods = plottingMethods(self.plottingMethods,
                                                 baseIndent = self.__getcodeindent("plottingMethods"))
        #Additional methods and overrides
        if self.additionalMethods is not None:
            self.additionalMethods = \
                functionCollection(self.additionalMethods,\
                                   baseIndent = self.__getcodeindent("additionalMethods"))
        else:
            self.additionalMethods = ""

        #Additional methods and overrides
        if self.ModelCall.ngrp is not None:
            assert self.ModelCall.ngrp > 1
            assert self.ModelMatrix.ExitFractions.ngrp == self.ModelCall.ngrp


            #call gencode for testing

        self.gencode()

    def __autoindent(self, code, marker):
        firstline = code.splitlines()[0]
        return textwrap.indent(code,\
                self.__getcodeindent(marker),lambda x: not firstline in x)

    def __isMarkerPresent (self,marker):
        p = re.compile(r"^(\s*)\$"+marker+r"[^\w]*$", flags = re.MULTILINE)
        result = p.findall(self.code_cfg)
        if len(result)==1:
            return True
        return False
    def __getcodeindent(self,marker):
        p = re.compile(r"^(\s*)\$"+marker+r"[^\w]*$",flags = re.MULTILINE)
        result = p.findall(self.code_cfg)
        assert len(result)==1,\
            "No single result for $"+ marker+": "+str(len(result)) +"\n"\
                +str(result)#+"\n\n"+self.code_cfg
        indent = str(result[0])
        return indent.replace("\n","")

    def gencode(self):
        """
        Create Class code

        Returns
        -------
        Raw Python code

        """
        fullcode = self.code_cfg
        variables = codeconfig_getvars(fullcode)
        if len(variables)>0:
            variables.sort(key = lambda x: len(x), reverse = True)
            for va in variables:
                if eval("self."+va[1:]) is None:
                    warnings.warn("Inserting None for "+va, UserWarning)
                fullcode = fullcode.replace(va, str(eval("self."+va[1:])))
        #testing for sintax errors
        tmpfile = os.path.join(gettempdir(),"TestModelClassCode." + self.complexity + ".py")
        with open(tmpfile,"w") as f:f.write(fullcode)
        compile(fullcode, tmpfile, 'exec')
        return(fullcode)

    def unusedVars(self):
        """Find defined, yet unused variables."""
        fullcode = self.code_cfg
        variables = set([x[1:] for x in codeconfig_getvars(fullcode)])
        exceptions = set(['ModelCall', 'ModelMatrix', 'additionalCallMethods', \
                           'baseModelVersion', 'code_cfg', 'complexity', 'moduleName'])
        nones = set(filter (lambda x: self.__dict__[x] is None, variables))
        nones = nones.union(set(filter(lambda x: len(str(self.__dict__[x])) == 0, variables)))
        variables = variables
        clsvars = set(vars(self).keys())
        nones = set(filter (lambda x: self.__dict__[x] is None, clsvars ))
        nones = nones.union(set(filter (lambda x: str(self.__dict__[x]) == "", clsvars)))
        unused = clsvars - variables - exceptions - nones
        # if unused:
        #     for v in unused:
        #         print(v+":")
        #         print(self.__dict__[v])
        return unused

    def writecode(self):
        """
        Create Class and Call python files, using the following structure::

            print ("creating class file: "+self.moduleName+".py")
            print ("creating call file: "+"call_"+self.moduleName+".py")

        Returns
        -------
        None.

        """
        print ("creating class file: "+self.moduleName+".py")
        with open(self.moduleName+".py","w") as f:f.write(self.gencode())
        if self.unusedVars():
            warnings.warn("Unused variables: "+str(self.unusedVars()))
        print("ok")
        print ("creating call file: "+"call_"+self.moduleName+".py")
        with open("call_"+self.moduleName+".py","w") as f:f.write(self.ModelCall.gencode())
        if self.ModelCall.unusedVars():
            warnings.warn("Unused variables: "+str(self.ModelCall.unusedVars()))
        print("ok")

    @classmethod
    def fromfilename(cls, filename):
        """
        Initialize instance from configuration file.

        Parameters
        ----------
        filename : str
            filename with valid hjson configuration.

        Returns
        -------
        ModelConfig

        """
        with open(filename,"r") as f:rawcfg = hjson.loads(f.read())
        return dacite.from_dict(data_class = cls, data = rawcfg)
