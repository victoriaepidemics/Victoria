#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:08:57 2020

"""

from dataclasses import dataclass
import collections
import dacite
import textwrap
import re
from .generic import codeconfig_getvars, prettifyvarname
from . import buildingblocks_dir, indentStr
import warnings
import typing

class EFractions (collections.OrderedDict):
    """
    Handles the Exit fractions.
    See :class:`victoria.AuxMatrix`
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.foundlist = False
        print ("Define: self.foundlist = "+str())
        for k in self.keys():
            frac = eval(str(self[k]))
            if isinstance(frac,float):
                assert frac<=1
                assert frac>=0
            elif isinstance(frac,typing.List):
                self.foundlist = True
                assert [x<=1 for x in frac] == [True] * len(frac), str(frac) + "<=1 not True"
                assert [x>=0 for x in frac] == [True] * len(frac), str(frac) + ">=0 not True"
            else:
                # Expecting a list, not leaving much options arbitrary code for now.
                assert isinstance(frac, typing.List), k + ": " + str(frac) + "not valid"
        if self.foundlist:
            lens = [len(self[k]) for k in self.keys()]
            assert lens == [lens[1]] * len(lens), "Exit fractions have different lengths: \n" + str({key:len(value) for (key,value) in self.items()})
            self.ngrp = lens[0]
        # TODO: post init asserts
    def __repr__(self):
        # FIX: unused?
        res="{----------------------------"
        for k in self.keys():
            res+=" '"+str(k)+"':"+str(self[k])+","
        res=res[:-1]+" }"
        return(res)
    def inlineargs(self):
        """ Inline argument definition.
        Outputs a string like this "f=5.0/100.0, g=0.03, h=0.4, i=0.5, "
        to be used for function definition
        """
        proportions_def=""
        for k in self.keys():
            proportions_def+=k+"="+str(self[k])+", "
        #proportions_def=proportions_def[:-2]
        return proportions_def
    def inlineargsequals(self):
        """ Inline argument definition.'
        Outputs a string like this "f=f, g=g, h=h, i=i, "
        to be used for Model_$className function calling
        if groups are present the output is like
        if ngroup
        f=self.f[grp], g=self.g[grp], h=self.h[grp], i=self.i[grp]
        """
        if self.foundlist:
            proportions_def=""
            for k in self.keys():
                proportions_def+=k+"=self."+k+"[grp], "
            #proportions_def=proportions_def[:-2]
            return proportions_def
        # in case exit probabilities are plain fractions
        proportions_def=""
        for k in self.keys():
            proportions_def+=k+"="+k+", "
        #proportions_def=proportions_def[:-2]
        return proportions_def
    def exitprobsDefinitions(self, baseIndent):
        """ Write exit probabilities definitions for Init_fm_matrix method.

            Outputs a string like::

                 f, g = exit_probs
                 self.f = f
                 self.g = g
        """
        res = ", ".join(self.keys()) + " = exit_probs\n"
        for k in self.keys():
            res += "self."+k+" = "+k+"\n"
        res = textwrap.indent(res, baseIndent)
        res = res.replace(baseIndent,"", 1)
        return res
    def ModelCall_exit_probs (self):
        res="["
        for k in self.keys():
            res+=str(self[k])+", "
        res=res[:-2]+"]"
        return(res)

@dataclass
class ModelMatrix:
    className: str
    """ See :any:`victoriaepi.configparsing.config.ModelConfig.className`"""
    complexity: str
    """ See :any:`victoriaepi.configparsing.config.ModelConfig.complexity`"""
    prnDefault: str
    """ prn -- to print or not the matrix after each operation, change default behavior."""
    names: str
    """ Names used in :any:`victoriaepi.victoria.AuxMatrix` """
    BaseVar: str
    """Define Base Variable. see :meth:`victoriaepi.victoria.AuxMatrix.BaseVar` """
    Exits: collections.OrderedDict
    """Dictionary defining connections between variables. The key represents the source variable and the value is a dict.

    Example::

        Exits:{
          #AuxMatrix.Exit( 'S', 'E')
          "S" : {"E":{}}
          #AuxMatrix.SplitExit( 'E', 'I^A', 'I^S', 1-f, prob_symb=['1-f','f'])
          "E" : {
              "I^A":{
                  "prob":1-f
                  "prob_symb":'1-f'
                  }
              "I^S":{
                  "prob":f
                  "prob_symb":'f'}
                  }
          #AuxMatrix.NoExit('R')
          "R" : {}
          #AuxMatrix.NoExit('D')
          "D" : {}
          }

    see also :any:`victoriaepi.victoria.AuxMatrix`
    """
    ExitFractions: collections.OrderedDict = None
    """Defines the exit fractions.

    Example::

        ExitFractions:{
            #Each can be an expression or a List of expressions
            # fraction of severe infections
            f:[0.05, 0.05, 0.05, 0.05]
            # fraction of severe infections that require hospitalization
            g:[0.004     , 0.0348    , 0.12333333, 0.2396    ]
            # fraction of hospitations that require ICU
            h:[0.05      , 0.0552    , 0.17266667, 0.5112    ]
            # fraction of ICU patients who die
            i:[0.5, 0.5, 0.5, 0.5]
            }

    see also :any:`victoriaepi.victoria.AuxMatrix`
    This is later cast as :class:`victoriaepi.configparsing.modelclass.EFractions`

    """
    def __post_init__(self):
        self.splitnames = self.names.split()
        # consistency checks
        assert self.BaseVar in self.splitnames, "BaseVar not in names"
        for k in set([kk for k in self.Exits for kk in [k]+list(self.Exits[k])]):
            assert k in self.splitnames, k + " used in connections but not in names"
        self.connections = self.create_connections()
        # Additional configuration based on complexity
        if self.ExitFractions is not None:
            if len(self.ExitFractions.keys()) == 0:
                self.ExitFractions = None
        if self.ExitFractions is not None:
            self.ExitFractions = EFractions(self.ExitFractions)
            self.inlineProbArgs = self.ExitFractions.inlineargs()



    def processcon(self, origin, dest):
        """
        Process individual connections.
        Creates connection code such as ``T.Exit( 'I^B', 'R')``

        """
        if(len(dest.keys())) == 0:
            return('T.NoExit("'+origin+'")\n')
        if(len(dest.keys())) == 1:
            return('T.Exit( "'+origin+'" , "'+list(dest.keys())[0]+'" )\n')
        if(len(dest.keys())) == 2:
            v1 = list(dest.keys())[0]
            v2 = list(dest.keys())[1]
            # TODO: assert probabilities?
            # assert dest[v1]["prob"]+dest[v2]["prob"]==1, "Proportion "+origin+"->"+v1+"+"+v2+"="+dest[v1]["prob"]+dest[v2]["prob"]+"!=1"
            prob = dest[v1]["prob"]
            prob_symb=[dest[v1]["prob_symb"], dest[v2]["prob_symb"]]
            return('T.SplitExit( "'+origin+'" , "'+v1+'" , "'+v2+'" , '+prob
                   +' , '+str(prob_symb)+' )\n')
        raise RuntimeError('Unable to handle more than two exit connections')
    def create_connections(self):
        """ Create connection code for all connections.

        Example output::

                T.BaseVar('S')
                T.Exit( 'S', 'E')
                T.SplitExit( 'E', 'I^A', 'I^S', 1-f, prob_symb=['1-f','f'])
                T.Exit( 'I^A', 'R')

        """
        code='T.BaseVar("'+self.BaseVar+'")\n'
        for origin in self.Exits:
            code+=self.processcon(origin, self.Exits[origin])
        code=textwrap.indent(code, indentStr)
        return code

    def gencode(self, indent=""):
        """Generate code for a function such as :meth:`victoriaepi.ama.Model_ama` or :meth:`victoriaepi.seid.Model_SEID` depending on complexity."""
        with open(buildingblocks_dir+self.complexity+"_ModelMatrix.py","r") as f:
            fullcode= f.read()
        variables=codeconfig_getvars(fullcode)
        if len(variables)>0:
            variables.sort(key=lambda x: len(x), reverse=True)
            for va in variables:
                if eval("self."+va[1:]) is None:
                    warnings.warn("Inserting None for "+va, UserWarning)
                fullcode=fullcode.replace(va, str(eval("self."+va[1:])))
        compile(fullcode, "<ModelMatrix Code>\n"+fullcode, 'exec')
        firstline=fullcode.splitlines()[0]
        return textwrap.indent(fullcode,
                               indent,
                               lambda x: not firstline in x)
    def genMaskcode(self, indent=""):
        """Generate code for the masks to select variables from list of state variables, for all variables.

        Example::


            self.mask_S   = self.T.SelectMask('S')
            self.mask_E   = self.T.SelectMask('E')
            self.mask_IA  = self.T.SelectMask('I^A')
            self.mask_ISs = self.T.SelectMask('I^S', E_range='all', as_col_vec=True)
            self.mask_IBs = self.T.SelectMask('I^B', E_range='all', as_col_vec=True)
            self.mask_ISs_flat = self.T.SelectMask('I^S', E_range='all')


        """
        code=''
        for var in self.splitnames:
            code+="self.mask_"+prettifyvarname(var)+" = self.T.SelectMask('"+var+"')\n"
            code+="self.mask_"+prettifyvarname(var)+"s_flat = self.T.SelectMask('"+var+"', E_range='all')\n"
            #Adding 's and _flat for splits
            if len(self.Exits[var].keys())==2:
                code+="self.mask_"+prettifyvarname(var)+"s = self.T.SelectMask('"+var+"', E_range='all', as_col_vec=True)\n"
        compile(code, "<ModelMatrix Mask Code>", 'exec')
        firstline=code.splitlines()[0]
        return textwrap.indent(code,
                               indent,
                               lambda x: not firstline in x)
    def genSelfExitProbscode(self, indent=""):
        return("")#TODO: genSelfExitProbscode
    def __repr__(self):
        """Return this object as raw code representation for a function such as :meth:`victoriaepi.ama.Model_ama` or :meth:`victoriaepi.seid.Model_SEID` depending on complexity."""
        return self.gencode()



class rawfuncbody():
    """

    Process a raw string as a function body and tests it for valid Python code with ``compile()``.
    """

    def __init__(self, strfun, indent=indentStr):
        self.strfun=strfun
        self.indent=indent
        strfun0=re.sub(r"(\s*)(return( |$|\n))", r"\1Return0=0#\3", strfun)
        #strfun0=strfun.replace("return ","Return0 = ")
        #strfun0=strfun.replace("return\n","#return\n")
        compile(strfun0, "<\n"+strfun0+"\n>", 'exec')
    def __repr__(self):
        return textwrap.indent(self.strfun,
                               self.indent)


@dataclass
class functiondef ():
    """Defining functions from ordered dicts."""
    defn: str
    body: str
    def __post_init__(self):
        assert re.search( "[:;]", self.defn) is None
        assert re.search( "^def ", self.defn) is None
        compile(self.defn, self.defn,"exec")
        self.defn ="def "+self.defn +":"
        self.body = rawfuncbody(self.body, indent = indentStr)
        self.baseIndent = ""
        compile(self.__repr__(), "<"+self.__repr__()+">", 'exec')
    def __repr__(self):
        return textwrap.indent(self.defn+"\n"+str(self.body),
                               self.baseIndent)
    @classmethod
    def fromodic(cls, odic, baseIndent=""):
        """
        Initialize instance from dictionary
        Parameters
        """
        assert re.search( r"[^\s]", baseIndent) is None
        cins = dacite.from_dict(data_class = cls, data = odic)
        cins.baseIndent = baseIndent
        return cins

class functionCollection ():
    """Handles a dictionary consisting of several :any:`victoriaepi.configparsing.modelclass.functiondef`.
    """
    def __init__(self, methodict, baseIndent):
        assert re.search( r"[^\s]", baseIndent) is None
        self.baseIndent = baseIndent
        self.methodict = methodict.copy()
        for k in self.methodict.keys():
            self.methodict[k] = functiondef.fromodic(self.methodict[k], self.baseIndent)
    def __repr__(self):
        """Generates Python code for the function"""
        res=""
        for k in self.methodict.keys():
            res+=str(self.methodict[k])+"\n\n"
        #remove indentation from first line as it is expected to be provided externaly
        res = res.replace(self.baseIndent,"", 1)
        return (res)

class solvingMethods (functionCollection):
    """ Making sure basic solving functions are present. This class is an extension of
    :any:`victoriaepi.configparsing.modelclass.functionCollection` that ensures that
    ``rhs``, ``solve``, ``llikelihood``, ``lprior``, ``support``, and ``sim_init`` are properly defined.
    """
    def __init__(self, methodict, baseIndent):
        super().__init__(methodict, baseIndent)
        funcs = self.methodict.keys()
        basicfuncs=['rhs', 'solve', 'llikelihood', 'lprior', 'support', 'sim_init']
        for bf in basicfuncs:
            assert bf in funcs, bf +" method expected in solvingMethods"
        r = re.compile("solve_plain.*")
        assert len(list(filter(r.match, self.methodict.keys() )))>0,"At least one 'solve_plain' method expected in solvingMethods"


class plottingMethods (functionCollection):
    """Making sure basic plotting functions are present. This class is an extension of
    :any:`victoriaepi.configparsing.modelclass.functionCollection` that ensures that at least one ``PlotEvolution.*`` method is present."""
    def __init__(self, methodict, baseIndent):
        super().__init__(methodict, baseIndent)
        r = re.compile("PlotEvolution.*")
        assert len(list(filter(r.match, self.methodict.keys() )))>0,"At least one 'PlotEvolution' method is expected in solvingMethods"