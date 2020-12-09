#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests
"""
import runpy

from victoriaepi.configparsing.config import ModelConfig

def func(conffile):
    fullmodel=ModelConfig.fromfilename(conffile)
    fullmodel.writecode()
    return(fullmodel)


def test_SEID():
    assert isinstance(func("SEIDconfig.hjson"),ModelConfig)

def test_AMA():
    assert isinstance(func("AMAconfig.hjson"),ModelConfig)

def test_GRP():
    assert isinstance(func("GRPconfig.hjson"),ModelConfig)

def test_SEIDclassrun():
   fglobals = runpy.run_path("autoseid.py")
   #exec(open("autoseid.py").read())

def test_AMAclassrun():
   fglobals = runpy.run_path("autoAMA.py")
   #exec(open("autoAMA.py").read())

def test_GRPclassrun():
   fglobals = runpy.run_path("autoGRP.py")
   #exec(open("autoGRP.py").read())

def test_SEIDcallrun():
   fglobals = runpy.run_path("call_autoseid.py")
