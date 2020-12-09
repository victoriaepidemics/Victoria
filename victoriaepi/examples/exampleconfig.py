#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model Config Exmple.
"""


from victoriaepi.configparsing.config import ModelConfig


fullmodel=ModelConfig.fromfilename("SEIDconfig.hjson")
fullmodel.writecode()


fullmodel = ModelConfig.fromfilename("AMAconfig.hjson")
fullmodel.writecode()

fullmodel = ModelConfig.fromfilename("GRPconfig.hjson")
fullmodel.writecode()
