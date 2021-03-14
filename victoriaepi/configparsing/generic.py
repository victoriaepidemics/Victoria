#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions.
"""

import re

#%%Generic functions
def codeconfig_getvars(code):
    """
    Extract configured parameters in template code.

    Parameters
    ----------
    code : str
        Raw template code.
    """
    res = re.findall(r'\$\w+', code)
    variables=list(set(res))
    variables.sort(key=lambda x: len(x),reverse=True)
    return(variables)

def prettifyvarname(var):
    """ Remove invalid characters"""
    return(re.sub (r"[^\w]","",var))