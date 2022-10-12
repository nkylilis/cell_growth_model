#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 12:33:23 2022

@author: nicolaskylilis
"""

import numpy as np

def model(parms,state_var = []):
    """
    cell size model based on a logistic function:
                
                       L
     f(x) = ----------------------
            1 + exp( -k * (x- m ))
    
    L: maximum value of the function
    k: steepnesness of the S-curve
    m: x-value of the mid-point of S-curve
    
    x: intracellular molecules state variable

    Parameters
    ----------
    par : 
        TYPE:     python list
        DESCRIPTION. list of parameters for the cell model function
    state_var : 
        TYPE, optional: python liast
        DESCRIPTION. The default is []. list of state variables for the cell model function

    Returns
    -------
    size : TYPE
        DESCRIPTION.

    """
 
    m = parms["mass_unit"]# m parameter
    k = parms["mass_k"] # k parameter
    mass_min = parms["mass_min"]
    mass_infl = parms["mass_infl"]
    
    a = state_var[0] # x state variable
    
    if (-k*(a-m)) < 700: # prevents overflow at > np.exp(700)
        size = mass_min +(mass_infl / (1 + np.exp(-k*(a-m))))
    else:
        size = 4E9
    
    return size
        