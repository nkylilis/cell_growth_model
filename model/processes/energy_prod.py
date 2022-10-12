#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 15:55:10 2022

@author: nicolaskylilis

-------------------------------------------------------------------------------
fwr rate constant : 
                            mod = Vmax*mod
                                where mod is activating:
                                          s0
                                mod = ----------
                                       Km  + s0
                                       
                                       s0:   parameter -  amount of nutrients extracellularly
                                       Km:  parameter -  half maximal activation threshold for nutrient utilisation


-------------------------------------------------------------------------------
reaction:
                         p_cat  --> a + p_cat
 

-------------------------------------------------------------------------------
"""


def fwr(p_cat, Vmax, act):
    
    
    s0 = act[0]
    Km = act[1]
    mod = s0/(Km + s0)
    
    dy  = + Vmax * p_cat * mod
    
    
    return dy




