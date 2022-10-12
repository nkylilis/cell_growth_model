#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:38:29 2022

@author: nicolaskylilis
"""
"""
-------------------------------------------------------------------------------
fwr rate constant : 
                                    gmax * 1/len_px * mod
                                    
                                        where mod is activating:
                                                  a
                                        mod = ----------
                                               a  + Kgamma
                                               
                                               a:      species -  energy levels
                                               Kgamma: parameter -  half maximal activatio threshold

-------------------------------------------------------------------------------
reaction:
                        ec + a*len_px  -->  p + ribo
 

-------------------------------------------------------------------------------
"""

def fwr(ec, gnorm, act):
    
    a = act[0]
    Kgamma = act[1]
    mod = a /(Kgamma + a)
    
    dy = +ec*gnorm*mod
    
    return dy