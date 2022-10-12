#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:23:56 2022

@author: nicolaskylilis
"""

"""
-------------------------------------------------------------------------------
fwr rate constant : 
                                    kc * mod
                                        where mod is activating:
                                                  a
                                        mod = ----------
                                               a  + Kmet
                                               
                                               a:      species -  energy levels
                                               Kmet: parameter -  half maximal activation threshold


-------------------------------------------------------------------------------
reaction:

                        ic --> ec + rbs -a
 

-------------------------------------------------------------------------------
"""

def fwr(ic, kc, act):
    
    a = act[0]
    Kmet = act[1]
    mod = a/(Kmet + a)
    
    mod_act = kc * mod
    
    dy = ic * mod_act
    
    return dy