#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:23:13 2022

@author: nicolaskylilis
"""

"""
-------------------------------------------------------------------------------
fwr rate constant : 
                                    kon


-------------------------------------------------------------------------------
reaction:
                        ribo + rbs  <--> ic
 

-------------------------------------------------------------------------------
rev rate constant : 
                                    koff


-------------------------------------------------------------------------------
"""


def fwr(ribo, rbs, kon, inh =[]):
    
    
    #### Modulation
    
    # ribosomal proteins innhibition
    if inh !=[]:
        p_rib = inh[0]
        Krepr = inh[1]
        mod_inh = 1/( 1 + (p_rib/Krepr))
    else: mod_inh =1
    
    
    #### Flux
    dy  = + kon * ribo * rbs * mod_inh
    
    
    return dy




def rev(ic, koff):
    
    
    dy  = +koff * ic
    
    return dy