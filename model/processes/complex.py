#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 13:37:46 2022

@author: nicolaskylilis
"""
def fwr(species, rate):
         
    """
    -------------------------------------------------------------------------------
    fwr rate constant : 
                                        k_form (rate constant)
    
    
    -------------------------------------------------------------------------------
    reaction:
                             p_rib (species)  -->  ribosome
     
    
    -------------------------------------------------------------------------------
    """
    
    
    #### Flux
    print(species,rate)
    dy= species * rate
    
    return dy