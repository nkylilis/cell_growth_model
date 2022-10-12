#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 12:48:02 2022

@author: nicolaskylilis
"""

def fwr(gene_copy,w,act,inh=1):

    """
    -------------------------------------------------------------------------------
    fwr rate constant : 
                                        w * mod_act * mod_inh
                                        
                                        where:
                                                  a
                                        mod_act = ----------
                                               theta + a
                                               
                                               a:      species -  energy levels
                                               theta:  parameter -  half maximal activation threshold
                                               
                                               
                                                     1
                                        mod_act = ----------
                                                 1 +  (p_hsk/Krepre)
                                               
                                               a:      species -  energy levels
                                               theta:  parameter -  half maximal activation threshold


    -------------------------------------------------------------------------------
    reaction:
                            gene_copy  --> mRNA
     

    -------------------------------------------------------------------------------
    """
    
    #### modulators
    
    # Energy dependance
    a = act[0]
    theta = act[1]
    mod_act = a/(theta + a)


    # hsk class inhibition
    if inh != 1:
        q =     inh[0]
        Kq = inh[1]
        nq =    inh[2]
    
        if (q/Kq) > 10:
            mod_inh= 0
        elif (q/Kq) < 0.01:
            mod_inh= 1
        else:
            mod_inh = (1/(1 + (q/Kq)**nq))  
    else:
        mod_inh =1
    
    
    #### flux
    dy = w* gene_copy * mod_act * mod_inh
    
    return dy
    