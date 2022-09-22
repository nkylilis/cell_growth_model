#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 14:33:06 2021

@author: nicolaskylilis


Description:
    - Cell growth model v08:
        - mass growth models:
            - Transcription
            - Translation Initiation Complex formation
            - Translation Elongation Complex maturation
            - Translation Protein Synthesis
            - mRNA degradation
        - cell size model:
            - phenomenological, intracellular energy dependent

"""

# packages
import pandas as pd
from scipy import integrate
import numpy as np
import cell_size



#%% simulation caller

def simulate( ns, s0=[], fpath_params = [], new_params = [] ):
    """
    solver caller. Includes species initial conditions for ode problem solution

    Parameters
    ----------
    particle : TYPE
        DESCRIPTION.
    ns : TYPE
        DESCRIPTION.

    Returns
    -------
    sim : TYPE
        DESCRIPTION.

    """
    
    ## initial conditions
    #energy
    a_0 = 5
    # mRNAs
    mr_0 = 1000
    mc_0 = 1000
    mq_0 = 1000
    # initiation complexes
    icr_0 = 100
    icc_0 = 100
    icq_0 = 100
    # ribocomplexes
    rmr_0 = 10000
    rmc_0 = 10000
    rmq_0 = 10000
    # proteins
    r_0 = 1E6
    em_0 = 1E6
    q_0 = 1E6
    # ribosomes
    ribo_0 = 1E4
    # m_ribo_0 = 100
    init = np.array([ a_0,mr_0, mc_0, mq_0,icr_0,icc_0,icq_0, rmr_0, rmc_0, rmq_0, r_0, em_0, q_0,ribo_0])
    
        
    #### modify default parameters
    
    ### modify cell parameters
    # filepath
    if len(fpath_params) != 0: fpath_params = fpath_params
    else: fpath_params = "cell_paramameters.csv"
    parms = pd.read_csv(fpath_params, index_col=("Parameter"))["Value"]
    # set values
    if len(new_params) != 0:
        for p in new_params.index: parms[p] = new_params[p]
    else:pass
    
    # nutrient amount
    if type(s0) == int:
        parms["s0"] =s0
        
        
    
    #### run simulation
    # simulation timespan 
    t_span = np.array([0,1e6])
    # solve 
    sim = integrate.solve_ivp(ode_system, t_span, init ,  method='BDF', args=(parms,ns), atol = 1e-2, rtol =1e-2) #Default values are 1e-3 for rtol and 1e-6 for atol.
    
    return sim



#%% ODE system definition

def ode_system(t, y, parms, ns):
    
    """
    Cell growth coarse grain model linking gene expression to growth rate phenotype

    Model features:
        - Simplified catabolism
        - Translation initiation processes
        - Ribosome assembly
        
        
    """
    
    # passed parms= [wr, wc, wq, thetar, thetax, Kq, Kgamma, Vmax, unit_mass, k_form, Krepr, wribo, mass_k]
    
    ## cell and reaction parameters
    #general
    lenR =      parms["lenR"]
    lenO =      parms["lenO"]
    lenC =      parms["lenC"]
    lenRibo=    parms["lenRibo"]
    dm=         parms["dm"]
    # catabolism
    s0=         parms["s0"]
    ns=         ns
    Vmax=       parms["Vmax"]
    Km=         parms["Km"]
    #transcription
    thetar=     parms["thetar"]*parms["thetax"]
    thetax=     parms["thetax"]
    wr=         parms["wr"]
    wc=         parms["wc"] 
    wq=         parms["wq"]
    Kq=         parms["Kq"]
    nq=         parms["nq"]
    # translation
    gmax=       parms["gmax"]
    Kgamma=     parms["Kgamma"]
    kc =        parms["kc"]
    # translation initiation {deltaG = 1180 #kcal/mol}
    ku=         parms["ku"]
    kb_ribo=    parms["kb_ribo"]
    kb_cat=     parms["kb_cat"]
    kb_others=  parms["kb_other"]
    # ribosome assembly
    k_form =    parms["k_form"]
    Krepr =     parms["K_repr"]
    # wribo =     parms["wribo"]

    

    ## cell species	
    # energy	
    a =         y[0] 	
    # mrnas	
    mr =        y[1]	
    mc =        y[2]	
    mq =        y[3]	
    # initiation complexes	
    icr =       y[4]	
    icc =       y[5]	
    icq =       y[6]	
    # ribocomplexes	
    rmr=        y[7]	
    rmc=        y[8]	
    rmq=        y[9]	
    # proteins	
    r=          y[10]	
    em=         y[11]	
    q=          y[12]	
    # ribosomes	
    ribo=       y[13]
    # m_ribo =    y[14]
    	
    	
    ### dynamic parameters	
    
    ## cell size
    state_var = [a]
    M = cell_size.model(parms,state_var)
    
    
    ## mass growth
    # catabolism	
    nucat= em * ((Vmax*s0)/(Km + s0))	
    # energy dependent process	
    gamma   = gmax * (a/(Kgamma + a))
    kc_eff  = kc   * (a/(Kgamma + a))		
    # translation effective rate constants
    tirate  = (icr + icc + icq) * kc_eff
    ttrate  = (rmr + rmc + rmq) * gamma
    # riboprotein repression	
    mod_rp = 1/(1 + r/Krepr)	
    # number of ribosomal proteins sets
    num_sets = r/(lenRibo/lenR)
    ## cell growth rate
    lam= (ttrate/M)	

	
    # ode system	
    dydt = np.zeros(14)	
    
    # energy rxns	
    dydt[0]= +ns*nucat -tirate -ttrate -lam*a;	
    # mRNAs rxns	
    dydt[1]= +(wr*a/(thetar + a))          -kb_ribo*ribo*mr*mod_rp  +ku*icr   +kc_eff*icr   -dm*mr -lam*mr	
    dydt[2]= +(wc*a/(thetax + a))          -kb_cat*ribo*mc          +ku*icc   +kc_eff*icc   -dm*mc -lam*mc	
    if (q/Kq) > 10:
        mod_fcn= 0
    elif (q/Kq) < 0.01:
        mod_fcn= 1
    else:
        mod_fcn = (1/(1 + (q/Kq)**nq))  
    dydt[3]= +(wq*a/(thetax + a)*mod_fcn)  -kb_others*ribo*mq       +ku*icq   +kc_eff*icq   -dm*mq -lam*mq	
    
    # initiation complexes	
    dydt[4]= +kb_ribo  *ribo *mr *mod_rp   -ku*icr    -kc_eff*icr   -lam*icr;	
    dydt[5]= +kb_cat   *ribo *mc           -ku*icc    -kc_eff*icc   -lam*icc;	
    dydt[6]= +kb_others*ribo *mq           -ku*icq    -kc_eff*icq   -lam*icq;	
    	
    # ribocomplexes rxns	
    dydt[7]= +kc_eff*icr -gamma/lenR*rmr -lam*rmr	
    dydt[8]= +kc_eff*icc -gamma/lenC*rmc -lam*rmc	
    dydt[9]= +kc_eff*icq -gamma/lenO*rmq -lam*rmq 	
    	
    # proteins processes:
    # ------------------
    # 1) translation elongation + termination
    # 2) dilution
    # *) for r proteins: ribosome assembly
    dydt[10]= +gamma/lenR*rmr           -lam*r          -k_form*num_sets*(lenRibo/lenR)
    dydt[11]= +gamma/lenC*rmc           -lam*em	
    dydt[12]= +gamma/lenO*rmq           -lam*q	
    	
    
    # ribosome processes:
    # -------------------
    # 1) ribosome assembbly 
    # 2) formation of initiation complexes 
    # 3) dissolution of initiation complexes 
    # 4) completion of protein translation 
    # 5) dilution
    dydt[13] = +k_form*num_sets                                                 \
               -kb_ribo*ribo*mr*mod_rp   -kb_cat*ribo*mc     -kb_others*ribo*mq \
               +ku*icr                   +ku*icc             +ku*icq            \
               +gamma/lenR*rmr           +gamma/lenC*rmc     +gamma/lenO*rmq    \
               -lam*ribo
    
    
    # # m_ribo
    # dydt[14] = (wribo*a/(thetar + a)) -(k_form*r*(lenR/lenRibo)*m_ribo) -lam*m_ribo
               
    
    return dydt


