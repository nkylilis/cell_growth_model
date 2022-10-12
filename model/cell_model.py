#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:04:52 2022

@author: nicolaskylilis
"""
# pacckages
import pandas as pd
from scipy import integrate
import numpy as np
# import cell processes modules
import processes.ic_form as ic_form
import processes.ec_form as ec_form
import processes.pro_syn as pro_syn
import processes.transcription as tx
import processes.degradation as dm
import processes.complex as cmplx
import processes.energy_prod as energy
import processes.cell_size as cell_size

def simulate(ns=[], s0=[], fpath_params = [], new_params = [] ):
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
    #### Load default cell parameters
    fpath_params = "cell_parameters.csv"
    parms = pd.read_csv(fpath_params, index_col=("Parameter"))["Value"]
    
    
    #### Modify default parameters
    
    ### modify cell parameters values
    # from filepath
    if len(fpath_params) != 0: 
        fpath_params = fpath_params
        parms = pd.read_csv(fpath_params, index_col=("Parameter"))["Value"]
    
    # set values
    if len(new_params) != 0:
        for p in new_params.index: parms[p] = new_params[p]
    else:pass
    
    # nutrient quality
    if type(s0) == int:
        parms["ns"] =ns
        
    # nutrient amount
    if type(s0) == int:
        parms["s0"] =s0
    
    #correctiions
    parms["thetar"] = parms["thetar"] * parms["thetax"]
    

        
    
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
    rmr_0 = 1000
    rmc_0 = 1000
    rmq_0 = 1000
    # proteins
    r_0 = 1000
    em_0 = 1000
    q_0 = 1000
    # ribosomes
    ribo_0 = 1E4
    # m_ribo_0 = 100
    init = np.array([ a_0,mr_0, mc_0, mq_0,icr_0,icc_0,icq_0, rmr_0, rmc_0, rmq_0, r_0, em_0, q_0,ribo_0])    
        
    
    #### run simulation
    # simulation timespan 
    t_span = np.array([0,1e6])
    # solve 
    sim = integrate.solve_ivp(ode_system, t_span, init ,  method='BDF', args=([parms]), atol = 1e-6, rtol =1e-3) #Default values are 1e-3 for rtol and 1e-6 for atol.
    

    a       = sim.y[0,-1]
    # print(["a: ",a])
    mr      = sim.y[1,-1]
    mc      = sim.y[2,-1]
    mq      = sim.y[3,-1]
    icr     = sim.y[4,-1]
    icc     = sim.y[5,-1]
    icq     = sim.y[6,-1]
    rmr     = sim.y[7,-1]
    rmc     = sim.y[8,-1]
    rmq     = sim.y[9,-1]
    r       = sim.y[10,-1]
    em      = sim.y[11,-1]
    q       = sim.y[12,-1]
    ribo    = sim.y[13,-1]
    

    
    ## cell size 
    state_var = [a]
    cell_mass = cell_size.model(parms,state_var)
    
    
    # growth rate @steady_state
    gamma = parms["gmax"]*a/(parms["Kgamma"] + a)
    ttrate= (rmr + rmc + rmq)*gamma
    lam= (ttrate/cell_mass)
    # protein mass in aa
    rp_mass = ((ribo + icr + icc +icq + rmr + rmc +rmq)*parms["lenRibo"]) + (r * parms["lenR"])
    em_mass = em * parms["lenC"]
    q_mass  = q * parms["lenO"]
    
    simulation = {}
    simulation["growth rate (hour-1)"]              = lam * 60
    simulation["mRNA ribosomal sector (ratio)"]     = (mr + icr) /(mr + mc + mq + icr + icc + icq)
    simulation["mRNA catabolic sector (ratio)"]     = (mc + icc) /(mr + mc + mq + icr + icc + icq)
    simulation["mRNA others sector (ratio)"]        = (mq + icq) /(mr + mc + mq + icr + icc + icq)
    simulation["proteome ribosomal sector (ratio)"] = (rp_mass) /(rp_mass + em_mass + q_mass)
    simulation["proteome catabolic sector (ratio)"] = (em_mass) /(rp_mass + em_mass + q_mass)
    simulation["proteome others sector (ratio)"]    = (q_mass)  /(rp_mass + em_mass + q_mass)
    simulation["cell mass (aa)"]                    = cell_mass #[ rp_mass + em_mass + q_mass ] 
    simulation["peptide chain elongation rate (aa/min)"] = gamma
    simulation["ribosomes (molecules)"]             = ribo + icr + icc +icq + rmr + rmc +rmq
    simulation["mRNA (molecules)"]                  = mr + mc + mq + icr + icc + icq
    simulation["ribo_mRNA ratio"]                   = (rmr + rmc +rmq) / (mr + mc + mq + icr + icc + icq)
    
    
    df = pd.DataFrame([simulation], index=[parms["ns"]])
    
    return sim, df


def ode_system(t,y, p):

    
    """
    cell growth model v0.9
    """
    
    s = {}
    s["a"]      = y[0]
    s["m_rib"]  = y[1]
    s["m_met"]  = y[2]
    s["m_hsk"]  = y[3]
    s["ic_rib"] = y[4]
    s["ic_met"] = y[5]
    s["ic_hsk"] = y[6]
    s["ec_rib"] = y[7]
    s["ec_met"] = y[8]
    s["ec_hsk"] = y[9]
    s["p_rib"]  = y[10]
    s["p_met"]  = y[11]
    s["p_hsk"]  = y[12]
    s["ribo"]   = y[13]
    
    
    

    
    dydt = np.zeros(14)
    #--------------------------------------------------------------------------    
    #### Energy, Cell Size & dilution
    # energy production
    prod = +energy.fwr(s["p_met"], p["Vmax"], act=[p["s0"],p["Km"]]) *p["ns"]
    # energy utilisationn
    ti_rate = +ec_form.fwr(s["ic_rib"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              +ec_form.fwr(s["ic_met"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              +ec_form.fwr(s["ic_hsk"], p["kc"], act=[s["a"], p["Kmet"]] )
    tt_rate = +p["lenR"] * pro_syn.fwr(s["ec_rib"], p["gmax"]/p["lenR"], act=[s["a"], p["Kgamma"]]) \
              +p["lenC"] * pro_syn.fwr(s["ec_met"], p["gmax"]/p["lenC"], act=[s["a"], p["Kgamma"]]) \
              +p["lenO"] * pro_syn.fwr(s["ec_hsk"], p["gmax"]/p["lenO"], act=[s["a"], p["Kgamma"]]) 
    # cell size
    M =  cell_size.model(p,state_var = [s["a"]])
    # Growth rate/dilution
    dilution = tt_rate/M
    
    # energy (a)
    dydt[0] = +prod \
              -tt_rate - ti_rate\
              -dilution*s["a"]
    
    #--------------------------------------------------------------------------    
    #### RBSs/mRNAs
    # m_rib 
    dydt[1] = +tx.fwr(1,p["wr"], act=[s["a"],p["thetar"]]) \
              -dm.fwr(s["m_rib"],p["dm"]) \
              -ic_form.fwr(s["ribo"], s["m_rib"], p["kb_ribo"], inh=[s["p_rib"],p["Krepr"]])  \
              +ic_form.rev(s["ic_rib"], p["ku"]) \
              +ec_form.fwr(s["ic_rib"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -dilution*s["m_rib"]
    # m_met 
    dydt[2] = +tx.fwr(1,p["wc"], act=[s["a"],p["thetax"]]) \
              -dm.fwr(s["m_met"],p["dm"]) \
              -ic_form.fwr(s["ribo"], s["m_met"], p["kb_cat"]) \
              +ic_form.rev(s["ic_met"], p["ku"]) \
              +ec_form.fwr(s["ic_met"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -dilution*s["m_met"]
    # m_hsk 
    dydt[3] = +tx.fwr(1,p["wq"], act=[s["a"],p["thetax"]], inh=[s["p_hsk"],p["Kq"],p["nq"]]) \
              -dm.fwr(s["m_hsk"],p["dm"]) \
              -ic_form.fwr(s["ribo"], s["m_hsk"], p["kb_other"]) \
              +ic_form.rev(s["ic_hsk"], p["ku"]) \
              +ec_form.fwr(s["ic_hsk"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -dilution*s["m_hsk"]
    #--------------------------------------------------------------------------   
    #### Initiation complexes
    # ic_rib
    dydt[4] = +ic_form.fwr(s["ribo"], s["m_rib"], p["kb_ribo"], inh=[s["p_rib"],p["Krepr"]])  \
              -ic_form.rev(s["ic_rib"], p["ku"]) \
              -ec_form.fwr(s["ic_rib"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -dilution*s["ic_rib"]
    # ic_met
    dydt[5] = +ic_form.fwr(s["ribo"], s["m_met"], p["kb_cat"]) \
              -ic_form.rev(s["ic_met"], p["ku"]) \
              -ec_form.fwr(s["ic_met"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -dilution*s["ic_met"]
    # ic_hsk
    dydt[6] = +ic_form.fwr(s["ribo"], s["m_hsk"], p["kb_other"]) \
              -ic_form.rev(s["ic_hsk"], p["ku"]) \
              -ec_form.fwr(s["ic_hsk"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -dilution*s["ic_hsk"]
    #--------------------------------------------------------------------------
    #### Elongation complexes
    # ec_rib
    dydt[7] = +ec_form.fwr(s["ic_rib"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -pro_syn.fwr(s["ec_rib"], p["gmax"]/p["lenR"], act=[s["a"], p["Kgamma"]]) \
              -dilution*s["ec_rib"]
    # ec_met
    dydt[8] = +ec_form.fwr(s["ic_met"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -pro_syn.fwr(s["ec_met"], p["gmax"]/p["lenC"], act=[s["a"], p["Kgamma"]]) \
              -dilution*s["ec_met"]
    # ec_hsk
    dydt[9] = +ec_form.fwr(s["ic_hsk"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -pro_syn.fwr(s["ec_hsk"], p["gmax"]/p["lenO"], act=[s["a"], p["Kgamma"]]) \
              -dilution*s["ec_hsk"]   
    #--------------------------------------------------------------------------              
    #### Proteins
    # p_rib
    dydt[10] = +pro_syn.fwr(s["ec_rib"], p["gmax"]/p["lenR"], act=[s["a"], p["Kgamma"]]) \
              -cmplx.fwr(s["p_rib"] * p["lenR"]/p["lenRibo"], p["k_form"]) * (p["lenRibo"]/p["lenR"]) \
              -dilution*s["p_rib"]
    # p_met
    dydt[11] = +pro_syn.fwr(s["ec_met"], p["gmax"]/p["lenC"], act=[s["a"], p["Kgamma"]]) \
              -dilution*s["p_met"]
    # p_hsk
    dydt[12] = +pro_syn.fwr(s["ec_hsk"], p["gmax"]/p["lenO"], act=[s["a"], p["Kgamma"]]) \
              -dilution*s["p_hsk"]
     #--------------------------------------------------------------------------         Note: free ribos goes very negative     
    #### ribosome complex 
    dydt[13]  = +cmplx.fwr(s["p_rib"] * p["lenR"]/p["lenRibo"], p["k_form"]) \
                -ic_form.fwr(s["ribo"], s["m_rib"], p["kb_ribo"], inh=[s["p_rib"],p["Krepr"]]) -ic_form.fwr(s["ribo"], s["m_met"], p["kb_cat"]) -ic_form.fwr(s["ribo"], s["m_hsk"], p["kb_other"]) \
                +ic_form.rev(s["ic_rib"], p["ku"]) +ic_form.rev(s["ic_met"], p["ku"]) +ic_form.rev(s["ic_hsk"], p["ku"]) \
                +pro_syn.fwr(s["ec_rib"], p["gmax"]/p["lenR"], act=[s["a"], p["Kgamma"]]) +pro_syn.fwr(s["ec_met"], p["gmax"]/p["lenC"], act=[s["a"], p["Kgamma"]]) +pro_syn.fwr(s["ec_hsk"], p["gmax"]/p["lenO"], act=[s["a"], p["Kgamma"]]) \
                -dilution*s["ribo"]
    #--------------------------------------------------------------------------    
    return dydt