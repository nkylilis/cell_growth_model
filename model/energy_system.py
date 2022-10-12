#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 16:02:54 2022

@author: nicolaskylilis
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import processes.ic_form as ic_form
import processes.ec_form as ec_form
import processes.pro_syn as pro_syn
import processes.energy_prod as energy_prod

def ode_system(t,y, parms):
    

    
    s = {}
    s["rbs"] = y[0]
    s["ribo"]  = y[1]
    s["ic"]   = y[2]
    s["ec"]   = y[3]
    s["a"]    = y[4]
    s["p"]    = y[5]
    
    p = parms
    
    
    ### ode system
    dydt = np.zeros(6)
    # rbs
    dydt[0] = -ic_form.fwr(s["ribo"], s["rbs"], p["kon"]) \
              +ic_form.rev(s["ic"], p["koff"]) \
              +ec_form.fwr(s["ic"], p["kc"], act=[s["a"], p["Kmet"]] )
    
    # ribo
    dydt[1] = -ic_form.fwr(s["ribo"], s["rbs"], p["kon"]) \
              +ic_form.rev(s["ic"], p["koff"]) \
              +pro_syn.fwr(s["ec"], p["gmax"]/p["lenX"], act=[s["a"], p["Kgamma"]])
    
    # ic
    dydt[2] = +ic_form.fwr(s["ribo"], s["rbs"], p["kon"]) \
              -ic_form.rev(s["ic"], p["koff"]) \
              -ec_form.fwr(s["ic"], p["kc"], act=[s["a"], p["Kmet"]] )
    
    # ec
    dydt[3] = +ec_form.fwr(s["ic"], p["kc"], act=[s["a"], p["Kmet"]] ) \
              -pro_syn.fwr(s["ec"], p["gmax"]/p["lenX"], act=[s["a"], p["Kgamma"]])
    
    
    # a
    a_prod  = +p["ns"] *    energy_prod.fwr(s["p"], p["Vmax"], act=[p["s0"],p["Km"]])
    ti_rate = +             ec_form.fwr(s["ic"], p["kc"], act=[s["a"], p["Kmet"]] )
    tt_rate = +p["lenX"]*   pro_syn.fwr(s["ec"], p["gmax"]/p["lenX"], act=[s["a"], p["Kgamma"]])
    dydt[4] = +a_prod - ti_rate - tt_rate
    
    # protein
    dydt[5] = 0
    # dydt[5] = +pro_syn.fwr(s["ec"], p["gmax"]/p["lenX"], act=[s["a"], p["Kgamma"]])
    
    
    return dydt


## initial values for state variables 
rbs_0 = 1000
ribo_0 = 10000
ic_0 = 0
ec_0 = 0
a_0 = 0
p_0 = 10

init = [rbs_0,ribo_0, ic_0,ec_0, a_0, p_0]

# model parameters
p = {}
p["kon"]    = 1    # min-1*mol-1
p["koff"]   = 0.1  # min-1
p["Kgamma"] = 10   # molecules
p["Kmet"]   = p["Kgamma"]*1.1
p["kc"]     = 120  # min-1
p["gmax"]   = 1320 # aa/min
p["lenX"]   = 50  # codons
p["fp"]     = 10   # codons
p["Vmax"]   = 1
p["s0"]     = 1E4
p["Km"]     = 1E2
p["ns"]     = 5

mp_lst =[]
a_lst=[]
prod_lst =[]
ti_lst = []
tt_lst = []

t_span = np.array([0, 1e6])
s0_list = np.geomspace(10,1000,10)
for s0 in s0_list:
    p["s0"]     = s0
        
    sim = integrate.solve_ivp(ode_system, t_span, init, args=([p]), method="BDF")
    
    
    #%% plot 
    
    t = sim.t
    y= sim.y
    
    rbs = y[0,:]
    ribo = y[1,:]
    ic = y[2,:]
    ec = y[3,:]
    a  = y[4,:]
    pro  = y[5,:]
    
    mp_lst += [ p["lenX"]*   pro_syn.fwr(ec[-1], p["gmax"]/p["lenX"], act=[a[-1], p["Kgamma"]])]
    
    prod_lst+=[+ p["ns"] *    energy_prod.fwr(pro[-1], p["Vmax"], act=[p["s0"],p["Km"]]) ]
    ti_lst += [-            ec_form.fwr(ic[-1], p["kc"], act=[a[-1], p["Kmet"]] ) ]
    tt_lst += [-p["lenX"]*   pro_syn.fwr(ec[-1], p["gmax"]/p["lenX"], act=[a[-1], p["Kgamma"]])]
    a_lst  += [a[-1]]
    
    
    # plt.plot(t, rbs)
    # plt.plot(t, ribo)
    # plt.plot(t, ic)
    # plt.plot(t, ec)
    # plt.plot(t, pro)
    # plt.plot(t, a)
    # plt.legend(["rbs","ribo","ic","ec", "a"])
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()
    
plt.plot(s0_list,mp_lst, color="black", linestyle="--")
plt.plot(s0_list,prod_lst, color="green")
plt.plot(s0_list,ti_lst, color="orangered")
plt.plot(s0_list,tt_lst, color="red")

labels = ["Bioproduction rate (amino acids)", "Energy production rate", "EC formation energy consumption rate", "Protein Synthesis energy consumption rate"]
plt.legend(labels, fontsize=10)
plt.xlabel("extracellular nutrients amount")
plt.ylabel("rate (molecules per unit of time)")
# plt.xscale("log")
plt.show()

# # energy levels
# plt.plot(s0_list,a_lst, color="yellow")
# plt.show()
