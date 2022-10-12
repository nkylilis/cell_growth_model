#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 13:00:20 2022

@author: nicolaskylilis
"""



#%% define system
    
def ode_system(t,y, parms):
    # import cell processes modules
    import processes.ic_form as ic_form
    import processes.ec_form as ec_form
    import processes.pro_syn as pro_syn
    
    """
    traslation model
    
    """
    
    s = {}
    s["rbs"] = y[0]
    s["ribo"]  = y[1]
    s["ic"]   = y[2]
    s["ec"]   = y[3]
    s["a"]    = y[4]
    s["p"]    = y[5]
    
    p = parms
    
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
    dydt[4] = 0
    
    # protein
    dydt[5] = +pro_syn.fwr(s["ec"], p["gmax"]/p["lenX"], act=[s["a"], p["Kgamma"]])
    
    
    return dydt

#%% solve - translation
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

rbs_0 = 50
ribo_0 = 1000
ic_0 = 0
ec_0 = 0
a_0 = 10
p_0 = 1

init = [rbs_0,ribo_0, ic_0,ec_0, a_0, p_0]

# model parameters
p = {}
p["kon"]    = 1    # min-1*mol-1
p["koff"]   = 0.1  # min-1
p["Kgamma"] = 10   # molecules
p["Kmet"]   = p["Kgamma"]*1.1
p["kc"]     = 120  # min-1
p["gmax"]   = 1320 # aa/min
p["lenX"]   = 100  # codons
p["fp"]     = 10   # codons


t_span = np.array([0, 1e6])
sim = integrate.solve_ivp(ode_system, t_span, init, args=([p]), method="BDF")

    


#plot

t = sim.t
y = sim.y 

rbs = y[0,:]
ribo = y[1,:]
ic = y[2,:]
ec = y[3,:]
a  = y[4,:]
pro  = y[5,:]
pos = p["lenX"]/p["fp"] * rbs_0

plt.semilogx(t, rbs)
plt.semilogx(t, ribo)
plt.semilogx(t, ic)
plt.semilogx(t, ec)
#plt.semilogx(t, pro)
plt.semilogx([1E-10,1E6],[pos,pos], linestyle="--")
# plt.ylim([1E-3,10000])
plt.ylabel("elong. complexes (molecules)")
# plt.xlim([1E-3,1E6])
plt.xlabel("simulation time")
plt.legend(["rbs","ribo","ic","ec", "ribopositions"], title="system species")
plt.title("Translation model - simulation")
plt.show()

#%% solve - energy var
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

rbs_0 = 50
ribo = 1000
ic_0 = 0
ec_0 = 0
a_0 = 10
p_0 = 1

init = [rbs_0,ribo_0, ic_0,ec_0, a_0, p_0]

# model parameters
p = {}
p["kon"]    = 1    # min-1*mol-1
p["koff"]   = 0.1  # min-1
p["Kgamma"] = 10   # molecules
p["Kmet"]   = p["Kgamma"]*1.1
p["kc"]     = 120  # min-1
p["gmax"]   = 1320 # aa/min
p["lenX"]   = 100  # codons
p["fp"]     = 10   # codons


sim_lst = []
t_span = np.array([0, 1e6])
a_lst = [1,10,100, 1000]
for a in a_lst:
    init[4] = a
    
    sim = integrate.solve_ivp(ode_system, t_span, init, args=([p]), method="BDF")
    sim_lst += [sim]
    


#plot

for sim in sim_lst:

    t = sim.t
    y = sim.y 
    
    rbs = y[0,:]
    ribo = y[1,:]
    ic = y[2,:]
    ec = y[3,:]
    a  = y[4,:]
    pro  = y[5,:]
    pos = p["lenX"]/p["fp"] * rbs_0
    
    # plt.semilogx(t, rbs)
    #plt.semilogx(t, ribo)
    #plt.semilogx(t, ic)
    plt.semilogx(t, ec)
    #plt.semilogx(t, pro)
    # plt.legend(["rbs","ribo","ic","ec"])

plt.semilogx([1E-10,1E6],[pos,pos], linestyle="--")
# plt.ylim([1E-3,10000])
plt.ylabel("elong. complexes (molecules)")
# plt.xlim([1E-3,1E6])
plt.xlabel("simulation time")
plt.legend(a_lst + ["ribopositions"], title="energy level")
plt.title("Parameter scan - energy")
plt.show()


#%% solve - ribo var
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

rbs_0 = 50
ribo = 1000
ic_0 = 0
ec_0 = 0
a_0 = 10
p_0 = 1

init = [rbs_0,ribo_0, ic_0,ec_0, a_0, p_0]

# model parameters
p = {}
p["kon"]    = 1    # min-1*mol-1
p["koff"]   = 0.1  # min-1
p["Kgamma"] = 10   # molecules
p["Kmet"]   = p["Kgamma"]*1.1
p["kc"]     = 120  # min-1
p["gmax"]   = 1320 # aa/min
p["lenX"]   = 100  # codons
p["fp"]     = 10   # codons

sim_lst = []
t_span = np.array([0, 1e6])
ribo_lst = [10,100,1000, 10000]
for ribo in ribo_lst:
    init[1] = ribo
    
    sim = integrate.solve_ivp(ode_system, t_span, init, args=([p]), method="BDF")
    sim_lst += [sim]
    


#plot

for sim in sim_lst:

    t = sim.t
    y = sim.y 
    
    rbs = y[0,:]
    ribo = y[1,:]
    ic = y[2,:]
    ec = y[3,:]
    a  = y[4,:]
    pro  = y[5,:]
    pos = p["lenX"]/p["fp"] * rbs_0
    
    # plt.semilogx(t, rbs)
    # plt.semilogx(t, ribo)
    # plt.semilogx(t, ic)
    plt.semilogx(t, ec)
    #plt.semilogx(t, pro)
    # plt.legend(["rbs","ribo","ic","ec"])

plt.semilogx([1E-10,1E6],[pos,pos], linestyle="--")
# plt.ylim([1E-3,10000])
plt.ylabel("elong. complexes (molecules)")
# plt.xlim([1E-3,1E6])
plt.xlabel("simulation time")
plt.legend(ribo_lst)
plt.title("Parameter scan - ribosomes")
plt.show()

#%% solve - CDS var
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

rbs_0 = 50
ribo = 1000
ic_0 = 0
ec_0 = 0
a_0 = 10
p_0 = 1

init = [rbs_0,ribo_0, ic_0,ec_0, a_0, p_0]

# model parameters
p = {}
p["kon"]    = 1    # min-1*mol-1
p["koff"]   = 0.1  # min-1
p["Kgamma"] = 10   # molecules
p["Kmet"]   = p["Kgamma"]*1.1
p["kc"]     = 120  # min-1
p["gmax"]   = 1320 # aa/min
p["lenX"]   = 100  # codons
p["fp"]     = 10   # codons

sim_lst = []
t_span = np.array([0, 1e6])
cds_lst = [50,250,500, 1000]
for l in cds_lst:
    p["lenX"] = l
    
    sim = integrate.solve_ivp(ode_system, t_span, init, args=([p]), method="BDF")
    sim_lst += [sim]
    


#plot

for sim in sim_lst:

    t = sim.t
    y = sim.y 
    
    rbs = y[0,:]
    ribo = y[1,:]
    ic = y[2,:]
    ec = y[3,:]
    a  = y[4,:]
    pro  = y[5,:]
    pos = p["lenX"]/p["fp"] * rbs_0
    
    # plt.semilogx(t, rbs)
    # plt.semilogx(t, ribo)
    # plt.semilogx(t, ic)
    plt.semilogx(t, ec)
    #plt.semilogx(t, pro)
    # plt.legend(["rbs","ribo","ic","ec"])

plt.semilogx([1E-10,1E6],[pos,pos], linestyle="--")
# plt.ylim([1E-3,10000])
plt.ylabel("elong. complexes (molecules)")
# plt.xlim([1E-3,1E6])
plt.xlabel("simulation time")
plt.legend(cds_lst + ["ribopositions"],title="# of codons")
plt.title("Parameter scan - CDS length")
plt.show()



#%% solve - Kmet

from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

rbs_0 = 50
ribo = 1000
ic_0 = 0
ec_0 = 0
a_0 = 10
p_0 = 1

init = [rbs_0,ribo_0, ic_0,ec_0, a_0, p_0]

# model parameters
p = {}
p["kon"]    = 1    # min-1*mol-1
p["koff"]   = 0.1  # min-1
p["Kgamma"] = 10   # molecules
p["Kmet"]   = p["Kgamma"]*1.1
p["kc"]     = 120  # min-1
p["gmax"]   = 1320 # aa/min
p["lenX"]   = 100  # codons
p["fp"]     = 10   # codons


t_span = np.array([0, 1e6])
kmet_lst = [1,1.1,1.5,2]
for k in kmet_lst:
    p["Kmet"] = p["Kgamma"]*k
    
    sim_lst = []
    for a in a_lst: 
        init[4] = a
    
        sim = integrate.solve_ivp(ode_system, t_span, init, args=([p]), method="BDF")
        sim_lst += [sim]
    


    #plot
    a_ =[]
    ec_ = []
    rbs_ = []
    ic_ =[]
    r_ = []
    for sim in sim_lst:
    
        t = sim.t
        y = sim.y 
        
        rbs = y[0,:]
        ribo = y[1,:]
        ic = y[2,:]
        ec = y[3,:]
        a  = y[4,:]
        pro  = y[5,:]
        pos = p["lenX"]/p["fp"] * rbs_0
        
        a_ +=[a[-1]]
        ec_ +=[ec[-1]]
        rbs_ +=[rbs[-1]]
        ic_ +=[ic[-1]]
        
        r_ += [ec[-1]/(rbs[-1]+ic_[-1])]
        
        # plt.semilogx(t, rbs)
        # plt.semilogx(t, ribo)
        # plt.semilogx(t, ic)
        
        #plt.semilogx(t, pro)
        # plt.legend(["rbs","ribo","ic","ec"])
    plt.semilogx(a_,r_)
plt.legend(kmet_lst,title="Kmet (relative to Kgamma)")
#plt.semilogx([1E-10,1E6],[pos,pos], linestyle="--")
# plt.ylim([1E-3,10000])
plt.ylabel("elong. complexes / mRNA molecules (ec/(rbs+ic))")
# plt.xlim([1E-3,1E6])
plt.xlabel("energy levels")

plt.title("Parameter scan - Kmet value")
plt.show()




