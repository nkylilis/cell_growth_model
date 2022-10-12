#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:02:44 2022

@author: nicolaskylilis


"""



# packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


#### figure
fig = plt.figure(constrained_layout=True,figsize=(16,22*0.8))
axd = fig.subplot_mosaic(
    """
    Aab
    BBC
    dDE
    FGH
    """)
    
#%% Panel A  
axd["a"].axis("off")
axd["b"].axis("off")
axd["A"].axis("off")
axd["A"].text(-0.15, 1.05, 'A',fontsize=20, fontweight=1000, transform=axd["A"].transAxes)
img = mpimg.imread('F1A_cell_model.png')
imagebox = OffsetImage(img, zoom=0.3)
xy = [0,0]
ab = AnnotationBbox(imagebox, xy, frameon=(False), box_alignment=(0, 0))
axd["A"].add_artist(ab)


  
#%% Panel B

axd["B"].axis("off")
axd["B"].text(-0.1, 1.05, 'B',fontsize=20, fontweight=1000, transform=axd["B"].transAxes)
img = mpimg.imread('F1B_translation_model.png')
imagebox = OffsetImage(img, zoom=0.55)
xy = [0,0]
ab = AnnotationBbox(imagebox, xy, frameon=(False), box_alignment=(0.05, 0))
axd["B"].add_artist(ab)



#%% Panels C-E - translattion model v09

#define system
    
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
  
## initial values for state variables 
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

y= sim.y

rbs = y[0,:]
ribo = y[1,:]
ic = y[2,:]
ec = y[3,:]
a  = y[4,:]
pro  = y[5,:]
pos = p["lenX"]/p["fp"] * rbs_0


#### - Panel C
axd["C"].plot(sim.t,rbs, color="teal")
axd["C"].plot(sim.t,ribo, color="red")
axd["C"].plot(sim.t,ic, color="purple")
axd["C"].plot(sim.t,ec, color="darkorange")
axd["C"].plot([1E-10,1E6],[pos,pos], linestyle="--", color = "black")

axd["C"].set_xscale("log")
#ax1.set_yscale("log")

labels = ["RBS/mRNAs",  "Ribosomes", "Init.Complexes", "Elong.Complexes", "ribopositions"]
axd["C"].legend(labels,fontsize="large", title="System Species")
axd["C"].set_xlabel("simulation time", size=14)
axd["C"].set_ylabel("molecules", size=14)
axd["C"].set_title('Translation model simulation', loc='center', size=14)
panel="C"
axd[panel].text(-0.1, 1.05, "C" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)


#### - Panel d enrgy levels scan


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
    
    axd["d"].semilogx(t, ec)


axd["d"].semilogx([1E-10,1E6],[pos,pos], linestyle="--", color="black")
axd["d"].set_ylabel("elongating complexes (molecules)", size=14)
axd["d"].set_xlabel("simulation time", size=14)
axd["d"].legend(a_lst + ["ribopositions"], title="system energy level")
axd["d"].set_title("Parameter scan - energy", loc='center', size=14)
axd["d"].set_xlim([1E-5,1E6])
panel="d"
axd[panel].text(-0.1, 1.05, "D" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)




#### ribosomes number variation

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
    axd["D"].semilogx(t, ec)

axd["D"].semilogx([1E-10,1E6],[pos,pos], linestyle="--", color="black")
axd["D"].set_ylabel("elongating complexes (molecules)", size=14)
axd["D"].set_xlabel("simulation time", size=14)
axd["D"].legend(ribo_lst, title="system ribosomes")
axd["D"].set_title("Parameter scan - ribosomes ", loc='center', size=14)
axd["D"].set_xlim([1E-5,1E6])
panel="D"
axd[panel].text(-0.1, 1.05, "E" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)


#### protein size variation
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

rbs_0 = 50
ribo = 100
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
c_lst = ["blue","orange","green","red"]
for sim,c,l in zip(sim_lst,c_lst,cds_lst):

    t = sim.t
    y = sim.y 
    
    rbs = y[0,:]
    ribo = y[1,:]
    ic = y[2,:]
    ec = y[3,:]
    a  = y[4,:]
    pro  = y[5,:]
    pos = l/p["fp"] * rbs_0
    axd["E"].semilogx(t, ec, color=c)
    axd["E"].semilogx([1E-10,1E6],[pos,pos], linestyle="--", color=c)


axd["E"].set_ylabel("elongating complexes (molecules)", size=14)
axd["E"].set_xlabel("simulation time", size=14)
axd["E"].legend([str(cds_lst[0])]+ ["ribopositions"] + [str(cds_lst[1])] + ["ribopositions"] + [str(cds_lst[2])]+ ["ribopositions"] + [str(cds_lst[3])] + ["ribopositions"],title="# of codons")
axd["E"].set_title("Parameter scan - CDS size (aa)", size=14, loc="center")
panel="E"
axd[panel].text(-0.1, 1.05, "F" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)


#### kmet

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
kmet_lst = [1,1.25,1.5,1.75,2]
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
    axd["F"].semilogx(a_,r_)
axd["F"].legend(kmet_lst,title="Kmet (relative to Kgamma)")
axd["F"].set_ylabel("EC / (rbs+ic)", size=14)
axd["F"].set_xlabel("system energy levels", size=14)
axd["F"].set_title("Ennergy funnction modulation", size=14, loc="center")
panel="F"
axd[panel].text(-0.1, 1.05, "G" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)

#%% Panel G - energy system

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
    
    
    #plot 
    
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
    
axd["G"].plot(s0_list,mp_lst, color="black", linestyle="--")
axd["G"].plot(s0_list,prod_lst, color="green")
axd["G"].plot(s0_list,ti_lst, color="orangered")
axd["G"].plot(s0_list,tt_lst, color="red")
axd["G"].set_title("Energy production/consumption model", loc='center', size=14)
labels = ["Bioproduction rate (amino acids)", "Energy production rate", "EC formation energy consumption rate", "Protein Synthesis energy consumption rate"]
axd["G"].legend(labels, fontsize=10)
axd["G"].set_xlabel("extracellular nutrients amount",size=14)
axd["G"].set_ylabel("rate (molecules per unit of time)",size=14)
panel="G"
axd[panel].text(-0.1, 1.05, "H" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)
#%% Panel H  cell size model

#### cell size data

data_size = [0.76E9, 1.19E9, 1.64E9, 2.15E9, 2.4E9]
data_gr = [0.6,1,1.5,2,2.5]

axd["H"].scatter(data_gr,data_size, color = "r")
axd["H"].set_xlabel("growth rate (h-1)", size=14, color="r")
axd["H"].set_ylabel("cell size (amino acids)", size=14, color="red")
limy =[1E8,4E9]
limx =[0.2,5]
axd["H"].set_ylim(limy)
axd["H"].set_xlim(limx)
axd["H"].tick_params(colors='red')
axd["H"].legend(["experimental data"],loc= 'lower right', fontsize=14)
panel="H"
axd[panel].text(-0.1, 1.05, "G" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)


#### celll size model

import processes.cell_size as cell_size

size_lst = []
a_lst = np.arange(0,100,1)
parms = {}
parms["mass_unit"]=35# m parameter
parms["mass_k"]=0.06 # k parameter
parms["mass_min"] =0.25E9
parms["mass_infl"] =3.25E9

for a in a_lst:
    y = cell_size.model(parms,state_var = [a])
    size_lst +=[y]


limx2 =[0,100]
ax = fig.add_subplot(4,3,12)
ax.plot(a_lst,size_lst, color = "k")
ax.patch.set_alpha(0)



#ax.set_xlabel("energy levels (molecules)", size=12)
# ax.set_ylabel("cell size (amino acids)", size=12)
ax.set_ylim(limy)
ax.set_xlim(limx2)
ax.set_title("energy levels (molecules)", size=14)
# ax.text(-0.15, 1.05, 'D',fontsize=20, fontweight=1000, transform=axd["D"].transAxes)
# ax.set_position([0.06324052130991283, 0.0329863005050505, 0.2747, 0.2015])  #axd["F"].get_position()
ax.patch.set_alpha(0)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.legend(["model: logistic function"],loc= 'upper left', fontsize=14)
ax.set_ylabel("cell size (amino acids)", fontsize=14)
ax.yaxis.set_label_position("right")



"""
#%% cell size model
#### biosynthesis model

def ode(t,y, params):
    
    # s0 -> a         (kin = (ns*Vmax))
    # a -> -          (kout)

    
    # species
    a = y[0]
    #print(a)
    
    
    # parameters
    kin = 10
    kout = 10
    s0=10
    ns = params[0]
    
    # model
    dydt= np.zeros(1)
    
    dydt[0] = +kin*ns*s0 -kout*a
    
    
    return dydt




a=0
init = [a]


t_span = np.array([0,1e3])


a_ss =[]
e_rate_prod=[]
e_rate_cons=[]
ns_array = [1,2,3,4,5]
ns_array = np.arange(0.1,10,0.1)
for ns in ns_array:

    params = [ns]

    sim = integrate.solve_ivp(ode, t_span, init, args= (params,), atol=1e-5,rtol=1e-5)
    
    t= sim.t
    a, = sim.y
    
    #plt.semilogx(t,a)
    
    a_ss+= [a[-1]]
    
    kin=10
    s0=10
    kout=10
    e_rate_prod += [kin*s0*ns]
    e_rate_cons += [a[-1]*kout]
    
    
    
### cell size data

data_size = [0.76E9, 1.19E9, 1.64E9, 2.15E9, 2.4E9]
data_gr = [0.6,1,1.5,2,2.5]


####cell size model - energy levels
size_extra=3.4E9
size_min=1E8
k = 0.06
m = 32

x = np.array(a_ss)

y = size_min + (size_extra / (1 + np.exp(-k * (x -m))))


#### - Panel G

limy =[1E8,4E9]
limx =[0.2,5]

axd["F"].scatter(data_gr,data_size, color = "r")
axd["F"].set_xlabel("growth rate (h-1)", size=14, color="r")
axd["F"].set_ylabel("cell size (amino acids)", size=14)
axd["F"].set_ylim(limy)
axd["F"].set_xlim(limx)
# axd["C"].set_title("experimental data", loc="center", size=14)
axd["F"].tick_params(colors='red')
axd["F"].legend(["experimental data"],loc= 'lower right', fontsize=14)
panel="F"
axd[panel].text(-0.1, 1.05, "G" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)


limx2 =[0,100]
ax = fig.add_subplot(4,3,10)
ax.plot(x,np.array(y), color = "k")
#ax.set_xlabel("energy levels (molecules)", size=12)
# ax.set_ylabel("cell size (amino acids)", size=12)
ax.set_ylim(limy)
ax.set_xlim(limx2)
ax.set_title("energy levels (molecules)", size=14)
# ax.text(-0.15, 1.05, 'D',fontsize=20, fontweight=1000, transform=axd["D"].transAxes)
ax.set_position([0.06324052130991283, 0.0329863005050505, 0.2747, 0.2015])  #axd["F"].get_position()
ax.patch.set_alpha(0)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.legend(["model: logistic function"],loc= 'upper left', fontsize=14)



#%% growth law monod

# importing sys
import sys
 
# adding Folder_2 to the system path
module_path = '/Users/nicolaskylilis/OneDrive - University of Cyprus/Employment 3 Fellowship UCY/project_cell_growth_model/experiments/model_v04/parameters_calibration_v04'
sys.path.insert(0,module_path)

 
import simulation
import pandas as pd

# cell parameters
fpath_cell_params = "/Users/nicolaskylilis/OneDrive - University of Cyprus/Employment 3 Fellowship UCY/project_cell_growth_model/experiments/model_v06/parameters_calibration_v04/logs/logs_2022.07.27_16.10.35/cell_paramameters.csv"
# fitted parameters
fname = "/Users/nicolaskylilis/OneDrive - University of Cyprus/Employment 3 Fellowship UCY/project_cell_growth_model/experiments/model_v06/parameters_calibration_v04/logs/logs_2022.07.27_16.10.35/log_file 2022.07.27_16.10.35.txt"
df = pd.read_csv(fname, sep="\t", index_col=("opt_run"))
df_best_cost = df[df["cost"] == df["cost"].min()]
particle_params = df_best_cost[df_best_cost.columns[1:]].squeeze()

x_g = [3E-1,1E0,3E0,1E1,3E1,1E2,3E2,1E3,3E3,1E4]
# x = np.geomspace(1E-1, 1E1,100)

y_g =[]
y_m = []
for s0 in x_g:
    s0=int(s0)
    sim, df = simulation.simulate(particle_params, ns=5 , s0 =s0, fpath_params=fpath_cell_params)
    y_g += [df['growth rate (hour-1)'].item()]
    y_m += [df['cell mass (aa)'].item()] 
    

    
axd["H"].semilogx(x_g,y_g, color="k", marker="*", linestyle="--")
axd["H"].set_ylabel("growth rate (h-1)", size=14)
axd["H"].set_xlabel("extracellular  nutrient (molecules)", size=14)
axd["H"].set_title("Monod's growth relation", size=14, loc="center")
panel="H"
axd[panel].text(-0.1, 1.05, "I" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)

#%% Schaechter growth Law

axd["G"].plot(y_g,y_m, color="k", marker="*", linestyle="--")
axd["G"].set_xlabel("growth rate (h-1)", size=14)
axd["G"].set_ylabel("cell size (aa)", size=14)
axd["G"].set_title("Schaechter growth relation", size=14, loc="center")
panel="G"
axd[panel].text(-0.1, 1.05, "H" ,fontsize=20, fontweight=1000, transform=axd[panel].transAxes)
"""