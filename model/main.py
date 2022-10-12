#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 11:06:30 2022

@author: nicolaskylilis
"""

# packages
import cell_model
import processes.cell_size as cell_size
import pandas as pd
import matplotlib.pyplot as plt

#%% simulate

ns=5
sim,_ = cell_model.simulate(ns)



t = sim.t

a       = sim.y[0,:]
mr      = sim.y[1,:]
mc      = sim.y[2,:]
mq      = sim.y[3,:]
icr     = sim.y[4,:]
icc     = sim.y[5,:]
icq     = sim.y[6,:]
rmr     = sim.y[7,:]
rmc     = sim.y[8,:]
rmq     = sim.y[9,:]
r       = sim.y[10,:]
em      = sim.y[11,:]
q       = sim.y[12,:]
ribo    = sim.y[13,:]

# cell parameters
fpath_params = "cell_parameters.csv"
parms = pd.read_csv(fpath_params, index_col=("Parameter"))["Value"]
gmax=       parms["gmax"]
Kgamma=     parms["Kgamma"]
lenR =      parms["lenR"]
lenO =      parms["lenO"]
lenC =      parms["lenC"]
lenRibo=    parms["lenRibo"]

## cell size 
cell_mass = []
for i in a:
    state_var = [i]
    cell_mass += [cell_size.model(parms,state_var)]

## biosynthesis
# growth rate @steady_state
ttrate= (rmr + rmc + rmq)*(gmax*a/(Kgamma + a))
lam= (ttrate/cell_mass)
# protein mass in aa
rp_mass = ((ribo[-1] + icr[-1] + icc[-1] +icq[-1] + rmr[-1] + rmc[-1] +rmq[-1])*lenRibo) + (r[-1] * lenR)
em_mass = em[-1] * lenC
q_mass  = q[-1] * lenO


#%% plot

fig, axs = plt.subplots(4, 3, tight_layout="tight", figsize=(8,12))

# proteome sectors
x = [rp_mass, em_mass, q_mass]
labels = ["rib", "met", "hsk"] 
axs[0, 1].pie(x, labels = labels, autopct='%1.1f%%')

axs[0, 1].set_title('Proteome sectors allocation')


# growth rate
gr = lam*60
axs[1, 2].plot(t, gr, color = "k")
axs[1, 2].set_xscale("log")

axs[1, 2].set_title('Growth rate')
axs[1, 2].set_xlabel('time')
axs[1, 2].set_ylabel('growth rate (h-1)')
axs[1, 2].legend(["ss=" + str(round(gr[-1],2))])


# cell mass 
axs[1, 0].plot(t, cell_mass, color = "k")
axs[1, 0].set_xscale("log")
axs[1, 0].set_yscale("log")
axs[1, 0].set_ylim([1e8,1e10])

axs[1, 0].set_title('Cell Size')
axs[1, 0].set_xlabel('time')
axs[1, 0].set_ylabel('cell mass (aa)')
axs[1, 0].legend(["ss=" + "{:.2e}".format(cell_mass[-1])])


# cell mass 
rate = ttrate*60
axs[1, 1].plot(t,rate , color = "k")
axs[1, 1].set_xscale("log")
axs[1, 1].set_yscale("log")

axs[1, 1].set_title('Biosynthesis rate')
axs[1, 1].set_xlabel('time')
axs[1, 1].set_ylabel('mass growth (aa/min)')
axs[1, 1].legend(["ss=" + "{:.2e}".format(rate[-1])])


# energy
x =2
y = 0
axs[x, y].plot(t,a , color = "yellow")
axs[x, y].set_xscale("log")
axs[x, y].set_yscale("log")
axs[x, y].set_xlabel('time')
axs[x, y].set_ylabel('molecules')
axs[x, y].set_title('energy')

# mrnas
x =2
y = 1
axs[x, y].plot(t,mr , color = "red")
axs[x, y].plot(t,mc , color = "orange")
axs[x, y].plot(t,mq , color = "green")
axs[x, y].set_xscale("log")
axs[x, y].set_xlabel('time')
axs[x, y].set_ylabel('molecules')
axs[x, y].set_title('mRNAs')



# ic
x =2
y = 2
axs[x, y].plot(t,icr , color = "red")
axs[x, y].plot(t,icc , color = "orange")
axs[x, y].plot(t,icq , color = "green")
axs[x, y].set_xscale("log")
axs[x, y].set_xlabel('time')
axs[x, y].set_ylabel('molecules')
axs[x, y].set_title('ICs')


# ribosome complexex
x =3
y = 0
axs[x, y].plot(t,ribo , color = "blue")
axs[x, y].plot(t,icr+icc+icq , color = "magenta")
axs[x, y].plot(t,rmr+rmc+rmq , color = "red")
axs[x, y].set_xscale("log")
axs[x, y].set_xlabel('time')
axs[x, y].set_ylabel('molecules')
axs[x, y].set_title('ribosomes')
axs[x, y].legend(["free","ICs","ECs"])


# ec
x =3
y = 1
axs[x, y].plot(t,rmr , color = "red")
axs[x, y].plot(t,rmc , color = "orange")
axs[x, y].plot(t,rmq , color = "green")
axs[x, y].set_xscale("log")
axs[x, y].set_xlabel('time')
axs[x, y].set_ylabel('molecules')
axs[x, y].set_title('ECs')


# p
x =3
y = 2
axs[x, y].plot(t,r , color = "red")
axs[x, y].plot(t,em , color = "orange")
axs[x, y].plot(t,q , color = "green")
axs[x, y].set_xscale("log")
axs[x, y].set_xlabel('time')
axs[x, y].set_ylabel('molecules')
axs[x, y].set_yscale("log")
axs[x, y].set_title('proteins')