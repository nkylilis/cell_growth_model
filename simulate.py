#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 11:51:15 2021

@author: nicolaskylilis
"""
#%% simulation

from scipy import integrate
import numpy as np

# experimental condotions
ns = 0.5
cl = 0

## initial conditions
#energy
a_0 = 100
# mRNAs
mr_0 = 100
mc_0 = 100
mq_0 = 100
# initiation complexes
icr_0 = 1000
icc_0 = 1000
icq_0 = 1000
# ribocomplexes
rmr_0 = 1000
rmc_0 = 10000
rmq_0 = 10000
# proteins
r_0 = 5000
em_0 = 5e4
q_0 = 1e6
# zombies
zmr_0= 0
zmc_0= 0
zmq_0= 0

init = np.array([ a_0,mr_0, mc_0, mq_0,icr_0,icc_0,icq_0, rmr_0, rmc_0, rmq_0, r_0, em_0, q_0,zmr_0, zmc_0,zmq_0])

# simulation timespan 
t_span = np.array([0,1e6])

# solve
import model_nk02 as model 
sim = integrate.solve_ivp(model.ode_system, t_span, init , method='BDF', args=([ns,cl]))

#%% plotting

# timecourse

t = sim.t

a = sim.y[0,:]
mr = sim.y[1,:]
mc = sim.y[2,:]
mq = sim.y[3,:]
icr = sim.y[4,:]
icc = sim.y[5,:]
icq = sim.y[6,:]
rmr = sim.y[7,:]
rmc = sim.y[8,:]
rmq = sim.y[9,:]
r= sim.y[10,:]
em= sim.y[11,:]
q= sim.y[12,:]
zmr = sim.y[13,:]
zmc = sim.y[14,:]
zmq = sim.y[15,:]

gmax= 1260.0
Kgamma= 7
M= 1.2e9
nx= 300.0
nr= 7549.0


# growth rate
ttrate= (rmr + rmc + rmq)*(gmax*a/(Kgamma + a))
lam= (ttrate/M)


# ribosomal mass fraction @steady_state
ribo_proteins = r + rmr + rmc + rmq + zmr + zmc +zmq + icr + icc + icq
other_proteins = em + q

fr = (nr*ribo_proteins)/((nr*ribo_proteins) +(nx*other_proteins))
fc = (nx*em)/((nr*ribo_proteins) +(nx*other_proteins))
fq = (nx*q)/((nr*ribo_proteins) +(nx*other_proteins))



import matplotlib.pyplot as plt

fig, ((ax1, ax2), (ax3, ax4), (ax5,ax6), (ax7,ax8)) = plt.subplots(4, 2, figsize = [6.4*1.5, 4.8*3], tight_layout = True)

fig.suptitle('E. coli cell growth model; ns=' + str(ns) + ", cl=" + str(cl), fontsize=12, fontweight = "semibold")

# proteome fractions
labels=["ribosomal", "catabolic",'others']
sizes = [round(fr[-1],2), round(fc[-1],2),round(fq[-1],2)]
ax1.pie(sizes,labels = labels, autopct='%1.1f%%')
ax1.set_title("proteome fractions")

# growth rate
ax2.semilogx(t, lam, color ='k')
ax2.set_title("growth rate @ss =" + str(round(lam[-1],3)) + " min-1")
ax2.set_xlabel("time")
ax2.set_ylabel("min-1")

# rnas
ax3.loglog(t,mr, color= "r")
ax3.loglog(t,mc, color= "m")
ax3.loglog(t,mq, color= "b")
ax3.loglog(t,rmr +rmc + rmq, color = 'k')
ax3.set_title("free mRNAs")
labels = ["mr","mc","mq","ribocomplex total"]
ax3.legend(labels)
ax3.set_xlabel("time")
ax3.set_ylabel("molecules number")
ax3.set_ylim(1,1e7)

ax5.loglog(t,rmr, color = 'r')
ax5.loglog(t,rmc, color = 'm')
ax5.loglog(t,rmq, color = 'b')
ax5.loglog(t,rmr +rmc + rmq, color = 'k')
ax5.set_title("active ribocomplexes ")
labels = ["rmr","rmc","rmq","ribocomplex total"]
ax5.legend(labels)
ax5.set_xlabel("time")
ax5.set_ylabel("molecules number")
ax5.set_ylim(1,1e7)

ax7.loglog(t,icr, color = 'r')
ax7.loglog(t,icc, color = 'm')
ax7.loglog(t,icq, color = 'b')
ax7.loglog(t,rmr +rmc + rmq, color = 'k')
ax7.set_title("initiation complexes")
labels = ["icr","icc","icq","ribocomplex total"]
ax7.legend(labels)
ax7.set_xlabel("time")
ax7.set_ylabel("molecules number")
ax7.set_ylim(1,1e7)

# energy
ax4.loglog(t,a, color = 'y')
ax4.set_title("energy @ss =" + str(round(a[-1],2)) + " molecules")
ax4.set_xlabel("time")
ax4.set_ylabel("molecules number")

# proteins
ax6.loglog(t,rmr +rmc + rmq + r, color = 'r')
ax6.loglog(t,em, color = 'm')
ax6.loglog(t,q, color = 'b')
ax6.set_title("cellular proteins")
labels = ["total ribosomal", "metabolic","other"]
ax6.legend(labels)
ax6.set_xlabel("time")
ax6.set_ylabel("molecules number")
ax6.set_ylim(1,1e7)

# ribosomal
ax8.loglog(t,r, color= "r")
ax8.loglog(t,rmr +rmc + rmq + icr + icc +icq, color = 'k')
ax8.loglog(t,zmr + zmc +zmq, color= "g")
ax8.set_title("ribosomes")
labels = ["free ribosomal","ribocomplex total","zombie total"]
ax8.legend(labels)
ax8.set_xlabel("time")
ax8.set_ylabel("molecules number")
ax8.set_ylim(1,1e7)










