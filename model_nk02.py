# -*- coding: utf-8 -*-

"""
Cell growth coarse grain model linking gene expression to growth rate phenotype

Model features:
    - simple catabolism  module lumping nutrient import and catabolism into a single reaction
    
"""

def ode_system(t,y, ns, cl):
    """
    ## cell and reaction parameters
    #general
    M= 1.2e9
    nx= 300.0
    nr= 7549.0
    dm= 0.1
    # catabolism
    s0= 1.0e4
    ns= 0.5
    Vmax= 1.5e+02
    Km= 1.0e3
    #transcription
    thetar= 7.8820000e+02 
    thetax= 1.752
    wr= 1.208e+03
    wc= 1e+02
    wq= 2.0220000e+04
    Kq= 1.2830720e+06
    nq= 4
    # translation
    gmax= 1260.0
    Kgamma= 7
    kb= 0.1
    ku= 1000
    # inhibition
    k_cm= 7.0163574e-03
    cl= 0
    kc = 180
    ni = 50
    fp = 10
    """
    w_str = "6.35155829e+02 706.29282 1.88724196e+02 8.13937502e+01 2.03250116e-01 8.82052732e+07 4.20496491e+01 0.00218189464 10.5075666 2.86976137e+02"
    w = [float(x) for x in w_str.split(" ")]
    
    M= 1.2e9
    nx= 300.0
    nr= 7549.0
    dm= 0.1
    # catabolism
    s0= 1.0e3
    ns= ns
    Vmax= w[8]
    Km= w[9] #Km= parms[0] #parms[0] #1.0e3
    #transcription
    thetar= w[4]* w[3] # 7.8820000e+02 #
    thetax= w[4] # 1.752 #
    wr= w[0]
    wc = w[1] # 1e+02 #
    wq= w[2] # 2.0220000e+04 #
    Kq=  w[5] # 1.2830720e+06 #
    nq= 4
    # translation
    gmax= 1260.0
    Kgamma=w[6] # 7 #
    kb= 0.1
    ku= 1000
    # inhibition
    k_cm= w[7] # 7.0163574e-03 #
    cl= cl
    kc = 180
    ni = 50
    fp = 10
    

    ## cell species
    # energy
    a = y[0] 
    # mrnas
    mr = y[1]
    mc = y[2]
    mq = y[3]
    # initiation complexes
    icr = y[4]
    icc = y[5]
    icq = y[6]
    # ribocomplexes
    rmr= y[7]
    rmc= y[8]
    rmq= y[9]
    # proteins
    r= y[10]
    em= y[11]
    q= y[12]
    # zombies
    zmr= y[13]
    zmc= y[14]
    zmq= y[15]

    
    
    ## dynamic parameters
    # catabolism
    nucat= em*Vmax*s0/(Km + s0)
    # energy & growth rate
    gamma= gmax*a/(Kgamma + a)
    ttrate= (rmr + rmc + rmq)*(gmax*a/(Kgamma + a))
    lam= (ttrate/M)
    # translation
    kr_eff = kc/(1+((rmr/((mr+icr)*(nr/fp)))**ni))
    kc_eff = kc/(1+((rmc/((mc+icc)*(nx/fp)))**ni))
    kq_eff = kc/(1+((rmq/((mq+icq)*(nx/fp)))**ni))
        

    # ode system
    import numpy as np
    dydt = np.zeros(16)
    
    # energy rxns
    dydt[0]= +ns*nucat -ttrate -lam*a;
    # mRNAs rxns
    dydt[1]= +(wr*a/(thetar + a))                 -kb*r*mr +ku*icr +kr_eff*icr -dm*mr -lam*mr
    dydt[2]= +(wc*a/(thetax + a))                 -kb*r*mc +ku*icc +kc_eff*icc -dm*mc -lam*mc
    dydt[3]= +(wq*a/(thetax + a)/(1 + (q/Kq)**nq)) -kb*r*mq +ku*icq +kq_eff*icq -dm*mq -lam*mq
    # initiation complexes
    dydt[4]= kb*r*mr -ku*icr -kr_eff*icr -lam*icr;
    dydt[5]= kb*r*mc -ku*icc -kc_eff*icc -lam*icc;
    dydt[6]= kb*r*mq -ku*icq -kq_eff*icq -lam*icq;
    # ribocomplexes rxns
    dydt[7]= +kr_eff*icr -gamma/nr*rmr -lam*rmr -k_cm*cl*rmr
    dydt[8]= +kc_eff*icc -gamma/nx*rmc -lam*rmc -k_cm*cl*rmc
    dydt[9]= +kq_eff*icq -gamma/nx*rmq -lam*rmq -k_cm*cl*rmq
    # proteins rxns
    dydt[10] = +ku*icr +ku*icc +ku*icq \
               -kb*r*mr -kb*r*mc -kb*r*mq \
               +gamma/nr*rmr \
               +gamma/nr*rmr +gamma/nx*rmc +gamma/nx*rmq \
               -lam*r
    dydt[11]= +gamma/nx*rmc -lam*em
    dydt[12]= +gamma/nx*rmq -lam*q
    # zombies rxns	
    dydt[13]= +k_cm*cl*rmr -lam*zmr
    dydt[14]= +k_cm*cl*rmc -lam*zmc
    dydt[15]= +k_cm*cl*rmq -lam*zmq
    
    return dydt

