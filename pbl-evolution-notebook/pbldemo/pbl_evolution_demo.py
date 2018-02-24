# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 21:21:08 2018

Purpose: 
    Demonstrate of the effects of boundary layer mixing on air pollutants
    
    Computes the diurnal evolution of the concentrations of a non-reactive
    or exponentially decaying tracer with constant emissions at the surface.
    
    Written for tropospheric chemistry lecture, ETH Zurich 2018
    
    Examples:
    pbl_evolution(nday=5)
    pbl_evolution(nday=5,lifetime=10)
    
@author: Dominik Brunner, Empa
"""
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

def pbl_evolution(nday=1,lifetime=-1):
    # define some constants
    dt = 0.1    # time step in hours
    h1 = 0.2    # height of nocturnal inversion layer in km
    h2 = 2.0    # height of daytime mixed PBL in km
    rate = 1.0  # mass emitted per hour per km3 (kg)
    
    # time evolution of PBL in hours local time
    t1 = 6.     # start of PBL growth
    t2 = 12.    # time when PBL has reached max. height
    t3 = 18.    # end of daytime PBL, start of PBL descent
    t4 = 21.    # end of PBL descent
    t5 = 24.    # of the the day
    
    # number of time steps to simulate
    ntime = int(nday*24*1./dt+1)
    time = np.arange(0,ntime)*dt
    Mc=np.zeros(ntime)  # mass of pollutant in surface layer C
    Vc=np.zeros(ntime)  # volume of layer C
    Mb=np.zeros(ntime)  # mass of pollutant in residual layer B
    Vb=np.zeros(ntime)  # volume of layer B
    
    # start values
    Mc[0]=0.
    Vc[0]=0.2
    Mb[0]=0.
    Vb[0]=2-0.2

    # outer loop: days
    for k in range(0,nday):
        kb=int(k*24./dt) # base index for day
        # nighttime 0 - 6 LT
        for i in range(kb+1,kb+int(t1/dt)+1):
            Mc[i]=Mc[i-1]+rate*dt
            Vc[i]=h1
            if lifetime > 0:
                Mc[i]=Mc[i]-dt/lifetime*Mc[i]
            Mb[i]=Mb[i-1]
            Vb[i]=Vb[i-1]
            if lifetime > 0:
                Mb[i]=Mb[i]-dt/lifetime*Mb[i]
        # growth of C on the expense of B: transfer mass from B to C
        for i in range(kb+int(t1/dt)+1,kb+int(t2/dt)+1):
            Mbtoc=Mb[i-1]*(h2-h1)/(t2-t1)/Vb[i-1]
            Mc[i]=Mc[i-1]+rate*dt+Mbtoc*dt
            Vc[i]=Vc[i-1]+(h2-h1)/(t2-t1)*dt
            if lifetime > 0:
                Mc[i]=Mc[i]-dt/lifetime*Mc[i]
            Mb[i]=Mb[i-1]-Mbtoc*dt
            Vb[i]=Vb[i-1]-(h2-h1)/(t2-t1)*dt
            if lifetime > 0:
                Mb[i]=Mb[i]-dt/lifetime*Mb[i]
        # afternoon max. PBL
        for i in range(kb+int(t2/dt)+1,kb+int(t3/dt)+1):
            Mc[i]=Mc[i-1]+rate*dt
            Vc[i]=h2
            if lifetime > 0:
                Mc[i]=Mc[i]-dt/lifetime*Mc[i]
        # evening PBL decay: growth of B on expense of C
        for i in range(kb+int(t3/dt)+1,kb+int(t4/dt)+1):
            Mctob=Mc[i-1]*(h2-h1)/(t4-t3)/Vc[i-1]
            Mc[i]=Mc[i-1]+rate*dt-Mctob*dt
            Vc[i]=Vc[i-1]-(h2-h1)/(t4-t3)*dt
            if lifetime > 0:
                Mc[i]=Mc[i]-dt/lifetime*Mc[i]
            Mb[i]=Mb[i-1]+Mctob*dt
            Vb[i]=Vb[i-1]+(h2-h1)/(t4-t3)*dt
            if lifetime > 0:
                Mb[i]=Mb[i]-dt/lifetime*Mb[i]
        # late evening to midnight, constant PBL
        for i in range(kb+int(t4/dt)+1,kb+int(t5/dt)+1):
            Mc[i]=Mc[i-1]+rate*dt
            Vc[i]=h1
            if lifetime > 0:
                Mc[i]=Mc[i]-dt/lifetime*Mc[i]
            Mb[i]=Mb[i-1]
            Vb[i]=h2-h1
            if lifetime > 0:
                Mb[i]=Mb[i]-dt/lifetime*Mb[i]
        
    conc = Mc/Vc
    try:
        conb = Mb/Vb
    except:
        print("Unexpected error:", sys.exc_info()[0])
    
    plt.rcParams['figure.figsize']=(12,8)
    plt.plot(time, conc,color='blue')
    plt.plot(time, conb,color='orange')
    plt.title('PBL evolution',fontsize=20)
    plt.xlabel('time (h)',fontsize=15)
    plt.ylabel('conc ($\mu$g m$^{-3}$)',fontsize=15)
    plt.text(1,27,'surface',color='blue',fontsize=20)
    plt.text(1,24,'residual',color='orange',fontsize=20)
    plt.grid()
    #fig.savefig("test.png")
    plt.show()