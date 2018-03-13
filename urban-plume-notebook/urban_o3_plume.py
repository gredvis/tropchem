# -*- coding: utf-8 -*-
"""
Created on Mon 12 Mar 23 11:21:08 2018

Purpose: 
    Simulate the formation of ozone in an urban plume given using
    a simple chemical box model with only a few reactions involving
    NOx and VOC and the QSSA numerical solver.
    
    The code is based on the chemical box model originally developed
    at IAC ETH by Joachim Orb and ported to php by Jörg Mäder, which
    has been used in various lectures on atmospheric chemistry.
    
    Written for tropospheric chemistry lecture, ETH Zurich 2018
    
    The following reactions are considered
    
    (1) ROG + OH -> HO2, k1
    (2) HO2 + NO -> OH + NO2, k2
    (3) NO + O3 -> NO2 + O2, k3
    (4) NO2 + O2 -> NO + O3, k4
    (5) NO2 + OH -> HNO3, k5
    (6) O3 + H2O -> O2 + 2 OH, k6
    (7) 2 HO2 -> O2 + H2O2, k7
    (8) O3 + HO2 -> OH + 2 O2, k8
    (9) OH + HO2 -> H2O + O2, k9
    (10) OH (+CH4) -> HO2 (+CH4), k10
    
    Examples:
    urban_o3_plume()
    
@author: Dominik Brunner, Empa
"""
import numpy as np
from math import exp
import matplotlib.pyplot as plt

def urban_plume(SNOX=0.01,SROG=0.08,CNOX=1.0,CROG=4.0,CO3=30.0,
                T=303.0,doPlot=True,savePlot=False,doCSV=False,
                showHOx=False):
    
    plume = simple_boxmodel(SNOX=SNOX,SROG=SROG,CNOX=CNOX,CROG=CROG,
                            CO3=CO3,T=T)
    
    colors = {}
    colors['O3'] = 'blue'
    colors['ROG'] = 'purple'
    colors['NO'] = 'blue'
    colors['NO2'] = 'darkblue'
    colors['OH'] = 'grey'
    colors['HO2'] = 'darkgrey'
    colors['HNO3'] = 'green'
    colors['H2O2'] = 'orange'
    
    if doPlot:        
        fs = 14
        xlim = (2.,175.)
        fig, (ax1, ax2) = plt.subplots(2)
        ax1.set_xlabel('Distance (km)', fontsize=fs)
        ax1.set_xlim(xlim)
        #ax1.set_ylim(0,160.)
        ax1.set_ylabel('conc ($\mu$g m$^{-3}$)', fontsize=fs)
        ax1.plot(plume['dist'], plume['conc']['O3'],color=colors['O3'],
                 linewidth=4)
        ax1.plot(plume['dist'], plume['conc']['ROG'],color=colors['ROG'])
        ax1.tick_params(axis='y')
        
        #ax2 = ax1.twinx()  # instantiate a second axis using the same x-axis
        color = 'red'
        ax2.set_xlabel('Distance (km)', fontsize=fs)
        ax2.set_xlim(xlim)
        #ax2.set_ylim(0,30.)
        ax2.set_ylabel('conc ($\mu$g m$^{-3}$)', fontsize=fs)
        ax2.plot(plume['dist'], plume['conc']['NO'], color=colors['NO'])
        ax2.plot(plume['dist'], plume['conc']['NO2'], color=colors['NO2'])
        ax2.plot(plume['dist'], plume['conc']['H2O2'], 
                 color=colors['H2O2'], linewidth=4)
        ax2.plot(plume['dist'], plume['conc']['HNO3'], 
                 color=colors['HNO3'], linewidth=4)
        ax2.tick_params(axis='y')
        
        if showHOx:
            ax3 = ax1.twinx()  # instantiate a third axis using the same x-axis
            color = 'green'
            ax3.set_ylim(0,0.2)
            ax3.set_ylabel('conc ($\mu$g m$^{-3}$)', color=color, fontsize=fs)
            ax3.plot(plume['dist'], plume['conc']['OH']*100,color=colors['OH'])
            ax3.plot(plume['dist'], plume['conc']['HO2'], color=colors['HO2'])
            ax3.tick_params(axis='y', labelcolor=color)

        plt.rcParams['figure.figsize']=(12,10)
        #fig.tight_layout()
        plt.show()
        
        if savePlot:
            plt.savefig("plume.png")
            
        # print minimum and maximum concentrations
        species = ['O3','NO','NO2','ROG','HNO3','H2O2','HO2','OH']
        for i in species:
            print("Min %4s:  %8.3f,  Max %4s:  %8.3f" %
                  (i,min(plume['conc'][i]),i,max(plume['conc'][i])))
        

def simple_boxmodel(SNOX=0.01,SROG=0.08,CNOX=1.0,CROG=4.0,CO3=30.0,T=303.0):
    
    # define some constants
    
    # wind speed (m/s) and max. transport distance (km)
    speed = 5.0
    distance = 180.0
    
    # reaction rates
    k1 = 15e-12
    k2 = 3.7e-12 * exp(240./T)
    k10 = 0.0
    k3 = 1.8e-12 * exp(-1370./T)
    k4 = 0.007
    k5 = 6.7e-11 * (T/300.)**-0.6
    k6 = 1.6e-5
    k7 = 2.3e-13 * exp(600./T)
    k8 = 1.1e-14 * exp(-500./T)
    k9= 4.6e-11 * exp(600./T)
    
    # total time (tt), time step (dt) and number of time steps (nt)
    tt = distance * 1000. / speed
    dt = 5.0
    nt = int(tt/dt)
    
    # concentrations as a dictionary of np arrays
    species = ['OH','HO2','NO','NO2','O3','ROG','HNO3','H2O2']
    c = {x: np.zeros(nt) for x in species}
    # destruction and production terms
    d = {x: 0.0 for x in species}
    p = {x: 0.0 for x in species}
        
    # factor eta for converting number density to concentration
    # assuming a pressure of 960 hPa
    eta = 1e-8 * 8.314 * T/(6.023*96000)
    
    # inital concentrations
    c['OH'][0] =  0.0001 / eta
    c['HO2'][0] = 0.1 / eta
    c['NO'][0] = 0.8 * CNOX / eta
    c['NO2'][0] = 0.2 * CNOX / eta
    c['O3'][0] = CO3 / eta
    c['ROG'][0] = CROG / eta
    c['HNO3'][0] = 0.
    c['H2O2'][0] = 0.

    # loop over time steps
    for j in range(0,nt-1):
        
        # add emissions between times 1000 s and 3000 s
        ct = j*dt
        if ct > 1000 and ct <3000:
            sNO = 0.8 * SNOX / eta
            sNO2 = 0.2 * SNOX / eta
            sROG = SROG / eta
        else:
            sNO = 0.
            sNO2 = 0.
            sROG = 0.
        
        d['OH'] = k1*c['ROG'][j]+k10+k5*c['NO2'][j]+k9*c['HO2'][j]
        p['OH'] = k2*c['HO2'][j]*c['NO'][j]+2*k6*c['O3'][j]+k8*c['HO2'][j]*c['O3'][j]
        
        d['HO2'] = k2*c['NO'][j]+2*k7*c['HO2'][j]+k8*c['O3'][j]+k9*c['OH'][j]
        p['HO2'] = k1*c['ROG'][j]*c['OH'][j]+k10*c['OH'][j]
        
        d['NO'] = k2*c['HO2'][j]+k3*c['O3'][j]
        p['NO'] = k4*c['NO2'][j]+sNO
        
        d['NO2'] = k4+k5*c['OH'][j]
        p['NO2'] = k2*c['HO2'][j]*c['NO'][j]+k3*c['O3'][j]*c['NO'][j]+sNO2
        
        d['O3'] = k3*c['NO'][j]+k6+k8*c['HO2'][j]
        p['O3'] = k4*c['NO2'][j]
        
        d['ROG'] = k1*c['OH'][j]
        p['ROG'] = sROG
        
        d['HNO3'] = 1e-30
        p['HNO3'] = k5*c['NO2'][j]*c['OH'][j]
        
        d['H2O2'] = 1e-30
        p['H2O2'] = k7*c['HO2'][j]*c['HO2'][j]

        # loop over species and integrate with QSSA
        for i in species:
            if d[i]*dt < 0.01:
                c[i][j+1] = c[i][j]+(p[i]-d[i]*c[i][j])*dt
            elif d[i]*dt > 10:
                c[i][j+1] = p[i]/d[i]
            else:
                c[i][j+1] = p[i]/d[i]+(c[i][j]-p[i]/d[i])*exp(-d[i]*dt)
    
    # convert to concentrations
    for i in species:
        c[i] = c[i]*eta
        
    #time = np.arange(0,nt)*dt    
    result = {}
    result['time'] = np.arange(0,nt) * dt  
    result['dist'] = result['time'] * speed / 1000. 
    result['conc'] = c
    return result
