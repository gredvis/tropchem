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
#from tkinter.filedialog import asksaveasfilename
#import tkinter as tk 
import matplotlib.pyplot as plt

def urban_plume(SNOx=0.01,SROG=0.08,CNOx=1.0,CROG=4.0,CO3=30.0,
                T=303.0,doPlot=True,savePlot=False,saveCSV=False,
                showHOx=True):
    
    plume = simple_boxmodel(SNOx=SNOx,SROG=SROG,CNOx=CNOx,CROG=CROG,
                            CO3=CO3,T=T)
    
    colors = {}
    colors['O3'] = 'tab:red'
    colors['ROG'] = 'tab:blue'
    colors['NO'] = 'tab:blue'
    colors['NO2'] = 'tab:cyan'
    colors['OH'] = 'tab:olive'
    colors['HO2'] = 'tab:green'
    colors['HNO3'] = 'tab:pink'
    colors['H2O2'] = 'tab:purple'
    
    if doPlot:        
        fs = 14 # font size for axis titles and labels
        lw = 2  # default line width
        xlim = (2.,175.)
        fig, (ax1, ax2) = plt.subplots(2)
        ax1.set_xlabel('Distance (km)', fontsize=fs)
        ax1.set_xlim(xlim)
        #ax1.set_ylim(0,160.)
        ax1.set_ylabel('O$_3$, ROG (ppb)', fontsize=fs)
        ax1.plot(plume['dist'], plume['conc']['O3'],color=colors['O3'],
                 linewidth=lw, label='O$_3$')
        ax1.plot(plume['dist'], plume['conc']['ROG'],color=colors['ROG'],
                 linewidth=lw, label = 'ROG')
        ax1.tick_params(axis='y')
        ax1.legend(bbox_to_anchor=(0.03, 0.86, 0.3, .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0., fontsize=fs)
        ax1.set_ylim(ymin=0)
        y_min, y_max = ax1.get_ylim()
        ax1.fill([5,15,15,5,5],[0,0,y_max,y_max,0],'lightgrey')
        ax1.text(8,y_max/2,'city',rotation='vertical',fontsize='x-large')
        
        #ax2 = ax1.twinx()  # instantiate a second axis using the same x-axis
        color = 'red'
        ax2.set_xlabel('Distance (km)', fontsize=fs)
        ax2.set_xlim(xlim)
        #ax2.set_ylim(0,30.)
        ax2.set_ylabel('NO, NO$_2$, HNO$_3$, H$_2$O$_2$ (ppb)', 
                       fontsize=fs)
        ax2.plot(plume['dist'], plume['conc']['NO'], color=colors['NO'],
                 linewidth=lw, label='NO')
        ax2.plot(plume['dist'], plume['conc']['NO2'], color=colors['NO2'],
                 linewidth=lw, label='NO$_2$')
        ax2.plot(plume['dist'], plume['conc']['H2O2'], color=colors['H2O2'],
                 linewidth=lw, label='H$_2$O$_2$')
        ax2.plot(plume['dist'], plume['conc']['HNO3'], color=colors['HNO3'],
                 linewidth=lw, label='HNO$_3$')
        ax2.tick_params(axis='y')
        ax2.legend(bbox_to_anchor=(0.03, 0.76, 0.3, .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0., fontsize=fs)
        ax2.set_ylim(ymin=0)
        y_min, y_max = ax2.get_ylim()
        ax2.fill([5,15,15,5,5],[0,0,y_max,y_max,0],'lightgrey')
        ax2.text(8,y_max/2,'city',rotation='vertical',fontsize='x-large')
        
        
        if showHOx:
            ax3 = ax1.twinx()  # instantiate a third axis using the same x-axis
            color = 'green'
            #ax3.set_ylim(0,0.2)
            ax3.set_ylabel('OH*100, HO$_2$ (ppt)', color=color, fontsize=fs)
            ax3.plot(plume['dist'], plume['conc']['OH']*1e3*100,
                     color=colors['OH'],linewidth=lw, label='OH*100')
            ax3.plot(plume['dist'], plume['conc']['HO2']*1e3, 
                     color=colors['HO2'],linewidth=lw, label='HO$_2$')
            ax3.tick_params(axis='y', labelcolor=color)
            ax3.legend(bbox_to_anchor=(0.67, 0.86, 0.3, .102), loc=3,
                       ncol=2, mode="expand", borderaxespad=0., fontsize=fs)
            ax3.set_ylim(ymin=0)

        plt.rcParams['figure.figsize']=(12,10)
        #fig.tight_layout()
        plt.show()
        
        if savePlot:
            #initfile = "SNOX-%6.2f_SROG-%6.2f.png" % (SNOx,SROG)
            #root = tk.Tk()
            #root.withdraw()
            #dlg = asksaveasfilename(title='Select file', 
            #    filetypes = (("csv files","*.csv"),("all files","*.*")),
            #    initialfile = initfile.replace(' ',''))
            #fname = dlg
            #fname = initfile.replace(' ','')
            fname = "plume.png"
            if fname != '':
                fig.savefig(fname)
            
        # print minimum and maximum concentrations
        species = ['O3','NO','NO2','ROG','HNO3','H2O2','HO2','OH']
        sstr = "           "
        minstr = "Min (ppb): "
        maxstr = "Max (ppb): "
        for i in species:
            sstr = sstr + "%8s" % i
            minstr = minstr+"%8.3f" % min(plume['conc'][i])
            maxstr = maxstr+"%8.3f" % max(plume['conc'][i])
        print(sstr)
        print(minstr)
        print(maxstr)
        
        if saveCSV:
            #initfile = "SNOX-%6.2f_SROG-%6.2f.csv" % (SNOx,SROG)
            #root = tk.Tk()
            #root.withdraw()
            #dlg = asksaveasfilename(title='Select file', 
            #    filetypes = (("csv files","*.csv"),("all files","*.*")),
            #    initialfile = initfile.replace(' ',''))
            #fname = dlg
            #fname = initfile.replace(' ','')
            fname = "plume.csv"
            if fname != '':
                # create a numpy array with columns 'dist' and species
                columns = sorted(plume['conc'].keys())
                d = np.transpose([plume['conc'][col] for col in columns])
                # add distance
                columns = ['dist']+columns
                d = np.c_[plume['dist'],d]
                np.savetxt(fname, d[0:7199:40,:], fmt='%10.6f', 
                           delimiter=';', comments='',
                           header=';'.join(s.rjust(8) for s in columns))
            #root.destroy()


def simple_boxmodel(SNOx=0.01,SROG=0.08,CNOx=1.0,CROG=4.0,CO3=20.0,T=303.0):
    
    # define some constants
    
    # wind speed (m/s) and max. transport distance (km)
    speed = 5.0
    distance = 180.0
    
    # reaction rates
    k1 = 15e-12
    k2 = 3.7e-12 * exp(240./T)
    k3 = 1.8e-12 * exp(-1370./T)
    k4 = 0.007
    k5 = 6.7e-11 * (T/300.)**-0.6
    k6 = 1.6e-5
    k7 = 2.3e-13 * exp(600./T)
    k8 = 1.1e-14 * exp(-500./T)
    k9= 4.6e-11 * exp(600./T)
    k10 = 0.0
    
    # total time (tt), time step (dt) and number of time steps (nt)
    tt = distance * 1000. / speed
    dt = 5.0
    nt = int(tt/dt)
    
    # number densities c as a dictionary of np arrays
    species = ['OH','HO2','NO','NO2','O3','ROG','HNO3','H2O2']
    c = {x: np.zeros(nt) for x in species}
    # destruction and production terms
    d = {x: 0.0 for x in species}
    p = {x: 0.0 for x in species}
        
    # factor eta for converting number density to ppb
    # assuming a pressure of 960 hPa
    eta = 1e-8 * 8.314 * T/(6.023*96000)
    
    # inital concentrations
    c['OH'][0] =  0.000 / eta
    c['HO2'][0] = 0.0 / eta
    c['NO'][0] = 0.4 * CNOx / eta
    c['NO2'][0] = 0.6 * CNOx / eta
    c['O3'][0] = CO3 / eta
    c['ROG'][0] = CROG / eta
    c['HNO3'][0] = 0.
    c['H2O2'][0] = 0.

    # loop over time steps
    for j in range(0,nt-1):
        
        # add emissions between times 1000 s and 3000 s
        ct = j*dt
        if ct > 1000 and ct <3000:
            sNO = 0.8 * SNOx / eta
            sNO2 = 0.2 * SNOx / eta
            sROG = SROG / eta
        else:
            sNO = 0.
            sNO2 = 0.
            sROG = 0.
        
        d['OH'] = k1*c['ROG'][j]+k10+k5*c['NO2'][j]+k9*c['HO2'][j]
        p['OH'] = k2*c['HO2'][j]*c['NO'][j]+2*k6*c['O3'][j]
        +k8*c['HO2'][j]*c['O3'][j]
        
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
    
    # convert to mole fractions (ppb)
    for i in species:
        c[i] = c[i]*eta
        
    #time = np.arange(0,nt)*dt    
    result = {}
    result['time'] = np.arange(0,nt) * dt  
    result['dist'] = result['time'] * speed / 1000. 
    result['conc'] = c
    return result
