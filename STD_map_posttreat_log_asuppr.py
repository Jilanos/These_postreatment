# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 15:45:10 2023

@author: PM263553
"""

import numpy as np
import time
# import picoscope
#import pyFFTw_pyQTgraph_Specgram as pyfftwspecgram
#import pyFFTw_pyQTgraph_Specgram_plot as pyfftwspecgram_plot
import math
import matplotlib.pyplot as plt
import matplotlib
import os
import BBBop_plots as Mpl
from classes import *

trajet = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Analyse_mouse_ajd\\"
nom = ['water_tank','souris_575','souris_588']
traj_0 = "F:\\20230220_PCDsig_TEST_gene_2\\"
trajet = traj_0

for ind, tra in enumerate([traj_0]):
    X= np.load(tra+"ControlArr.npy")
    
    Clock=np.transpose(X[:,0,:])
    i_max=np.argmax(np.max(Clock,axis=1))
    t_sequence=Clock[i_max,-65]/1e6
    Clock_vec = np.max(Clock[:i_max,:],axis=1)
    pulse_size=len(X[:,0,0])
    seq_size=i_max
        
    Clock_mat=np.transpose(X[:,0,:i_max])
    Output2_mat=np.transpose(X[:,2,:i_max])
    Output_mat=np.transpose(X[:,3,:i_max])
    UH1_mat=np.transpose(X[:,4,:i_max])
    UH2_mat=np.transpose(X[:,5,:i_max])
    
    BN1_mat=np.transpose(X[:,7,:i_max])
    

    BN1_mat_std = np.subtract(np.transpose(BN1_mat) , np.mean(BN1_mat[:,0:40],axis = 1))/np.std(BN1_mat[:,0:40],axis = 1)
    UH1_mat_std = np.subtract(np.transpose(UH1_mat) , np.mean(UH1_mat[:,0:40],axis = 1))/np.std(UH1_mat[:,0:40],axis = 1)
    UH2_mat_std = np.subtract(np.transpose(UH2_mat) , np.mean(UH2_mat[:,0:40],axis = 1))/np.std(UH2_mat[:,0:40],axis = 1)
    BN1_mat_std_c = np.array(BN1_mat_std)[2:-7,:]
    UH1_mat_std_c = np.array(UH1_mat_std)[2:-7,:]
    UH2_mat_std_c = np.array(UH2_mat_std)[2:-7,:]
    
    pas = 0.002
    nombre = 10
    
    
    UH_STD_MAt = np.array([UH1_mat_std_c,UH2_mat_std_c])
    std_UH_max= np.max(UH_STD_MAt)
    val_std_uh = np.arange(np.median(UH_STD_MAt),std_UH_max,pas)
    std_y_uh = []
    for std_x_uh in val_std_uh:
        std_y_uh.append(np.count_nonzero(UH_STD_MAt > std_x_uh))
        
    
    BB_STD_MAt = np.array(BN1_mat_std_c)
    std_BB_max= np.max(BB_STD_MAt)
    val_std_BB = np.arange(np.median(BB_STD_MAt),std_BB_max,pas)
    std_y_BB = []
    for std_x_BB in val_std_BB:
        std_y_BB.append(np.count_nonzero(BB_STD_MAt > std_x_BB))
        
    std_opt_UH = np.max(val_std_uh)
    std_opt_BB = np.max(val_std_BB)
    for j in range(len(val_std_uh)-1):
        if std_y_uh[-(j+1)]>5:
            break
        std_opt_UH = val_std_uh[-(j+1)]
        if j == len(val_std_uh)-2:
            print("error bout de liste atteint pour les UH")
    for j in range(len(val_std_BB)-1):
        if std_y_BB[-(j+1)]>5:
            break
        std_opt_BB = val_std_BB[-(j+1)]
        if j == len(val_std_BB)-2:
            print("error bout de liste atteint pour le BB")
    print("std opti UH = {}".format(std_opt_UH))    
    print("std opti BB = {}".format(std_opt_BB))    



    
    

    Val0=np.array(BN1_mat_std_c).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN1_mat_std_c= Val

    Val0=np.array(UH1_mat_std_c).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    UH1_mat_std_c= Val

    Val0=np.array(UH2_mat_std_c).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    UH2_mat_std_c= Val

         

    UH_mat_std = (UH2_mat_std + UH1_mat_std)/2.
    list_born=[np.nanmin(20 * np.log10(np.abs(BN1_mat_std_c))),np.nanmax(20 * np.log10(np.abs(BN1_mat_std_c))),np.nanmin(20 * np.log10(np.abs([UH1_mat_std_c,UH2_mat_std_c]))),np.nanmax(20 * np.log10(np.abs([UH1_mat_std_c,UH2_mat_std_c])))]
    Thr_BN=20 * np.log10(std_opt_BB)
    Thr_UH=20 * np.log10(std_opt_UH)
    Thr_born=[Thr_BN,Thr_UH]
    print(Thr_born)
    fig_scatter,ax = plt.subplots(1,1,figsize=(15,10))
    BNarr=[BN1_mat_std_c,BN1_mat_std_c]
    iUHarr=[UH1_mat_std_c,UH2_mat_std_c]
    Mpl.Scatter_simple_Paul(ax,BNarr,iUHarr,list_born)
    ax.plot([list_born[0],list_born[1]],[Thr_UH,Thr_UH],color='k')
    ax.plot([Thr_BN,Thr_BN],[list_born[2],list_born[3]],color='k')
    ax.set_xlabel('ICD/ std')
    ax.set_ylabel('UCD/ std')
    #ax.set_xlim([-list_born[1],list_born[1]*1.1])
    #ax.set_ylim([-list_born[3],list_born[3]*1.1])
    
    plt.savefig(trajet+nom[ind]+'events_map_baseline.png')
    
    