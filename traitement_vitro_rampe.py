# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:01:38 2022

@author: PM263553
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from classes import *
import sys 
import gc

gc.collect(generation=2)

#VITRO
tra = 'C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_13\\'
traj='C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_13\\comparaison_bubbles_FFT_new\\'
traj="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\Analyse_RAMPE\\"
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\"
path(traj)
     
doss=["bubbles_0_75","bubbles_80_75","bubbles_50_75"]  
doss=["RAMP_0_75","RAMP_240_75","bubbles_80_75","bubbles_27_75"]
doss=["RAMP_0_80","RAMP_666_40","bubbles_80_75","bubbles_27_75"]
legend=["Eau pure","Sonovue dilué 240 fois","Sonovue dilué 80 fois","Sonovue dilué 27 fois"]
# =============================================================================
# start,end = 1100,23437 + 1100     #zone rouge
# start,end = 23437 + 1100 ,312800     #zone Vide  
# start,end = 1100,31280     #zone verte    
# =============================================================================

start,end = 0,334000
bit_max_shot = 40


test_m_exp=experiment_mult(25000000,1500000,start=start,end=end)
nexp=2
fit = np.array([4.66745471*2, 5.80567673])
fit = np.array([7.92060316*2, 2.42161125])
pression_max = 500
pression_min = 1
bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
bitpress = int(np.round(((400-fit[1])/fit[0])))
press_max = bit_max_shot*fit[0]+fit[1]
nbit=[1,bitmax]


(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]



for j in range(nexp):
    test_m_exp.creat_expe()

for i in range(0,nexp):
    
    print("\nloading data : "+doss[i])
    dossier=tra+doss[i]
    traj_comp = tra+doss[i]
    path(traj_comp)
    data = np.load(dossier+'\\data.npy')

    if i ==0 :
        n_pulse = np.shape(data)[0]
        x_press = valeurs(data[:,:end], press_max)
        
    test_exp=experiment(25000000,1500000,start=start,end=end)
    
    print("adding pulses exp")
    test_exp.add_pulses(data, spacer =100e3)
    test_exp.plot_indice_RAMP(legend[i],traj,x_press)

    print("adding pulses multi exp")
    test_m_exp.add_pulses(data, i, spacer =100e3)
    print('done')
    del data
    gc.collect(generation=2)

nom = "test_ramp"
dossier = traj+nom+"\\"
test_m_exp.plot_indice_RAMP(nom,dossier,x_press,legend)    
# x_n_fenetre = [i for i in range(len(x_press))]  
# nom = "ramp_fenêtre"
# dossier = traj+nom+"\\"
# test_m_exp.plot_indice_RAMP(nom,dossier,x_n_fenetre,legend)   

#%%
nom_doss = "cartes_de_pression\\"
traj_carte = traj + nom_doss
path(traj_carte)
traj1= traj_carte
path(traj1)
print("plot cartes de pression raw")
fit = [1,0]
nbit= [1,n_pulse]
test_m_exp.plot_windowed(traj1,nbit,fit,10,100,1,legend,True)

sys.exit()
    
# x_n_fenetre = [i for i in range(len(x_press))]  
# nom = "ramp_fenêtre"
# dossier = traj+nom+"\\"
# test_m_exp.plot_indice_RAMP(nom,dossier,x_n_fenetre,legend)    
# nom = "ramp_pression"
# dossier = traj+nom+"\\"
# test_m_exp.plot_indice_RAMP(nom,dossier,x_press,legend)    


#%%
# =============================================================================
# print("plotting pulses")
# rep = 30
# bit = 50
# destination = "C:\\Users\\PM263553\\Desktop\\plot\\new_pulses"
# path(destination)
# destination += "\\"
# for a in range(183):
#     test_m_exp.exp[0].pulses[a].plot(destination+"bit _{}".format(a+2))
# =============================================================================
    
sys.exit()
path(traj)
traj1= traj+"raw\\"
path(traj1)
legend=["no Mbs","Mbs=0.25mL","Mbs=0.50mL","Mbs=0.75mL","Mbs=1.00mL"]
# =============================================================================
# print("plot UH")
# test_m_exp.plot_UH_windowed(traj1,nbit,20,99,1,legend)
# print("plot BB")
# test_m_exp.plot_BB_windowed(traj1,nbit,25,99,1,legend)
# print("plot H")
# test_m_exp.plot_H_windowed(traj1,nbit,15,99,1,legend)
# =============================================================================

sys.exit()


#%%

plt.figure(figsize=(19,9))
for j in range(1,30):
    plt.clf()
    i = j*100
    p = b_to_p(fit,i//30)
    
    plt.plot(np.arange(len(data[i]))/25000.,data[i])
    plt.title("Signal temporel d'un bit tiré à {} KPa".format(p),fontsize=30, fontweight = 'bold')
    plt.xlabel('Temps (ms)',fontsize=20)
    plt.ylabel("Amplitude (mV)",fontsize=20, fontweight = 'bold')
    #plt.ylim([Ymin, Ymax])  
    plt.grid(True)
    plt.tight_layout()
    
    
    traj="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter5\\plot_spectre_pc_paul\\pulses_temp_0\\"
    path(traj)
    plt.savefig(traj+'_{}_log.png'.format(i),bbox_inches='tight')


#%%





#VITRO
traj="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_4\\comparaison_bubbles_FFT\\"

doss=["bubbles_0","bubbles_025_50","bubbles_05_50","bubbles_075_50","bubbles_10_50"]       
start,end = 1100,312800
dossier="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_4\\"+doss[1]
data=np.load(dossier+'\\data.npy')
# = 40000
nexp=5
nbit=[1,51]
 
for i in range(nexp):
    dossier="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_4\\"+doss[i]
    print(doss[i])
    test_exp=experiment(31250000,1500000,start=start,end=end) 
    data=np.reshape(np.load(dossier+'\\data.npy'),(1500,524288))
    test_exp.add_pulses(data, spacer =100e3)
    test_exp.plot_indice_component(traj+doss[i])


#%%

Y = data[-100]
X = np.arange(len(Y))/25000.
plt.plot(X,Y)
plt.xlabel("Temps (ms)", fontsize = 20)
plt.ylabel("Amplitude (mV)", fontsize = 20)
plt.title("Pulse ultrasonore tiré à 1.2MPa", fontsize = 30, fontweight="bold")

