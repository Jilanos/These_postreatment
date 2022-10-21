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
tra = 'C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_6\\'
traj='C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_6\\comparaison_bubbles_FFT_new\\'
tra = 'C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_13\\'
traj='C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_13\\comparaison_bubbles_FFT_new\\'



path(traj)

doss=["bubbles_0_75","bubbles_240_75","bubbles_80_75","bubbles_27_75"]      
# =============================================================================
# start,end = 1100,23437 + 1100     #zone rouge
# start,end = 23437 + 1100 ,312800     #zone Vide  
# start,end = 1100,31280     #zone verte    
# =============================================================================

start,end = 27000,200000
start,end = 10000,236000

test_m_exp=experiment_mult(25000000,1500000,start=start,end=end)
nexp=4
fit = np.array([4.66745471*2, 5.80567673])
fit = np.array([7.92060316*2, 2.42161125])
pression_max = 1000
pression_min = 1
bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
bitpress = int(np.round(((400-fit[1])/fit[0])))
nbit=[1,bitmax]




# press = [200,450,700,1000]
# val = [p_to_b(fit,elt) for elt in press ]
# val_m = []
# for elt in val:
#     for i in range(5):
#         val_m.append(elt*30+i)

# i=0
# print("\nloading data : "+doss[i])
# dossier="C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\EXPE_RESULTS\\"+doss[i]
# dossier="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter5\\small_data\\"+doss[i]#plop
# traj_comp = traj
# path(traj_comp)
# data = np.load(dossier+'\\data_of.npy')
# amp = np.load(dossier+'\\amp_of.npy')
# test_exp=experiment(25000000,1500000,start=start,end=end)
# print("adding pulses exp")
# test_exp.add_pulses(data[:30*(bitmax-1)], spacer =100e3)


# print("starting plot")

# plt.figure(figsize=(19,10))
# for elt in val:
#     print("plotting  " + str(elt))
#     for i in range(5):
#         chemin = traj_comp + "\\pulse_{}".format(30*elt + i)
#         path(chemin)
#         chemin += "\\"
#         test_exp.pulses[30*elt+i].plot_W(chemin)
# plt.close("all")
# sys.exit()



for j in range(nexp):
    test_m_exp.creat_expe()

for i in range(nexp):
    print("\nloading data : "+doss[i])
    dossier=tra+doss[i]
    #dossier="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter5\\small_data\\"+doss[i]#plop
    traj_comp = tra+doss[i]
    path(traj_comp)
    data = np.load(dossier+'\\data_o.npy') #
    #amp = np.load(dossier+'\\amp.npy') #_order

    
    
    test_exp=experiment(25000000,1500000,start=start,end=end)
    print("adding pulses exp")
    test_exp.add_pulses(data[:20*(bitmax-1)], spacer =100e3)
    test_exp.plot_indice_component(traj+doss[i]+"\\",nbit,fit,20)
    del test_exp
    
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:20*(bitmax-1)], i, spacer =100e3)
    print('done')
    del data
    gc.collect(generation=2)
sys.exit()
#%%

# traj="C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\comparaison_bubbles_FFT_goodKPa\\"
# traj="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter5\\comparaison_bubbles_FFT_goodKPa\\"
nom_doss = "pression_carte\\"
traj += nom_doss
path(traj)

traj1= traj+"raw\\"
path(traj1)
legend=["no Mbs","Dilution 240","Dilution 80","Dilution 27"]
print("plot UH")
test_m_exp.plot_UH_windowed(traj1,nbit,fit,10,100,1,legend)
test_m_exp.plot_UH_norm_windowed(traj1,nbit,fit,10,100,1,legend)
print("plot BB")
test_m_exp.plot_BB_windowed(traj1,nbit,fit,10,100,1,legend)
print("plot H")
test_m_exp.plot_H_windowed(traj1,nbit,fit,10,100,1,legend)

traj2= traj+"moy\\"
path(traj2)
test_m_exp.plot_UH_windowed(traj2,nbit,fit,10,100,30,legend)
test_m_exp.plot_UH_norm_windowed(traj2,nbit,fit,10,100,30,legend)
test_m_exp.plot_H_windowed(traj2,nbit,fit,10,100,30,legend)
test_m_exp.plot_BB_windowed(traj2,nbit,fit,20,100,30,legend)
sys.exit()
#%%

nom = "vitro"
dossier = traj+nom+"\\"
legend=["no Mbs","Dilution 240","Dilution 80","Dilution 27"]
test_m_exp.plot_indice_together_grp(nom,dossier,nbit,20,legend,fit = list(fit))


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

tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\rats\\expe_postraitement\\20211125_PCDsig_Rat2\\"
traj="C:\\Users\\PM263553\\Desktop\\bis"
start=0
end=-1

i=40
Fe= 2*15.625e6

data = np.load(tra+'PCD_20211125_Bubbles.npy') #
Y = data[:,200]
#data = np.fromfile(tra+"PCD_10112021_Bubbles_{}.dat".format(i),dtype=float)
temp_samp = temp_sample(25000000,1500000,Y,start=start,end=end)
temp_samp.plot(traj)
Y = data[:,200]/2.5/10.
X = np.arange(len(Y))/31250.
plt.plot(X,Y)
plt.xlabel("Temps (ms)", fontsize = 20)
plt.ylabel("Amplitude (mV)", fontsize = 20)
plt.title("Tir ultrasonore tiré à 750kPa", fontsize = 30, fontweight="bold")
plt.title("Exemple de signal ultrasonore enregistré", fontsize = 30, fontweight="bold")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#%%




Y = data[1500,:40000]
X = np.arange(len(Y))/25000.
plt.plot(X,Y)
plt.xlabel("Temps (ms)", fontsize = 20)
plt.ylabel("Amplitude (mV)", fontsize = 20)
plt.title("Tir ultrasonore tiré à 750kPa", fontsize = 30, fontweight="bold")
plt.title("Tir ultrasonore tiré à 750kPa", fontsize = 30, fontweight="bold")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

