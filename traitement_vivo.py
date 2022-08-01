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
sys.exit()


#%%
#VITRO
tra = "C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\"
traj="C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\comparaison_bubbles_FFT_UHchange\\"
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\monkey_fevrier\\"
traj= tra + "Post_treatment\\"


path(traj)
start,end = 0,-1

exp_vivo=experiment(25000000,500000,start=start,end=end)



print("\nloading data.... ")
dossier=tra + "Bubbles.npy"
data = np.load(dossier)
amp = np.load(tra + "amplitude_history.npy")
data = np.transpose(data)[:-16]
print("data loaded! \n\nadding pulses....")
traj_plot = traj + "temporel\\"
path(traj_plot)
# x = np.arange(len(data[0]))/25000.
# for j in range(100):
#     plott(x,data[j],traj_plot+"pulse_{}".format(j),color='blue',titre = "Pulse temporel nr {}".format(j),Xlabel="Temps (ms)",Ylabel="Magnitude (V)")
# sys.exit()

exp_vivo.add_pulses(data, spacer =30e3)
print('done')
del data
gc.collect(generation=2)

nom="UH_H_BB"
print("plotting map")
exp_vivo.plot_indice(nom,traj,list(amp))
exp_vivo.plot_UH_windowed(traj,10,100)
exp_vivo.plot_UH_windowed2(traj,10,100)
exp_vivo.plot_BB_windowed(traj,10,100)
exp_vivo.plot_H_windowed(traj,10,100)
#%%

#VITRO
tra = "C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\"
traj="C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\comparaison_bubbles_FFT_UHchange\\"
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\s3\\"
traj= tra + "Post_treatment\\"


path(traj)
start,end = 14000,234000

exp_vivo=experiment(25000000,1500000,start=start,end=end)



print("\nloading data.... ")
dossier=tra + "Data.dat"

data_raw = np.fromfile(dossier,dtype=np.int32)
amp_mat = np.load(tra+'amp_bit.npy')
amp = []
for i in amp_mat:
    amp.append(max(i))
    
    
data = np.reshape(data_raw[2:],(data_raw[0],data_raw[1]))
data = data[:-2]

print("data loaded! \n\nadding pulses....")

exp_vivo.add_pulses(data[:-2], spacer =30e3)
print('done')
# del data
# gc.collect(generation=2)

nom="UH_H_BB"
print("plotting map")
exp_vivo.plot_indice(nom,traj,amp)
exp_vivo.plot_UH_windowed(traj,10,100)
exp_vivo.plot_UH_windowed2(traj,10,100)
exp_vivo.plot_BB_windowed(traj,10,100)
exp_vivo.plot_H_windowed(traj,10,100)

#%%

#VITRO
tra = "C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\"
traj="C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\comparaison_bubbles_FFT_UHchange\\"
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\s4\\"
traj= tra + "Post_treatment\\"


path(traj)
start,end = 14000,234000

exp_vivo=experiment(25000000,1500000,start=start,end=end)



print("\nloading data.... ")
dossier=tra + "Data.dat"

data_raw = np.fromfile(dossier,dtype=np.int32)
amp_mat = np.load(tra+'amp_bit.npy')
amp = []
for i in amp_mat:
    amp.append(max(i))
    
    
data = np.reshape(data_raw[2:],(data_raw[0],data_raw[1]))

data = data[:-2]

print("data loaded! \n\nadding pulses....")

exp_vivo.add_pulses(data[:-2], spacer =30e3)
print('done')
# del data
# gc.collect(generation=2)

nom="UH_H_BB"
print("plotting map")
exp_vivo.plot_indice(nom,traj,amp)
exp_vivo.plot_UH_windowed(traj,10,100)
exp_vivo.plot_UH_windowed2(traj,10,100)
exp_vivo.plot_BB_windowed(traj,10,100)
exp_vivo.plot_H_windowed(traj,10,100)

#%%

#VITRO
tra = "C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\"
traj="C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\EXPE_RESULTS\\comparaison_bubbles_FFT_UHchange\\"
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\s5\\"
traj= tra + "Post_treatment\\"


path(traj)
start,end = 14000,234000

exp_vivo=experiment(25000000,1500000,start=start,end=end)



print("\nloading data.... ")
dossier=tra + "Data.dat"

data_raw = np.fromfile(dossier,dtype=np.int32)
amp_mat = np.load(tra+'amp_bit.npy')
amp = []
for i in amp_mat:
    amp.append(max(i))
    
    
data = np.reshape(data_raw[2:],(data_raw[0],data_raw[1]))

data = data[:-2]

print("data loaded! \n\nadding pulses....")

exp_vivo.add_pulses(data[:-2], spacer =30e3)
print('done')
# del data
# gc.collect(generation=2)

nom="UH_H_BB"
print("plotting map")
exp_vivo.plot_indice(nom,traj,amp)
exp_vivo.plot_UH_windowed(traj,10,100)
exp_vivo.plot_UH_windowed2(traj,10,100)
exp_vivo.plot_BB_windowed(traj,10,100)
exp_vivo.plot_H_windowed(traj,10,100)
#%%
exp_vivo.pulses[100].plot(traj)


chemin=traj
path(traj)
chemin_log=chemin+'\\log\\'
chemin_lin=chemin+'\\linear\\'
path(chemin_log)
path(chemin_lin)
plot=np.zeros((4,exp_vivo.n_pulse))
n_plot = 4
temp = np.zeros((2,n_plot))
posi=0

temp = np.zeros((2,n_plot))
for j in range(exp_vivo.n_pulse):
    plot[0,j] = np.mean(exp_vivo.pulses[j].indice_harm[1:-1]) #np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
    plot[1,j] = np.mean(exp_vivo.pulses[j].indice_Uharm[1:4]) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
    plot[2,j] = np.mean(exp_vivo.pulses[j].indice_BB) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
    plot[3,j] = 30 * plot[1,j] - 1 * plot[2,j] 

plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit de bande","Indices représentant les 3xUH-30xBB","Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
nom_img=["harm","U_harm","BB","ratio_"+nom+"_all","ratio_"+nom+"_start","ratio_"+nom+"_end"]
colors=['black','red','peru','forestgreen','dodgerblue','gold']
pression = np.arange(exp_vivo.n_pulse)
fit = np.array([1, 0])
redblue=['b','g','r','teal','black','black','black']
# Stable cavitation dose

for i in range(0,4):
    fig, host = plt.subplots(figsize=(20,11)) # (width, height) in inches
    host.plot(pression,plot[i,:],  c=redblue[i])
    if amp != False:
        par1 = host.twinx()
        par1.plot(pression,amp,color='gray')
    
    plt.title(plot_legend[i],color=redblue[i],fontsize=30, fontweight = 'bold')
    host.set_xlabel('Pression (KPa)',fontsize=20)
    host.set_ylabel('Magnitude [a.u.]', color=redblue[i],fontsize=20)
    par1.set_ylabel('Pulse amplitude [bit IGT]', color=redblue[i],fontsize=20)
    Ymax = 1.01*np.amax(plot[i,:])
    Ymin = 0.99*np.amin(plot[i,:])
    host.set_ylim([Ymin, Ymax])  
    par1.set_ylim([0.99*np.amin(amp),1.01*np.amax(amp)])    
    host.grid(True)
    fig.tight_layout()
    host.set_yscale('linear')
    plt.savefig(chemin_lin+nom_img[i]+'_linear.png',bbox_inches='tight')  
    plt.close()
    # host.set_yscale('log')
    # plt.savefig(chemin_log+nom_img[i]+'_linear.png',bbox_inches='tight') 






