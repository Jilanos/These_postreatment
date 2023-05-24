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


def plot(fft_1, fft_2, chemin):
    
        freq_emiss = 400000.
        freq_ech = 25000000.
        num_BB = 5
        n_harm = 4
        freq_harm=np.array([freq_emiss*i for i in range(1,4+2)])
        freq_Uharm=np.array([freq_emiss*(i-0.5) for i in range(1,4+2)])
        freq_BB=np.array([freq_emiss*(i/2.+1.75) for i in range(num_BB)])
        size_window = 2048
        delta_harm = 30e3
        delta_Uharm = 30e3
        spacer = 30e3
        
        delta_BB=0.25*freq_emiss-max(delta_harm,delta_Uharm)-spacer
        frequencies= np.arange(size_window//2) * freq_ech/size_window
        plt.figure(figsize=(20,11))
        spec_1=fft_1[:size_window//2]#size_window
        spec_2=fft_2[:size_window//2]#size_window
        Ymax = 1.01*np.amax([spec_1,spec_2])
        Ymin = 0.99*np.amin([spec_1,spec_2])
        
        plt.plot(frequencies/1e6, spec_1, 'maroon', label="Without skull", linewidth=2.0)
        plt.plot(frequencies/1e6, spec_2, 'darkviolet', label="With skull", linewidth=2.0)
        plt.plot((freq_emiss/1e6 , freq_emiss/1e6), (0, Ymax), 'r-')
        for n in range(num_BB):
            plt.axvspan((freq_BB[n]-delta_BB)/1e6,(freq_BB[n]+delta_BB)/1e6, color='r',alpha=0.2,lw=0)
        plt.plot(((freq_harm-delta_harm)/1e6 , (freq_harm-delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        plt.plot(((freq_harm+delta_harm)/1e6 , (freq_harm+delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        for n in range(n_harm+1):
            plt.axvspan((freq_harm[n]-delta_harm)/1e6,(freq_harm[n]+delta_harm)/1e6, color='b',alpha=0.2,lw=0)
        plt.plot(((freq_Uharm-delta_Uharm)/1e6 , (freq_Uharm-delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        plt.plot(((freq_Uharm+delta_Uharm)/1e6 , (freq_Uharm+delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        for n in range(n_harm+1):
            plt.axvspan((freq_Uharm[n]-delta_Uharm)/1e6,(freq_Uharm[n]+delta_Uharm)/1e6, color='g',alpha=0.2,lw=0)
        plt.grid(True, which='major')
        plt.ylabel('Magnitude',fontsize=20)
        plt.yscale('log')
        plt.xlabel('Frequency [MHz]',fontsize=20)
        plt.legend()
        plt.title('Averaged spectrum of the pulse with visualization of the different frequency bands',fontsize=20, fontweight = 'bold')#Spectre moyenné du pulse avec visualisation des différentes bandes fréquentielles
        plt.xlim([freq_emiss/1e6/5, (freq_harm[n_harm-1]+delta_harm)/1e6*1.01])#(harm_f[-2]-2*harm_df)/1e6])
        plt.ylim([Ymin, Ymax])
        
        plt.savefig(chemin+'spec_comparison.png',bbox_inches='tight')
        plt.close("all")
sys.exit()
gc.collect(generation=2)

#VITRO
tra = 'D:\\code_UH_long\\GENE_MOD\\iter_20\\'
traj='D:\\code_UH_long\\GENE_MOD\\iter_20\\Analyse_PULSE\\'
tra = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\'
traj='C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\Analyse_PULSE\\'


path(traj)

doss=["PULSE_0_75","PULSE_240_75","PULSE_80_75","PULSE_27_75"]  
doss=["PULSE_0_40","PULSE_666_40","TRI_80_75","TRI_27_75"]     
legend=["Eau pure","Sonovue dilué 240 fois","Sonovue dilué 80 fois","Sonovue dilué 27 fois"]

start,end = 10000,236000
start,end = 0,-1

test_m_exp=experiment_mult(25000000,1500000,start=start,end=end)
nexp=2
fit = np.array([4.66745471*2, 5.80567673])
fit = np.array([7.92060316*2, 2.42161125])
pression_max = 500
pression_min = 1
bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
bitpress = int(np.round(((400-fit[1])/fit[0])))
nbit=[1,bitmax]

rep = 6
order = True

for j in range(nexp):
    test_m_exp.creat_expe()

for i in range(0,nexp):
    print("\nloading data : "+doss[i])
    dossier=tra+doss[i]
    traj_comp = tra+doss[i]
    path(traj_comp)
    if not(order):
        data = np.load(dossier+'\\data.npy') #
        amp = np.load(dossier+'\\amp.npy') #_order
        print("Reordering")
        data, amp = reorder(data,amp)
        print("Saving ordered datas")
        np.save(dossier+'\\data_o.npy',data)
        np.save(dossier+'\\amp_o.npy',amp)
    else :
        data = np.load(dossier+'\\data_o.npy') #
    sys.exit()
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:rep*(bitmax-1)], i, spacer =100e3)
    print('done')
    del data
    gc.collect(generation=2)
#%%

tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_20\PULSE_0_75\\"
doss=["PULSE_0_40","PULSE_666_40","TRI_80_75","TRI_27_75"]    
doss=["PULSE_80_75","PULSE_27_75"]   
legend=["Sonovue dilué 80 fois","Sonovue dilué 27 fois"]

start,end = 10000,236000
start,end = 0,-1

test_m_exp=experiment_mult(25000000,1500000,start=start,end=end)
test_m_exp.creat_expe()
nexp=1
fit = np.array([4.66745471*2, 5.80567673])
fit = np.array([7.92060316*2, 2.42161125])
pression_max = 1000
pression_min = 1
bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
bitpress = int(np.round(((400-fit[1])/fit[0])))
nbit=[1,bitmax]
rep = 20
data = np.load(tra+'data_o.npy')
test_m_exp.add_pulses(data[:rep*(bitmax-1)], 0, spacer =100e3)
#%%
nom_doss = "cartes_de_pression\\"
traj_carte = traj + nom_doss
path(traj_carte)
traj1= traj_carte+"raw\\"
path(traj1)

print("plot cartes de pression raw")
test_m_exp.plot_windowed(traj1,nbit,fit,10,100,1,legend)

traj2= traj_carte+"moy\\"
path(traj2)
print("plot cartes de pression moyennées")
test_m_exp.plot_windowed(traj2,nbit,fit,10,100,rep,legend)

#%%
from classes import *
start,end = 10000,236000

test_m_exp=experiment_mult(25000000,1500000,start=start,end=end)

for j in range(nexp):
    test_m_exp.creat_expe()

for i in range(0,nexp):
    print("\nloading data : "+doss[i])
    dossier=tra+doss[i]
    traj_comp = tra+doss[i]
    path(traj_comp)
    data = np.load(dossier+'\\data_o.npy') #

    test_exp=experiment(25000000,1500000,start=start,end=end)
    print("adding pulses exp")
    test_exp.add_pulses(data[:rep*(bitmax-1)], spacer =100e3)
    test_exp.plot_indice_component(traj+doss[i]+"_Components\\",nbit,fit,rep)
    test_exp.plot_indice_bis(legend[i],traj,rep,nbit,legend,fit)
    del test_exp
    
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:rep*(bitmax-1)], i, spacer =100e3)
    print('done')
    del data
    gc.collect(generation=2)

nom = "différents_indices"
dossier = traj+nom+"\\"
test_m_exp.plot_indice_together_grp(nom,dossier,nbit,rep,legend,fit = list(fit))

sys.exit()
#%%

import numpy as np
import matplotlib.pyplot as plt
import os
from classes import *
import sys 
import gc

gc.collect(generation=2)

#VITRO
tra = 'D:\\code_UH_long\\GENE_MOD\\iter_20\\'
traj='D:\\code_UH_long\\GENE_MOD\\iter_20\\Analyse_PULSE\\'
tra = 'D:\\data_vitro\\POS_2\\'
traj='C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\Analyse_Crane\\POS_2\\'



path(traj)

doss=["PULSE_0_75","PULSE_240_75","PULSE_80_75","PULSE_27_75"]  
doss=["20230307_PCDsig_PULSE_crane_water_1","20230307_PCDsig_PULSE_crane_bubbles_1","20230307_PCDsig_PULSE_nocrane_bubbles_1","TRI_27_75"]     
legend=["Crâne + eau","Crâne + MBs","Pas de crâne + MBs","Sonovue dilué 27 fois"]

start,end = 80415,249290
start,end = 0,-1

test_m_exp=experiment_mult(25000000.,400000.,start=start,end=end,delta_harm =30e3,delta_Uharm =30e3,spacer = 30e3)
nexp=3
fit = np.array([7.92060316*2, 2.42161125])
fit = np.array([1, 0])
pression_max = 500
pression_min = 1
bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
bitpress = int(np.round(((400-fit[1])/fit[0])))

bitmax = 30
nbit=[1,bitmax]

rep = 30

order = True

for j in range(nexp):
    test_m_exp.creat_expe()

for i in range(0,nexp):
    print("\nloading data : "+doss[i])
    dossier=tra+doss[i]
    traj_comp = tra+doss[i]
    path(traj_comp)
    if not(order):
        data = np.load(dossier+'\\ControlArr.npy') #
        amp = np.load(dossier+'\\amp.npy') #_order
        print("Reordering")
        data, amp = reorder(data,amp)
        print("Saving ordered datas")
        np.save(dossier+'\\data_o.npy',data)
        #np.save(dossier+'\\amp_o.npy',amp)
    else :
        data = np.transpose(np.load(dossier+'\\Control.npy')) #
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:rep*(bitmax-1)], i, spacer =30e3)
    print('done')
    del data
    gc.collect(generation=2)

nom_doss = "cartes_de_pression\\"
traj_carte = traj + nom_doss
path(traj_carte)
traj1= traj_carte+"raw\\"
path(traj1)

print("plot cartes de pression raw ")
test_m_exp.plot_windowed(traj1,nbit,fit,10,100,1,legend)

traj2= traj_carte+"moy\\"
path(traj2)
print("plot cartes de pression moyennées")
test_m_exp.plot_windowed(traj2,nbit,fit,10,100,rep,legend)

#%%
import numpy as np
import matplotlib.pyplot as plt
import os
from classes import *
import sys 
import gc

gc.collect(generation=2)

#VITRO
for n_azkrjhakh in range(1,6):
    if n_azkrjhakh==1:
        bitmax, rep = 20, 30
    else:
        bitmax, rep = 30, 40
    tra = 'D:\\data_vitro\\POS_{}\\'.format(n_azkrjhakh)
    traj='C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\Analyse_Crane\\POS_{}\\'.format(n_azkrjhakh)
    
    
    
    path(traj)
    print("\n POS_{}".format(n_azkrjhakh))
    
    doss=["PULSE_0_75","PULSE_240_75","PULSE_80_75","PULSE_27_75"]  
    doss=["20230307_PCDsig_PULSE_crane_water_{}".format(n_azkrjhakh-1),"20230307_PCDsig_PULSE_crane_bubbles_{}".format(n_azkrjhakh-1),"20230307_PCDsig_PULSE_nocrane_bubbles_{}".format(n_azkrjhakh-1),"TRI_27_75"]     
    legend=["Crâne + eau","Crâne + MBs","Pas de crâne + MBs","Sonovue dilué 27 fois"]
    
    start,end = 80415,249290
    
    test_m_exp=experiment_mult(25000000.,400000.,start=start,end=end,delta_harm =30e3,delta_Uharm =30e3)
    nexp=3
    fit = np.array([7.92060316*2, 2.42161125])
    fit = np.array([1, 0])
    pression_max = 500
    pression_min = 1
    bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
    bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
    bitpress = int(np.round(((400-fit[1])/fit[0])))
    
    bitmax = 20
    nbit=[1,bitmax]
    
    rep = 30
    
    order = True
    
    
    for j in range(nexp):
        test_m_exp.creat_expe()
    
        
    for i in range(0,nexp):
        print("\nloading data : "+doss[i])
        dossier=tra+doss[i]
        traj_comp = tra+doss[i]
        path(traj_comp)
        data = np.transpose(np.load(dossier+'\\Control.npy'))
    
        test_exp=experiment(25000000.,400000.,start=start,end=end,delta_harm =30e3,delta_Uharm =30e3)
        print("adding pulses exp")
        test_exp.add_pulses(data[:rep*(bitmax-1)], spacer =30e3)
        test_exp.plot_indice_component(traj+doss[i]+"_Components\\",nbit,fit,rep)
        test_exp.plot_indice_bis(legend[i],traj,rep,nbit,legend,fit)
        del test_exp
        
        print("adding pulses multi exp")
        test_m_exp.add_pulses(data[:rep*(bitmax-1)], i, spacer =30e3)
        print('done')
        del data
        gc.collect(generation=2)
    
    nom = "différents_indices"
    dossier = traj+nom+"\\"
    test_m_exp.plot_indice_together_grp(nom,dossier,nbit,rep,legend,fit = list(fit))
    
    
    nom_doss = "cartes_de_pression\\"
    traj_carte = traj + nom_doss
    path(traj_carte)
    traj1= traj_carte+"raw\\"
    path(traj1)
    
    print("plot cartes de pression raw ")
    test_m_exp.plot_windowed(traj1,nbit,fit,10,100,1,legend)
    
    traj2= traj_carte+"moy\\"
    path(traj2)
    print("plot cartes de pression moyennées")
    test_m_exp.plot_windowed(traj2,nbit,fit,10,100,rep,legend)

sys.exit()


#%%
import numpy as np
import matplotlib.pyplot as plt
import os
from classes import *
import sys 
import gc

gc.collect(generation=2)

for n_azkrjhakh in range(1,6):
    if n_azkrjhakh==1:
        bitmax, rep = 20, 30
    else:
        bitmax, rep = 30, 40
    #VITRO
    tra = 'D:\\data_vitro\\POS_{}\\'.format(n_azkrjhakh)
    traj='C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\Analyse_Crane\\POS_{}\\'.format(n_azkrjhakh)
    path(traj)
    
    doss=["PULSE_0_75","PULSE_240_75","PULSE_80_75","PULSE_27_75"]  
    doss=["20230307_PCDsig_PULSE_crane_bubbles_{}".format(n_azkrjhakh-1),"20230307_PCDsig_PULSE_nocrane_bubbles_{}".format(n_azkrjhakh-1),"TRI_27_75"]     
    legend=["Crâne + MBs","Pas de crâne + MBs","Sonovue dilué 27 fois"]
    
    start,end = 80415,249290
    
    nexp=2
    fit = np.array([1, 0])
    nbit=[1,bitmax]
    

    harm = np.zeros((2))
    Uharm = np.zeros((2))
    bb = np.zeros((2))
    ffttt = []
    
    print("\n POS_{}".format(n_azkrjhakh))
    for i in range(0,nexp):
        dossier=tra+doss[i]
        traj_comp = tra+doss[i]
        path(traj_comp)
        data = np.transpose(np.load(dossier+'\\Control.npy'))
        test_exp=experiment(25000000.,400000.,start=start,end=end,delta_harm =30e3,delta_Uharm =30e3)
        print("adding pulses exp")
        test_exp.add_pulses(data[:rep*(bitmax-1)], spacer =30e3)
        ffttt.append(np.mean([test_exp.pulses[-ik].fftsimple for ik in range(1,21)], axis = 0))
        
        
        harm[i] = 20 * np.log10(np.mean([test_exp.pulses[-ik].indice_harm[1] for ik in range(1,21)]))
        Uharm[i] = 20 * np.log10(np.mean([test_exp.pulses[-ik].indice_Uharm[2] for ik in range(1,21)]))
        bb[i] = 20 * np.log10(np.mean([test_exp.pulses[-ik].indice_BB for ik in range(1,21)]))
        #print('harm value : {}'.format(harm[i]))

        
        print('done')
        del data
        gc.collect(generation=2)
    plot(ffttt[1], ffttt[0], dossier)
    print('DB  absoption crane: {}'.format(harm[1]-harm[0]))
    print('DB  diff  BB: {}'.format(bb[1]-bb[0]))
    print('DB  diff  UH: {}'.format(Uharm[1]-Uharm[0]))



#%%

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

def plot_spec_out(fft_1, chemin):
    
        freq_emiss = 1500000.
        freq_ech = 25000000.
        num_BB = 5
        n_harm = 4
        freq_harm=np.array([freq_emiss*i for i in range(1,4+2)])
        freq_Uharm=np.array([freq_emiss*(i-0.5) for i in range(1,4+2)])
        freq_BB=np.array([freq_emiss*(i/2.+1.75) for i in range(num_BB)])
        size_window = 2048
        delta_harm = 50e3
        delta_Uharm = 50e3
        spacer = 100e3
        
        delta_BB=0.25*freq_emiss-max(delta_harm,delta_Uharm)-spacer
        frequencies= np.arange(size_window//2) * freq_ech/size_window
        plt.figure(figsize=(10,4))
        spec_1=fft_1[:size_window//2]#size_window
        Ymax = 1.01*np.amax([spec_1])
        Ymin = 0.99*np.amin([spec_1])
        
        plt.plot(frequencies/1e6, spec_1, 'k', label="Without skull", linewidth=2.0)
        plt.plot((freq_emiss/1e6 , freq_emiss/1e6), (0, Ymax), 'r-')
        for n in range(num_BB):
            plt.axvspan((freq_BB[n]-delta_BB)/1e6,(freq_BB[n]+delta_BB)/1e6, color='r',alpha=0.2,lw=0)
        plt.plot(((freq_harm-delta_harm)/1e6 , (freq_harm-delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        plt.plot(((freq_harm+delta_harm)/1e6 , (freq_harm+delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        for n in range(n_harm+1):
            plt.axvspan((freq_harm[n]-delta_harm)/1e6,(freq_harm[n]+delta_harm)/1e6, color='b',alpha=0.2,lw=0)
        plt.plot(((freq_Uharm-delta_Uharm)/1e6 , (freq_Uharm-delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        plt.plot(((freq_Uharm+delta_Uharm)/1e6 , (freq_Uharm+delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        for n in range(n_harm+1):
            plt.axvspan((freq_Uharm[n]-delta_Uharm)/1e6,(freq_Uharm[n]+delta_Uharm)/1e6, color='g',alpha=0.2,lw=0)
        plt.grid(True, which='major')
        plt.ylabel('Magnitude',fontsize=20)
        plt.yscale('log')
        plt.xlabel('Frequency [MHz]',fontsize=20)
        plt.title('Spectrum',fontsize=20, fontweight = 'bold')#Spectre moyenné du pulse avec visualisation des différentes bandes fréquentielles
        plt.xlim([freq_emiss/1e6/5, (freq_harm[n_harm-1]+delta_harm)/1e6*1.01])#(harm_f[-2]-2*harm_df)/1e6])
        plt.ylim([Ymin, Ymax])
        plt.xticks([1.5,3,4.5,6],["f\u2080",'2f\u2080',"3f\u2080","4f\u2080"],fontweight='demi', fontsize=15)
        plt.savefig(chemin+'spec_comparison.png',bbox_inches='tight')
        plt.close("all")
spec = test_m_exp.exp[0].pulses[1166].fftsimple
plot_spec_out(spec,"C:\\Users\\PM263553\\Desktop\\These\\presentation\\ISTU\\fig\\")

#%%
#VITRO
tra = 'D:\\code_UH_long\\GENE_MOD\\iter_20\\'
traj='D:\\code_UH_long\\GENE_MOD\\iter_20\\Analyse_PULSE\\'
tra = 'D:\\data_vitro\\iter_21\\PULSE_0_75\\'
traj='C:\\Users\\PM263553\\Desktop\\These\\presentation\\ISTU\\fig\\'


path(traj)

data = np.load(tra+'\\data_o.npy') #
# start,end = 10000,236000
# start,end = 0,-1

# test_m_exp=experiment_mult(25000000,1500000,start=start,end=end)
# nexp=2
# fit = np.array([4.66745471*2, 5.80567673])
# fit = np.array([7.92060316*2, 2.42161125])
# pression_max = 500
# pression_min = 1
# bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
# bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
# bitpress = int(np.round(((400-fit[1])/fit[0])))
# nbit=[1,bitmax]
np.shape(data)
#%%
plt.close('all')
plt.figure(figsize = (24,6))
plt.subplot(1,2,1)
import random
def window(arr):
    arr_out = np.ones(len(arr))
    n = int(0.5*25000.)
    for j in range(n+6250):
        if j<6250:
            arr_out[-j] = 0
            arr_out[j] = 0
        else :
            
            arr_out[j] = float(j-6250)/n
            arr_out[-j] = float(j-6250)/n
    return arr_out*100
  
def window2(arr):
    arr_out = np.ones(len(arr))
    n = int(13*25000.)
    for j in range(n+6250):
        if j<6250:
            arr_out[-j] = 0
            arr_out[j] = 0
        else :
            
            arr_out[j] = float(j-6250)/n
    return arr_out*100
  
x = np.arange(250000 + 12500)/25000.
Y = np.sin(x*1500*2*np.pi)
win = window(Y)
win_bruité = [elt + (random.randint(0,70)-35)/10 for elt in win]
Y_mod = np.array(Y)*np.array(win_bruité)

plt.plot(x,Y_mod, color = 'navy')
plt.title("A) Signal émis trapézoïdal",fontsize=25, fontweight = 'bold')
plt.xlabel("Temps (ms)",fontsize=20, fontweight = 'demi')
plt.ylabel("Amplitude (mV)",fontsize=20, fontweight = 'demi')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True)

    
x = np.arange(25000 * 15 + 12500)/25000.
Y = np.sin(x*1500*2*np.pi)
win = window2(Y)
win_bruité = [elt + (random.randint(0,70)-35)/10 for elt in win]
Y_mod = np.array(Y)*np.array(win_bruité)

plt.subplot(1,2,2)
plt.plot(x,Y_mod, color = 'navy')
plt.title("B) Signal émis rampe",fontsize=25, fontweight = 'bold')
plt.xlabel("Temps (ms)",fontsize=20, fontweight = 'demi')
plt.ylabel("Amplitude (mV)",fontsize=20, fontweight = 'demi')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid(True)
#%%
num = 4 * 4.85 * 0.49 / 3 * np.pi
l = 1483/400000.
D = 70e-3
f = 60e-3
v = num *l**3 *(f/D)**4
print(v)
