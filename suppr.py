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
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000

gc.collect(generation=2)

#VITRO
tra = 'D:\\code_UH_long\\GENE_MOD\\iter_20\\'
traj='D:\\code_UH_long\\GENE_MOD\\iter_20\\Analyse_PULSE\\'
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\monkey_mircen\\20221026_PCDsig_SHOT6_NHP_SIMIO\\"
#C:\Users\PM263553\Desktop\These\big_projects\data_vivo\s5\
#C:\Users\PM263553\Desktop\These\big_projects\in_vitro\data_vivo\s5
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\s3\\"
traj= tra + "Post_treatment\\"


path(traj)
start,end = 0,-1
fe = 25000000
exp_vivo=experiment(fe,1500000,start=start,end=end)



print("\nloading data.... ")
dossier=tra + "Data.dat"

data_raw = np.fromfile(dossier,dtype=np.int32)
data = np.reshape(data_raw[2:],(data_raw[0],data_raw[1]))
data = data[:-2]
#%%
Y = data_bbb[400]
X = [i / 25000. for i in range(len(Y))]

plott(X,Y,traj+"pulse",color='blue',titre = "Example of in vivo ultrasound signal",Xlabel="Time (ms)",Ylabel="Magnitude (mV)",color_leg = 'black')
ts = temp_sample(fe, 1500000, Y,start, end)
ts.plot(traj)

#%%
sys.exit()
sys.exit()
amp = np.load(tra + "amplitude_history.npy")
n_pulses = len(amp)
data = np.transpose(data)[:len(amp)]
print("data loaded! \n\nadding pulses....")
traj_plot = traj + "pulse_a_pulse\\"
nom="UH_H_BB"
path(traj_plot)
exp_vivo.add_pulses(data, spacer =30e3)

# for j in range(10):
#     i = int(len(amp)/10*j)
#     X = [t/float(fe)*1000 for t in range(len(data[i]))]
#     plott(X,data[i],traj_plot+"{}_pulse_temp.png".format(i),color='black',titre = "Pulse temporel {}".format(i),Xlabel="Temps (ms)",Ylabel="Magnitude")
#     exp_vivo.pulses[i].plot(traj_plot+"{}_pulse_".format(i))
fit = [1,0]
nbit= [1,n_pulses]
# exp_vivo.plot_windowed(traj,nbit,fit,10,99,1,ramp = True)
# exp_vivo.plot_indice(nom,traj,list(amp))
exp_vivo.plot_indice_bis("04 02 2022 CI",traj,1,'aza','aza',False)
sys.exit()

path(traj)

doss=["PULSE_0_75","PULSE_240_75","PULSE_80_75","PULSE_27_75"]  
doss=["PULSE_0_40","PULSE_666_40","TRI_80_75","TRI_27_75"]     
legend=["Eau pure","Sonovue dilué 240 fois","Sonovue dilué 80 fois","Sonovue dilué 27 fois"]

start,end = 10000,236000
start,end = 0,-1

test_m_exp=experiment_mult(25000000,500000,start=start,end=end)
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
    
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:rep*(bitmax-1)], i, spacer =100e3)
    print('done')
    del data
    gc.collect(generation=2)



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




Y = data[1500,:40000]
X = np.arange(len(Y))/25000.
plt.plot(X,Y)
plt.xlabel("Temps (ms)", fontsize = 20)
plt.ylabel("Amplitude (mV)", fontsize = 20)
plt.title("Tir ultrasonore tiré à 750kPa", fontsize = 30, fontweight="bold")
plt.title("Tir ultrasonore tiré à 750kPa", fontsize = 30, fontweight="bold")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

