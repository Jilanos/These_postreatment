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
general_path = 'F:\\data_vitro\\EXPE_aout_2023\\20230808__PHD_PAUL_control_'
general_path = 'F:\\data_vivo\\'
name_lst = ["20230601__604", '27', '80', '240']



traj_plot = 'C:\\Users\\PM263553\\Desktop\\images_control_anthony\\'
path(traj_plot)


legend = ["souris 604", 'Dilution 1 27', 'Dilution 1 80', 'Dilution 1 240']

#VIVO_Mouse587_221122_45
for i in range(1):
    print(legend[i])
    tra = general_path + name_lst[i] + '\\'
    traj = traj_plot + legend[i]+ "\\"
    path(traj)
    start,end = 0,-1
    
    fit = np.array([7.92060316*2, 2.42161125])
    pression_max = 1400
    pression_min = 1
    bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
    bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
    bitpress = int(np.round(((400-fit[1])/fit[0])))
    max_bit = int(name_lst[i][-2:])*2
    press_max = max_bit*fit[0]+fit[1]
    nbit=[1,max_bit//2]
    
    print("\nLoading pulses......")
    data = np.transpose(np.load(tra+'\\Bubbles.npy'))
    
    
    print("Loading done!")
    test_exp=experiment(25000000,1500000,start=start,end=end,size_decal=1024)
    print("\nAdding pulses exp")
    test_exp.add_pulses(data, spacer =100e3)
    print("Adding done!")
    x_press = [indi for indi in range(test_exp.pulses[0].n_window)]

    print("\nPLot des différents indices")
    std = 7
    test_exp.plot_indice_RAMP_std(legend[i],traj,still_wind = 40, std_tresh = std, true_harm = False, plot_true = True, BB_mult=0.4)
    #test_exp.plot_windowed(traj,nbit,fit,60,100,1,'pour : ' + legend[i], ramp = True)
    #test_exp.plot_windowed(traj,nbit,fit,40,100,1,'pour : ' + legend[i], ramp = True)
sys.exit()



#%%
#VITRO
general_path = 'F:\data_vitro\\EXPE_aout_2023\\ITER_22\\RAMP_27_75'
traj_plot = 'C:\\Users\\PM263553\\Desktop\\images_control_anthony\\'
path(traj_plot)
i = 0

legend = ['RAMP_27_75']
print(legend[i])
tra = general_path + '\\'
traj = traj_plot + "Pulse_27_75" + "\\"
path(traj)
start,end = 0,-1

fit = np.array([7.92060316*2, 2.42161125])
pression_max = 1000
pression_min = 1
bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
bitpress = int(np.round(((400-fit[1])/fit[0])))
max_bit = int(name_lst[i][-2:])*2
press_max = max_bit*fit[0]+fit[1]
nbit=[1,max_bit//2]

print("\nLoading pulses......")
data = np.load(tra+'data.npy')

print("Loading done!")
test_exp=experiment(25000000,1500000,start=start,end=end)
print("\nAdding pulses exp")
test_exp.add_pulses(data, spacer =100e3)
print("Adding done!")
x_press = [indi for indi in range(test_exp.pulses[0].n_window)]

print("\nPLot des différents indices")
std = 5
test_exp.plot_indice_RAMP_std(legend[i],traj,still_wind = 40, std_tresh = std, true_harm = False, plot_true = True)



#%%
from classes import *
test_exp=experiment(25000000,1500000,start=start,end=end)
test_exp.add_pulses(data, spacer =100e3)

#%%

for i in range(len(animal_name_lst)):
    #test_exp.plot_indice_RAMP("treatment",traj,x_press, vivo= True, all_plot= True)
    for std in range(5,10,2):
            print("STD = ", std)
            test_exp.plot_indice_RAMP_std(legend[i],traj+"MOY_std_{}\\".format(std),still_wind = 40, std_tresh = std, true_harm = False, plot_true = False)
       
    print("plot cartes de pression raw")
    n_pulse = np.shape(data)[0]
    fit = [1,0]
    nbit= [1,n_pulse]
    test_exp.plot_windowed(traj,nbit,fit,20,99,1,legend, ramp = True)
    test_exp.plot_windowed(traj,nbit,fit,40,100,1,legend, ramp = True)
    
    for j in range(0,n_pulse,5):
        print(j)
        test_exp.pulses[j].plot(traji+"spec_{}".format(j))
        
    
rep = 1
order = False
p_deb = 11000


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
    
    if i == 0 :
        n_fen, pos_max = valeurs_rampe(data, rep, p_deb)
        start = p_deb
        end = p_deb + 2* n_fen*2048 +1
        test_m_exp.end = end
        test_m_exp.start = start
    
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:rep*(bitmax-1)], i, spacer =100e3)

    print('done')
    del data
    gc.collect(generation=2)
    

nom_doss = "cartes_de_pression\\"
traj1 = traj + nom_doss
path(traj1)
test_m_exp.plot_windowed(traj1,nbit,fit,10,98,rep,legend)    

nom_doss = "cartes_de_pression_différentielles\\"
traj1 = traj + nom_doss
path(traj1)
test_m_exp.plot_windowed_RAMP(traj1,nbit,fit,10,98,rep,legend)

sys.exit()
#%%

nom = "test_ramp"
dossier = traj+nom+"\\"
legend=["no Mbs","Mbs=0.50mL","Mbs=1.00mL"]
test_m_exp.plot_indice_RAMP(nom,dossier,nbit,legend,fit = False)

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

Y = data[100]
X = np.arange(len(Y))/25000.
plt.plot(X,Y)
plt.xlabel("Temps (ms)", fontsize = 25)
plt.ylabel("Amplitude (mV)", fontsize = 25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("Pulse ultrasonore tiré à 1MPa chez la souris", fontsize = 30, fontweight="bold")

