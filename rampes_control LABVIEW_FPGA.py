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

start, start_detec, end = 55, 62, 150
end_detec = 160

plt.close('all')
trajet = r'C:\Users\PM263553\Desktop\These\report\Manuscript\img\3_img\plot_python'
path(trajet)
tra = r"C:\Users\PM263553\Desktop\These\big_projects\in_vitro\data_vivo\control a plot"
name = ["no retour.npy"] + ["avec control{}.npy".format(i+1) for i in range(8)] + ["mice_conv.npy"]
legend = ["sans retour"] + ["control reduce "+str(i) for i in range(8)] + ["mice convergence"]
for i in range(1): #, tra1[tra0, tra1, tra2, tra3], tra44, tra444
    traj = trajet+ "\\" + legend[i]+ "\\"
    path(traj)
    # fit = np.array([7.92060316*2, 2.42161125])    
    # start,end = 0,357000
    # pression_max = 1000
    # pression_min = 1
    # bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
    # bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
    # bitpress = int(np.round(((400-fit[1])/fit[0])))
    # press_max = 30*fit[0]+fit[1]
    # nbit=[1,30]

    print("\nLoading pulses......")
    data_raw = np.load(tra+"\\"+name[i]) #Bubbles #Control
    print(np.shape(data_raw))
    data = data_raw #(np.transpose(data_raw[:,:]))
    # on récupère l'array avec les dimensions suivantes : (314, 32, 605)
    # je veux changer l'ordre des axes pour avoir (605, 32, 314)
    data = np.transpose(data, (2, 0, 1))
    print("Loading done!")    
    n_plot = 50
    plt.subplots(figsize=(10,7))
    for j in range(0,len(data), len(data)//n_plot):
        plt.clf()
        #Values
        DataArrays = data[j]
        ultra_32=(DataArrays[:,4])/16777216.
        ultra_52=(DataArrays[:,5])/16777216.
        BB=(DataArrays[:,16])/16777216.
        C1F0=DataArrays[:,18]

        #limites
        ultra_52_lim=(DataArrays[:,11])/16777216.
        ultra_32_lim=(DataArrays[:,10])/16777216.
        BB_lim=(DataArrays[:,13])/16777216.

        #Dépassement
        CUT_52=-(DataArrays[:,14])/1.
        CUT_32=-(DataArrays[:,15])/1.
        CUT_BB=-(DataArrays[:,8])/1.
        intra_c=-(DataArrays[:,30])/1.

        val = -1
        val_l = []

        for z in range(start_detec,end_detec):
            if (ultra_32[z] > ultra_32_lim[z] or  ultra_52[z] > ultra_52_lim[z] or BB[z] > BB_lim[z]) and val == -1:
                val = z
                #val_l.append(z)
            if val != -1 and z-val > 3:
                C1F0[z] = C1F0[val-1]
        #val_l = [98, 118, 119, 120]
        if val == -1 or float(np.max(C1F0))<=0:
            continue
        #val = val_l[0]       
        plt.plot(C1F0/float(np.max(C1F0)), c = 'blue', label = "Pression émise", linewidth = 3, marker = 'o', markersize = 10)
        plt.scatter(val, C1F0[val]/float(np.max(C1F0)), c = 'red', s=450, label = "Pression non stable", marker='x', linewidths = 5)
        plt.scatter(val-1, C1F0[val-1]/float(np.max(C1F0)), c = 'lime', s=450 , label = "Dernière pression stable", marker='x', linewidths = 5)
        plt.plot([val,end],[C1F0[110]/float(np.max(C1F0)), C1F0[110]/float(np.max(C1F0))], c = 'orange', label = "Pression de stabilisation", linestyle = '--', linewidth = 1.8)          
        # for elt in val_l[1:] :
        #     plt.plot([elt, elt], [0,1], c = 'black', linestyle = '--', linewidth = 1.8)
        
        plt.plot([val, val], [0,1], c = 'black', label = "Détection de cavitation", linestyle = '--', linewidth = 1.8)
        plt.legend(fontsize=15, loc = 'lower left')
        plt.ylim(0,1.1)
        plt.grid(True, which='major')
        plt.ylabel("Valeurs normalisées", fontsize=20)
        plt.xlabel("Fenêtres", fontsize=20)
        plt.xticks(fontsize= 18)
        plt.yticks(fontsize= 18)
        #plt.title('Tir numéro {}'.format(j) ,fontsize=22, fontweight = 'bold')
        plt.xlim(start,end)
        plt.savefig(traj+"pulse_{}".format(j),bbox_inches='tight')
        plt.xlim(val-4,val+7)
        #plt.xlim(val-4,128)
        plt.legend('',frameon=False)
        plt.ylim(C1F0[val-4]/float(np.max(C1F0))*0.97 ,1.03)
        plt.savefig(traj+"Zoom_pulse_{}".format(j),bbox_inches='tight')
    plt.close("all")
sys.exit()
#%%VItro control
trajet = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\transmission crane\\transmission_crane\\pk_control_fuckedup\\"
path(trajet)
tra4 = "D:\\mesure_cavitation_crane_exp_2\\20230426__Control_crane_POS1_30\\"
tra44 = "D:\\mesure_cavitation_crane_exp_2\\20230426__Control_crane_POS1_30_bis\\"
tra444 = "D:\\mesure_cavitation_crane_exp_2\\20230426__Control_crane_POS2_30\\"
legend = ["control crane pos1 30","control crane pos1 30 bis","control crane pos2 30"]
for i, tra in enumerate([tra4, tra44, tra444]): #, tra1[tra0, tra1, tra2, tra3]
    if i == 0:
        continue

    traj = trajet + legend[i]+ "\\"
    path(traj)
    # traj = tra + "Analyse_trueharm\\"
    # path(traj)
    traji = traj + "spectres\\"
    path(traji)
    
    start,end = 0,357000
    
    fit = np.array([4.66745471*2, 5.80567673])
    fit = np.array([7.92060316*2, 2.42161125])
    pression_max = 1000
    pression_min = 1
    bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
    bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
    bitpress = int(np.round(((400-fit[1])/fit[0])))
    press_max = 30*fit[0]+fit[1]
    nbit=[1,30]
    
    print("\nLoading pulses......")
    data_raw = np.load(tra+'\\Bubbles.npy') #Bubbles #Control
    data = (np.transpose(data_raw[:,:]))
    #sys.exit()
    print("Loading done!")
    test_exp=experiment(25000000,400000,start=start,end=end, size_decal = 1271)#1271
    print("\nAdding pulses exp")
    test_exp.add_pulses(data, spacer =30e3)
    print("Adding done!")
    x_press = [indi for indi in range(test_exp.pulses[0].n_window)]

    print("\nPLot des différents indices")
    print(legend[i])
    test_exp.plot_indice_RAMP("treatment",legend[i],traj,x_press)
    for std in range(5,10,2):
            print("STD = ", std)
            test_exp.plot_indice_RAMP_std(legend[i],traj+"MOY_std_{}\\".format(std),still_wind = 40, std_tresh = std, true_harm = False, plot_true = True)
       
    print("plot cartes de pression raw")
    n_pulse = np.shape(data)[0]
    fit = [1,0]
    nbit= [1,n_pulse]
    test_exp.plot_windowed(traj,nbit,fit,20,99,1,legend, ramp = True)
    test_exp.plot_windowed(traj,nbit,fit,40,100,1,legend, ramp = True)
    
    # for j in range(0,n_pulse,25):
    #     print(j)
    #     test_exp.pulses[j].plot(traji+"spec_{}".format(j))
        
#%%VIVO posttreatment
tra4 = "F:\\20230601__603\\"
tra44 = "F:\\20230601__604\\"
tra444 = "D:\\data_vivo\\20230322__mouse_588\\"
legend = ["xp 4 mouse 603","xp 4 mouse 604","xp 3 mouse 588"]
#VIVO_Mouse587_221122_45
for i, tra in enumerate([tra4,tra44]): #, tra1[tra0, tra1, tra2, tra3]

    traj = trajet + legend[i]+ "\\"
    path(traj)
    # traj = tra + "Analyse_trueharm\\"
    # path(traj)
    traji = traj + "spectres\\"
    path(traji)
    
    start,end = 0,357000
    
    fit = np.array([4.66745471*2, 5.80567673])
    fit = np.array([7.92060316*2, 2.42161125])
    pression_max = 1000
    pression_min = 1
    bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
    bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
    bitpress = int(np.round(((400-fit[1])/fit[0])))
    press_max = 100*fit[0]+fit[1]
    nbit=[1,bitmax]
    
    print("\nLoading pulses......")
    data_raw = np.load(tra+'\\Bubbles.npy') #Bubbles #Control
    data = (np.transpose(data_raw[:,:1104]))
    #sys.exit()
    print("Loading done!")
    test_exp=experiment(25000000,1500000,start=start,end=end, size_decal = 1271)#1271
    print("\nAdding pulses exp")
    test_exp.add_pulses(data, spacer =100e3)
    print("Adding done!")
    x_press = [indi for indi in range(test_exp.pulses[0].n_window)]

    print("\nPLot des différents indices")
    print(legend[i])
    test_exp.plot_indice_RAMP("treatment",legend[i],traj,x_press)
    # for std in range(5,10,2):
    #         print("STD = ", std)
    #         test_exp.plot_indice_RAMP_std(legend[i],traj+"MOY_std_{}\\".format(std),still_wind = 40, std_tresh = std, true_harm = False, plot_true = True)
       
    print("plot cartes de pression raw")
    n_pulse = np.shape(data)[0]
    fit = [1,0]
    nbit= [1,n_pulse]
    test_exp.plot_windowed(traj,nbit,fit,20,99,1,legend, ramp = True)
    test_exp.plot_windowed(traj,nbit,fit,40,100,1,legend, ramp = True)
    
    # for j in range(0,n_pulse,25):
    #     print(j)
    #     test_exp.pulses[j].plot(traji+"spec_{}".format(j))
        
#%%
sys.exit()
rep = 20
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

