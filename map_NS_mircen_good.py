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
tra = 'D:\\code_UH_long\\GENE_MOD\\iter_20\\'
traj='D:\\code_UH_long\\GENE_MOD\\iter_20\\Analyse_PULSE\\'
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\monkey_mircen\\20221026_PCDsig_SHOT6_NHP_SIMIO\\"
tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221109_PCDsig_BBB_NHP\\old\\"
tra0 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\04_02_2022\\"
tra00 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\11_06_2021\\"
tra000 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20220602_PCDsig_MONKEY\\"

tra1 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221130_PCDsig_NHP_BBB_2\\old\\"
tra2 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221130_PCDsig_NHP_BBB_2\\"
tra3 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221130_PCDsig_NHP_BBB\\old\\"
tra4 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221130_PCDsig_NHP_BBB\\"
tra5 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221116_PCDsig_BBB_NHP\\"
tra6 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221109_PCDsig_BBB_NHP\\"
tra7 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221109_PCDsig_BBB_NHP\\old\\"
tra8 = "D:\\data_vivo\\20230210_PCDsig_NHP_BB8846_pente\\"
tra9 = "D:\\data_vivo\\20230310__NHP_BY824\\"

name = ["Mircen_juin_21 ouverture d'après fiche manip", "Mircen_juin_22 Ouverture profonde + trajet", "Mircen_fev + belle ouverture", 
        "NS 30/11/22 bis converge 25bits", "NS 30/11/22 bis forcé 70bits", "NS 30/11/22 converge 16bits", "NS 30/11/22 converge 22bits", 
        "NS 16/11/22 converge 55bits", "NS 09/11/22 converge 65bits", "NS 09/11/22 converge 35bits", "Mircen 10/02/23 Ouverture"]
seuil = [1.07, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 1.8, 1.7, 1.5, 1.5]

def more_events(val,seuil = 0.8):
    return 1 + (val-1)*seuil


all_event_log = []
all_event = []
start_i = [0] + [ 10615 for i in range(10)]
end_i = [312500] + [ 248490 for i in range(10)]
fe_i =  [2*15.625e6] + [ 25000000 for i in range(10)]
factor =  [1] + [ 1 for i in range(10)]


#%%
for ind, tra in enumerate([ tra9]):#tra00, tra000, tra0, tra1, tra2, tra3, tra4, tra5, tra6, tra7, tra8
        
    

    # traj= tra + "Post_NS_mircen\\"
    # path(traj)
    start,end = start_i[ind],end_i[ind]
    fe = fe_i[ind]
    average = 5
    exp_vivo=experiment(fe,500000,start=start,end=end)
    
    
    
    print("\nloading data.... " + name[ind])
    dossier=tra + "Bubbles.npy"
    dossier_base=tra + "Control.npy"
    data_bbb = np.load(dossier)
    data_base = np.load(dossier_base)
    amp = np.load(tra + "amplitude_history.npy")
    n_pulses = len(amp)
    data_bbb = np.transpose(data_bbb)[:len(amp)]
    data_base = np.transpose(data_base)
    
    if ind==0 :
        dataaa = np.zeros((396, 524288))
        dataaa[:395] = data_base
        data_base=dataaa
        del dataaa
        
    
    
    print("\nCreating inertial baselin index")
    inertial_dose = []
    for j in range((len(data_base)-1)//average):
        exp_vivo=experiment(fe,500000,start=start,end=end)
        exp_vivo.add_pulses(data_base[j*average:(j+1)*average], spacer =30e3)
        ic = [puls.indice_BB for puls in exp_vivo.pulses]
        inertial_dose.append(np.mean(ic))
    print("Done")
        
    
    print("\nCreating inertial BBB list")
    inertial_bbb_raw = []
    inertial_bbb_norm = []
    events_pente = []
    events_pente_log = []
    for j in range(len(amp)):
        exp_vivo=experiment(fe,500000,start=start,end=end)
        exp_vivo.add_pulse(data_bbb[j], spacer =30e3)
        inertial_bbb_raw.append(exp_vivo.pulses[0].indice_BB)
        inertial_bbb_norm.append((inertial_bbb_raw[-1]/inertial_dose[int(amp[j])-1])*factor[ind])
        if inertial_bbb_norm[-1]>more_events(seuil[ind]) :
            arr = exp_vivo.pulses[0].indice_BB_sliced
            a,b = np.polyfit([i for i in range(len(arr))], arr,1)
            events_pente.append([inertial_bbb_norm[-1],a])
            arr = np.log10(arr)
            a,b = np.polyfit([i for i in range(len(arr))], arr,1)
            events_pente_log.append([inertial_bbb_norm[-1],a])
    all_event.append(events_pente)
    all_event_log.append(events_pente_log)
    np.save(tra + 'event_BB_pente_log.npy',np.array(events_pente_log)  )
    print(len(events_pente), " events")

#%%  
  
nom = ["image_log", "image_lin"]
plt.figure(figsize = (20,11))
for itera, allev in enumerate([all_event_log]):
    plt.clf()
    color = ['c', 'navy', 'blue', 'forestgreen', 'lime', 'purple', 'magenta', 'dimgrey', 'crimson', 'tomato', 'darkorchid']
    style = ['8', 'v', 'o', 'D', 's', 'p', 'x', '^' , '>', '<', 'o']

    for i, exp in enumerate(allev):
        for j, even in enumerate(exp) :
            if j==0:
                plt.scatter(even[0],even[1], c = color[i], marker=style[i], label = name[i])
            else :
                plt.scatter(even[0],even[1], c = color[i], marker=style[i])
    plt.plot([1,1],[-0.22, 0.09],c = 'black', label = 'Baseline')       
    plt.plot([1.5,1.5],[-0.22, 0.09],c = 'red')            
    plt.plot([1,18],[-0.05, -0.05],c = 'sienna', label = "Limite de pente à implémenter")        
    plt.legend(fontsize=17)
    plt.title("Différents événements des exp à Mircen et NS",fontsize=30, fontweight = 'bold')
    plt.xlabel('Valeurs icd de l event (normalisée)',fontsize=20)
    plt.ylabel('Pente du BB',fontsize=20)
    plt.grid(True)
    plt.tight_layout()
    plt.yscale('linear')
    plt.savefig("C:\\Users\\PM263553\\Desktop\\muscles_monkey\\ouplayop" + nom[itera] +".png",bbox_inches='tight')    
    
plt.close("all")
sys.exit()






print("data loaded! \n\nadding pulses....")

exp_vivo.add_pulses(data, spacer =30e3)
print("pulses added!! ")



#%%
tra9 = "D:\\data_vivo\\20230310__NHP_BY824\\"

name = ["Mircen_session_25"]
seuil = [1.5]


all_event_log = []
all_event = []
start_i = [ 10615 for i in range(10)]
end_i = [ 248490 for i in range(10)]
fe_i =   [ 25000000 for i in range(10)]
factor =  [ 1 for i in range(10)]

for ind, tra in enumerate([tra9]):#tra00, tra000, tra0, tra1, tra2, tra3, tra4, tra5, tra6, tra7, tra8

    start,end = start_i[ind],end_i[ind]
    fe = fe_i[ind]
    average = 5
    exp_vivo=experiment(fe,500000,start=start,end=end)
    
    
    
    print("\nloading data.... " + name[ind])
    dossier=tra + "Bubbles.npy"
    dossier_base=tra + "Control.npy"
    data_bbb = np.load(dossier)
    data_base = np.load(dossier_base)
    amp = np.load(tra + "amplitude_history.npy")
    n_pulses = len(amp)
    data_bbb = np.transpose(data_bbb)[:len(amp)]
    data_base = np.transpose(data_base)
    
    print("\nCreating inertial baselin index")
    inertial_dose = []
    for j in range((len(data_base)-1)//average):
        exp_vivo=experiment(fe,500000,start=start,end=end)
        exp_vivo.add_pulses(data_base[j*average:(j+1)*average], spacer =30e3)
        ic = [puls.indice_BB for puls in exp_vivo.pulses]
        inertial_dose.append(np.mean(ic))
    print("Done")
        
    
    print("\nCreating inertial BBB list")
    inertial_bbb_raw = []
    inertial_bbb_norm = []
    events_pente = []
    events_pente_log = []
    for j in range(len(amp)):
        exp_vivo=experiment(fe,500000,start=start,end=end)
        exp_vivo.add_pulse(data_bbb[j], spacer =30e3)
        inertial_bbb_raw.append(exp_vivo.pulses[0].indice_BB)
        inertial_bbb_norm.append((inertial_bbb_raw[-1]/inertial_dose[int(amp[j])-1])*factor[ind])
        if inertial_bbb_norm[-1]>more_events(seuil[ind]) :
            arr = exp_vivo.pulses[0].indice_BB_sliced
            a,b = np.polyfit([i for i in range(len(arr))], arr,1)
            events_pente.append([inertial_bbb_norm[-1],a])
            arr = np.log10(arr)
            a,b = np.polyfit([i for i in range(len(arr))], arr,1)
            events_pente_log.append([inertial_bbb_norm[-1],a])
    all_event.append(events_pente)
    all_event_log.append(events_pente_log)
    np.save(tra + 'event_BB_pente_log.npy',np.array(events_pente_log)  )
    print(len(events_pente), " events")

#%%  
  
nom = ["image_log", "image_lin"]
plt.figure(figsize = (20,11))
for itera, allev in enumerate([all_event_log]):
    plt.clf()
    color = ['c', 'navy', 'blue', 'forestgreen', 'lime', 'purple', 'magenta', 'dimgrey', 'crimson', 'tomato', 'darkorchid']
    style = ['8', 'v', 'o', 'D', 's', 'p', 'x', '^' , '>', '<', 'o']

    for i, exp in enumerate(allev):
        for j, even in enumerate(exp) :
            if j==0:
                plt.scatter(even[0],even[1], c = color[i], marker=style[i], label = name[i])
            else :
                plt.scatter(even[0],even[1], c = color[i], marker=style[i])
    plt.plot([1,1],[-0.22, 0.09],c = 'black', label = 'Baseline')       
    plt.plot([1.5,1.5],[-0.22, 0.09],c = 'red')            
    plt.plot([1,18],[-0.05, -0.05],c = 'sienna', label = "Limite de pente à implémenter")        
    plt.legend(fontsize=17)
    plt.title("Différents événements des exp à Mircen et NS",fontsize=30, fontweight = 'bold')
    plt.xlabel('Valeurs icd de l event (normalisée)',fontsize=20)
    plt.ylabel('Pente du BB',fontsize=20)
    plt.grid(True)
    plt.tight_layout()
    plt.yscale('linear')
    plt.savefig("C:\\Users\\PM263553\\Desktop\\muscles_monkey\\session25" + nom[itera] +".png",bbox_inches='tight')    
    
plt.close("all")
sys.exit()






print("data loaded! \n\nadding pulses....")

exp_vivo.add_pulses(data, spacer =30e3)
print("pulses added!! ")











#%%









sys.exit()




for tra in [tra9]: #, tra3, tra4
    traj= tra + "Post_treatment_same\\"
    
    
    path(traj)
    start,end = 10615,248490
    fe = 25000000
    exp_vivo=experiment(fe,500000,start=start,end=end)
    
    
    
    print("\nloading data.... ")
    dossier=tra + "Bubbles.npy"
    data = np.load(dossier)
    amp = np.load(tra + "amplitude_history.npy")
    n_pulses = len(amp)
    data = np.transpose(data)[:len(amp)]
    print("data loaded! \n\nadding pulses....")
    traj_plot = traj + "pulse_a_pulse\\"
    nom="UH_H_BB"
    path(traj_plot)
    exp_vivo.add_pulses(data, spacer =30e3)
    print("pulses added!! ")
    
    n_plot = 30
    for j in range(n_plot):
        i = int(len(amp)/n_plot*j)
        X = [t/float(fe)*1000 for t in range(len(data[i]))]
        #plott(X,data[i],traj_plot+"{}_pulse_temp.png".format(i),color='black',titre = "Pulse temporel {}".format(i),Xlabel="Temps (ms)",Ylabel="Magnitude")
        exp_vivo.pulses[i].plot(traj_plot+"{}_pulse_".format(i))
    fit = [1,0]
    nbit= [1,n_pulses]
    exp_vivo.plot_windowed(traj,nbit,fit,10,99,1,ramp = True)
    exp_vivo.plot_indice_bis("   ",traj,1,'aza','aza',False)
    exp_vivo.plot_indice("16 11 2022",traj,amp)
sys.exit()
#%%
tra1 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\monkey_mircen\\04_02_2022\\"
tra2 = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\20221130_PCDsig_NHP_BBB_2\\"
for tra in [tra1, tra2]: #, tra3, tra4
    traj= tra + "Post_treatment_specta_comparison\\"
    
    
    path(traj)
    start,end = 10615,248490
    fe = 25000000
    exp_vivo=experiment(fe,500000,start=start,end=end)
    
    
    
    print("\nloading data.... ")
    dossier=tra + "Bubbles.npy"
    data = np.load(dossier)
    amp = np.load(tra + "amplitude_history.npy")
    n_pulses = len(amp)
    data = np.transpose(data)[:len(amp)]
    print("data loaded! \n\nadding pulses....")
    traj_plot = traj + "pulse_a_pulse\\"
    nom="UH_H_BB"
    path(traj_plot)
    exp_vivo.add_pulses(data, spacer =30e3)
    print("pulses added!! ")
    
    n_plot = 100
    for j in range(n_plot):
        exp_vivo.pulses[j].plot(traj+"{}_pulse_bit_{}".format(j,amp[j]))
#%%

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
Y = data[100000:100070,200]/2.5/10.
X = np.arange(len(Y))/31250. + 100000/31250.
plt.plot(X,Y)
plt.xlabel("Time (ms)", fontsize = 25)
plt.ylabel("Magnitude (mV)", fontsize = 25)
plt.title("Tir ultrasonore tiré à 750kPa", fontsize = 30, fontweight="bold")
plt.title("Example of in vivo ultrasound signal", fontsize = 30, fontweight="bold")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#%%

tra = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_20\\PULSE_80_75\\"
traj="C:\\Users\\PM263553\\Desktop\\bis"
start=0
end=-1

i=40
Fe= 2*15.625e6

data = np.load(tra+'data_o.npy')
print(np.shape(data)) #

Y = data[200,50000:56000] - np.mean(data[200,50000:52048])
#data = np.fromfile(tra+"PCD_10112021_Bubbles_{}.dat".format(i),dtype=float)
temp_samp = temp_sample(25000000,1500000,Y,start=start,end=end)
temp_samp.plot(traj)

#Y = data[100000:100070,200]
X = np.arange(len(Y))/31250. + 50000/31250.
fen = np.hanning(2048)
plt.plot(X,Y*fen,color = 'forestgreen')
plt.xlabel("Time (ms)", fontsize = 25)
plt.ylabel("Magnitude (mV)", fontsize = 25)
plt.title("Tir ultrasonore tiré à 750kPa", fontsize = 30, fontweight="bold")
plt.title("Windowing with Hann window",fontsize = 30, fontweight="bold")
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

