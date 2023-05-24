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

def Mat_Pr_Gen(Output_mat,maxPressure,val_min,val_max):
    b=maxPressure/(1-val_min/val_max)
    a=b/(-val_max)
    Output_IGT_mat=a*np.array(Output_mat)+b
    return Output_IGT_mat


def glissante(array, n):
    out = np.copy(array)
    for j in range(n-1,len(array)):
        out[j] = np.mean([array[j-k] for k in range(n)])
    return out



tra1 = 'D:\\data_vivo\\20230221_PCDsig_mouser_574\\'
tra2 = 'D:\\data_vivo\\20230221_PCDsig_mouser_575\\'
tra3 = 'D:\\data_vivo\\20230221_PCDsig_mouser_588\\'


tra1 = 'D:\\data_vivo\\20230215_PCDsig_Paul_Manip_souris_374\\'
tra2 = 'D:\\data_vivo\\20230215_PCDsig_Paul_Manip_souris_375\\'
tra3 = 'D:\\data_vivo\\20230215_PCDsig_Paul_Manip_souris_588\\'
trajet = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Analyse_FPGA\\"
trajet_base = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Tirs multicanaux\\"
path(trajet_base)
tra1 = 'F:\\20230228_PCDsig_WATER_test_0\\'

legend = ["xp 2 mouse 74","xp 2 mouse 575","xp 2 mouse 588"]
legend = ["water_0","mouse 375","mouse 588"]
plt.close("all")
plt.figure(figsize=(15,8))
#VIVO_Mouse587_221122_45
for ind, tra in enumerate([tra1, tra2, tra3]) :
    trajet = trajet_base+legend[ind] + "\\"
    path(trajet)
    i=0
    #data_Con = np.load(tra+'\\ControlArr.npy') #ControlArr
    data_arr = np.load(tra+'\\BubblesArr.npy') #ControlArr #BubblesArr
    rat = 150/16777216.
    amplitude = 170
    fit = np.array([7.92060316E-3, 2.42161125E-3])
    Pr_max =(amplitude*fit[0]+fit[1] )*1000
    
    
    Clock=np.transpose(data_arr[:,0,:])
    i_max=np.argmax(np.max(Clock,axis=1))
    t_sequence=Clock[i_max,-65]/1e6
    Output2_mat=np.transpose(data_arr[:-5,2,:i_max])
    C1F0_mat=np.transpose(data_arr[:-5,18,:i_max])/16777216.  
    C2F0_mat=np.transpose(data_arr[:-5,19,:i_max])/16777216.  
    C3F0_mat=np.transpose(data_arr[:-5,20,:i_max])/16777216.  
    C4F0_mat=np.transpose(data_arr[:-5,21,:i_max])/16777216.  
    
    UH32_mat=np.transpose(data_arr[:,4,:])/16777216.
    UH52_mat=np.transpose(data_arr[:,5,:])/16777216.
    
    BB = np.transpose(data_arr[:,7,:])/16777216.
    
    lim_UH32=np.transpose(data_arr[:-5,10,:i_max])/16777216. 
    std_UH32=np.transpose(data_arr[:-5,9,:i_max])/16777216. 
    CUT_32=-np.transpose(data_arr[:-5,15,:i_max])
    CUT_52=-np.transpose(data_arr[:-5,14,:i_max])
    CUT_BB=-np.transpose(data_arr[:-5,8,:i_max])
    mean_UH32=np.transpose(data_arr[:-5,12,:i_max])/16777216. 
    lim_UH52=np.transpose(data_arr[:-5,11,:i_max])/16777216. 
    mean4_UH32=np.transpose(data_arr[:-5,4,:i_max])/16777216. 
    mean4_UH52=np.transpose(data_arr[:-5,5,:i_max])/16777216. 
    lim_BB=np.transpose(data_arr[:-5,13,:i_max])/16777216. 
    mean4_BB=np.transpose(data_arr[:-5,7,:i_max])/16777216. 
    
    
    UHorBB=-np.transpose(data_arr[:-5,31,:i_max])
    
    b_fen=[25,260]    
    # Clock_control=np.transpose(data_Con[:,0,:])
    # i_max_Con=np.argmax(np.max(Clock_control,axis=1))
    # Output2_mat_con=np.transpose(data_Con[:,2,:i_max_Con])
    # C1F0_mat_con=np.transpose(data_Con[:,18,:i_max_Con])
    
    Output2_IGT_mat=Mat_Pr_Gen(Output2_mat, amplitude*fit[0]+fit[1], -16000., 14000.)
    
    Output_IGT_mat_R=Output2_IGT_mat[:,b_fen[0]:b_fen[1]]
    arrayPr=np.max(Output_IGT_mat_R,axis=1)
    # Output2_IGT_mat_con=Mat_Pr_Gen(Output2_mat_con, amplitude*fit[0]+fit[1], -16000., 14000.)
    
    #%%

    trajet_retour = trajet

    
    path(trajet_retour)
    delai = []
    deb, fin = 40, 280
    convergence = []
    convergence_nosat = []
    cut_pulse = []
    n_cut = 0    
    pulses_cuted = []
    time_cuted = []
    n_satur = 0
    for j in range(np.shape(CUT_BB)[0]):
        convergence.append(np.max(Output2_IGT_mat[j][deb:fin]))
        val_end = Output2_IGT_mat[j][260]
        if convergence[-1] >= 0.999*np.max(Output2_IGT_mat[:][deb:fin]):
            n_satur += 1
        else :
            convergence_nosat.append(np.max(Output2_IGT_mat[j][deb:fin]))
        if val_end <= 0.2*convergence[-1]:
            cut_pulse.append(True)
            n_cut+=1
            pulses_cuted.append(j)
            seuil = (np.max(Output2_IGT_mat[j][deb:fin])-np.min(Output2_IGT_mat[j][deb:fin]))*0.1 + np.min(Output2_IGT_mat[j][deb:fin])
            for k in range(45, 280):
                if Output2_IGT_mat[j][k] <= seuil and np.max(Output2_IGT_mat[j][deb : k]) == np.max(Output2_IGT_mat[j][deb : fin]):
                    # plt.plot(Output2_IGT_mat[j][deb:fin])
                    # plt.plot([0,fin-deb],[seuil,seuil])
                    # sys.exit()
                    time_cuted.append(k)
                    break
            
                    
        else :
            cut_pulse.append(False)
        
    print(" Il ya eu {} coupure : {}%".format(n_cut, np.round(n_cut/np.shape(CUT_BB)[0]*100, decimals = 2)))
    print(" Il ya eu {} saturation : {}%".format(n_satur, np.round(n_satur/np.shape(CUT_BB)[0]*100, decimals = 2)))  
    plt.clf()
    plt.title("Différentes pressions de stabilisation",fontsize=30, fontweight = 'bold')
    plt.xlabel('Différents pulses', fontsize=20)
    plt.ylabel('Pression (kPa)', fontsize=20)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.plot(convergence, c = 'blue', label = "Valeur de stabilisation de pression")
    for ind_p, cut_v in enumerate(cut_pulse):
        if cut_v:
            plt.scatter(ind_p,convergence[ind_p], marker = 'x' , color = 'red', s = 100, linewidths = 4, label = "Pulse coupé")
    plt.savefig(trajet + "pression de stabilisation.png",bbox_inches='tight')
    plt.clf()
    plt.hist(convergence,bins =120)
    plt.xlabel('Pression (kPa)', fontsize=20)
    plt.ylabel('Nombres de pulses', fontsize=20)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.savefig(trajet + "histo pression.png",bbox_inches='tight')
    plt.clf()
    plt.hist(convergence_nosat,bins =120)
    plt.xlabel('Pression (kPa)', fontsize=20)
    plt.ylabel('Nombres de pulses', fontsize=20)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.savefig(trajet + "histo pression no sat.png",bbox_inches='tight')
    plt.clf()
    plt.hist(pulses_cuted,bins =8)
    plt.xlabel('Numéro de pulses', fontsize=20)
    plt.ylabel('Nombres de pulses', fontsize=20)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.savefig(trajet + "histo pulses cuted.png",bbox_inches='tight')
    plt.clf()
    plt.hist(time_cuted,bins =8)
    plt.xlabel('fenêtre de coupure', fontsize=20)
    plt.ylabel('Nombres de pulses', fontsize=20)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.savefig(trajet + "histo window du cut.png",bbox_inches='tight')
    plt.clf()
    
    # plt.savefig(trajet_retour + "pulse_{}.png".format(j),bbox_inches='tight')
     
plt.close("all")
sys.exit()
#%%
plt.close('all')
#TEST de la map de pression et somme des events :
CUT_UH1 =  CUT_32   
CUT_UH2 =  CUT_52   
fsize_single=[14,12]

def plot_events(arrayPr,Pr_max,matComp1,matComp2,matComp3,arrayComp,titles,t_sequence,b_fen,born,fsize):

    t=np.linspace(0,t_sequence,len(arrayPr))
    
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1,5, gridspec_kw={'width_ratios': [1,2,2,2,1]}, figsize=(fsize[0],fsize[1]))
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    ax1.set_ylim(t_sequence,0)
    # ax1.set_xlim(1.1*Pr_max,0)
    ax1.fill_betweenx(t,arrayPr)
    ax1.plot([Pr_max,Pr_max],[t[0],t[-1]],color='r')
    ax1.set_ylabel('SEQUENCE TIME (s)')
    ax1.set_xlabel(titles[0])
    ax1.grid('on')
    im = ax2.imshow(matComp1,aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax2.set_title(titles[1][0])
#    ax2.scatter()
    #xlabel_choice(ax2)    
    im = ax3.imshow(matComp2,aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax3.set_title(titles[1][1])
#    ax2.scatter()
    #xlabel_choice(ax2)    
    im = ax4.imshow(matComp3,aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax4.set_title(titles[1][2])
#    ax2.scatter()
    #xlabel_choice(ax2)
    ax5.plot((arrayComp),t)
    # ax3.plot(20*np.log10(array1),t)
    ax5.invert_yaxis()
    ax5.set_ylim(t_sequence,0)
    # ax3.set_xlim(0,15)
    ax5.set_xlabel(titles[2])
    ax5.grid('on')
    # fig.colorbar(im,ax=ax2)
    # fig.tight_layout()
    # ##
    return fig, (ax1, ax2, ax3, ax4, ax5)


b_abs_lower=[0, 9e-6 * 1.3]
titles=['PNP (MPa)',['Events UH1', 'Events UH2', 'Events BB'],'<BN>']
matComp1=np.array(CUT_UH1)
matComp2=np.array(CUT_UH2)
matComp3=np.array(CUT_BB)
arrayComp=np.mean(matComp1,axis=1) + np.mean(matComp2,axis=1) +np.mean(matComp2,axis=1)
fig, (ax1, ax2, ax3, ax4, ax5)=plot_events(arrayPr,Pr_max,matComp1,matComp2,matComp3,arrayComp,titles,t_sequence,b_fen,b_abs_lower,fsize_single)



#%%
start,end = 0,393021



print("\nLoading pulses......")
data = np.transpose(np.load(tra+'\\Bubbles.npy'))[:-3] #

print("Loading done!")
test_exp=experiment(25000000,1500000,start=start,end=end, size_decal = 1271)#1271
sys.exit()
print("\nAdding pulses exp")
test_exp.add_pulses(data, spacer =100e3)
print("Adding done!")
x_press = [indi for indi in range(test_exp.pulses[0].n_window)]

#%%
indice = 8
plt.figure(figsize=(15,8))
plt.title("BB", fontsize=20)
comparaison_pulse = test_exp.pulses[indice].indice_BB_w
comparaison_lissée = glissante(comparaison_pulse, 4)
# comparaison_pulse1 = test_exp.pulses[indice].indice_harm_w[2]
# comparaison_pulse2 = test_exp.pulses[indice].indice_harm_w[3]
comparaison_mat = BB[indice]
plt.plot(comparaison_lissée, 'magenta',label = 'a partir pulse')    
#plt.plot((np.transpose(Output2_IGT_mat[indice])+50)/350, c = 'black',label = 'commande IGT')# * np.max(C1F0_mat[indice])/np.max(Output2_IGT_mat[indice])

comparaison_mat = BB[indice+i]
plt.plot(np.transpose(comparaison_mat)* np.max(comparaison_lissée)/np.max(comparaison_mat), label = 'a partir array {}'.format(i))    

# plt.plot(comparaison_pulse1,label = 'a partir pulse2')
# plt.plot(comparaison_pulse2,label = 'a partir pulse3')
plt.legend()
plt.show()

# print("plot cartes de pression raw")
# n_pulse = np.shape(data)[0]
# fit = [1,0]
# nbit= [1,n_pulse]
# test_exp.plot_windowed(traj,nbit,fit,20,99,1,legend, ramp = True)
# test_exp.plot_windowed(traj,nbit,fit,40,100,1,legend, ramp = True)

#%%
time_dec = []
for j in range(np.shape(UHorBB)[0]):
    for k in range(np.shape(UHorBB)[1]):
        if UHorBB[j,k] ==1 :
            time_dec.append(k*50-3000)
            break
plt.figure(figsize=(15,8))
plt.title("Différents temps de déclenchement de la coupure",fontsize=30, fontweight = 'bold')
plt.xlabel('Différents pulses', fontsize=20)
plt.ylabel('Temps de déclenchement (µs)', fontsize=20)
plt.xticks(fontsize= 15)
plt.yticks(fontsize= 15)
plt.plot(time_dec, c = 'blue', label = 'UH 5f0/2')
plt.show()  
plt.figure(figsize=(15,8))
plt.hist(time_dec,bins =120)
plt.xlabel('Temps de déclenchement (µs)', fontsize=20)
plt.ylabel('Nombres de pulses', fontsize=20)
plt.xticks(fontsize= 15)
plt.yticks(fontsize= 15)
plt.show()

#%%
plt.close('all') 
deb = 7    
indice = 0
xabs = np.arange(len(data[10]))/25000.
plt.figure(figsize=(15,8))
plt.title("Commande IGT", fontsize=20)
plt.plot(xabs,data[20] )

plt.xlabel('Temps (ms)', fontsize=20, fontweight = 'bold')
plt.xticks(fontsize= 15)

#%%
trajet_retour = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Analyse_VItro_control\Paul bulles ampli\\temps de rep'
trajet_retour = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Analyse_VItro_control\Paul bulles ampli\\'
trajet_retour +='temps de rep 1 bloc - 6 feed _ pour slides\\'
path(trajet_retour)
delai = []
deb, fin = 20, 220
indiceee = []
for decal in [ 0]:
    chem = trajet_retour + "decal_new={}\\".format(decal)
    path(chem)
    for j in range(33,53):
        plt.figure(figsize=(18,10))
        event_ind = 0
        stop_ramp = 0
        #print(CUT_BB[j])
        
        for ind in range(67,160):
            if (CUT_BB[j][ind] == 1 or CUT_52[j][ind] == 1 or CUT_32[j][ind] == 1) and event_ind == 0 :
                event_ind  = ind
                #print("propre")
            if (Output2_IGT_mat[j+decal][ind] == Output2_IGT_mat[j+decal][ind-1]) and stop_ramp == 0 :
                stop_ramp = ind
                #print("salutation")
                delai.append(stop_ramp-event_ind)
                indiceee.append(j)
                break
        commande = Output2_IGT_mat[j+decal][deb:fin] 
        commande_0 = commande - np.min(commande)
        plt.clf()
        
        plt.plot(C1F0_mat[j+decal][deb:fin] /np.max(C1F0_mat[2:200][deb:fin]- np.min(C1F0_mat)), c = 'orange', marker = '+')
        plt.title("Différents cut pulse {}  ".format(j), fontsize=20)
        plt.plot(commande_0 /np.max(Output2_IGT_mat[2:200][deb:fin]- np.min(commande)), c = 'blue', marker = '+')
        plt.plot(np.transpose(CUT_BB[j][deb:fin]) , c = 'red', label = 'BB')
        plt.plot(np.transpose(CUT_52[j][deb:fin]) , c = 'green', label = 'UH 52')
        plt.plot(np.transpose(CUT_32[j][deb:fin]) , c = 'cyan', label = 'UH 32')
        plt.plot(np.transpose(UHorBB[j][deb:fin]) , c = 'black', label = 'UH or BB events cut')
        # if event_ind > 0:
        #     plt.scatter(event_ind-deb,np.max(commande_0 /np.max(Output2_IGT_mat[2:200][deb:fin])), marker = 'x' , color = 'k', s = 100, linewidths = 3, label = "premiere detection d'event")
        # if stop_ramp > 0:
        #     plt.scatter(stop_ramp-deb,np.max(commande_0 /np.max(Output2_IGT_mat[2:200][deb:fin])), marker = 'x' , color = 'red', s = 150, linewidths = 2, label  = "fin de la rampe")
        plt.legend()
        #plt.xlim([event_ind-deb-6,event_ind-deb+1])
        
        #plt.ylim([-0.05,1.05])
        plt.grid()
        plt.show()
        plt.savefig(chem + "pulse_{}.png".format(j),bbox_inches='tight')
            
sys.exit()        
plt.clf()
plt.plot(indiceee ,delai, marker = 'x')
plt.title("Temps de réponse en fenêtres de l'algo", fontsize=20) 

#%%
plt.close("all")
trajet_retour = trajet
plt.figure(figsize=(16,9))

path(trajet_retour)
delai = []
deb, fin = 50, 280
for n_cut in [3,5]:
    indiceee = []
    n_pulsess=0
    
    for decal in [ 0]:
        for j in range(1,1120):
            BB_cons = 0
            UH_cons = 0
            
            
            for elt in  CUT_BB[j]:
                if elt ==1 : 
                    BB_cons+=1
                else :
                    BB_cons = 0
                if BB_cons >=n_cut:
                    n_pulsess+=1
                    break
            if BB_cons<n_cut:
                for k_w in  range(len(CUT_52[j])):
                    if CUT_52[j][k_w] ==1 or CUT_32[j][k_w] ==1 : 
                        UH_cons+=1
                    else :
                        UH_cons = 0
                    if UH_cons >=n_cut:
                        n_pulsess+=1
                        break
            if n_cut == 3:
                plt.clf()
                plt.plot(C1F0_mat[j+decal][deb:fin] /np.max(C1F0_mat[2:200][deb:fin]- np.min(C1F0_mat)), c = 'orange', marker = '+')
                plt.title("Différents cut pulse {}  ".format(j), fontsize=20)
                #plt.plot(commande_0 /np.max(Output2_IGT_mat[2:200][deb:fin]- np.min(commande)), c = 'blue', marker = '+')
                plt.plot(np.transpose(CUT_BB[j][deb:fin]) , c = 'red', label = 'BB')
                plt.plot(np.transpose(CUT_52[j][deb:fin]) , c = 'green', label = 'UH 52')
                plt.plot(np.transpose(CUT_32[j][deb:fin]) , c = 'cyan', label = 'UH 32')
                plt.plot(np.transpose(UHorBB[j][deb:fin]) , c = 'black', label = 'UH or BB events cut')
                # if event_ind > 0:
                #     plt.scatter(event_ind-deb,np.max(commande_0 /np.max(Output2_IGT_mat[2:200][deb:fin])), marker = 'x' , color = 'k', s = 100, linewidths = 3, label = "premiere detection d'event")
                # if stop_ramp > 0:
                #     plt.scatter(stop_ramp-deb,np.max(commande_0 /np.max(Output2_IGT_mat[2:200][deb:fin])), marker = 'x' , color = 'red', s = 150, linewidths = 2, label  = "fin de la rampe")
                plt.legend()
                #plt.xlim([event_ind-deb-6,event_ind-deb+1])
                
                #plt.ylim([-0.05,1.05])
                plt.grid()
                plt.savefig(trajet_retour + "pulse_{}.png".format(j),bbox_inches='tight')
    print(" pour {} consécutifs, il y a {} pulses qui auraient coupés".format(n_cut,n_pulsess))    

#%%
plt.close("all")
trajet_retour = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Analyse_VItro_control\Paul bulles ampli\\fin de array'
path(trajet_retour)
delai = []
deb, fin = 50, 100
plt.figure(figsize=(18,10))
chem = trajet_retour + "decal={}\\".format(decal)
path(chem)
res = []
for j in range(2,150):
    event_ind = 0
    #print(CUT_BB[j])
    ind_res = [-10,-10,-10]
    trouve = [False, False, False]
    for ind in range(100,0,-1):
        if (lim_BB[j][ind-1] != lim_BB[j][ind]) and not(trouve[0]) :
            ind_res[0] = ind
            trouve[0] = True
        if (lim_UH52[j][ind-1] != lim_UH52[j][ind]) and not(trouve[2]) :
            ind_res[2] = ind
            trouve[2] = True
        if (lim_UH32[j][ind-1] != lim_UH32[j][ind]) and not(trouve[1]) :
            ind_res[1] = ind
            trouve[1] = True
        if trouve == [True, True, True]:
            break
    res.append(ind_res)
res = np.array(res)

plt.clf()
plt.plot(res[:,0], c = 'red', label = 'BB')
plt.plot(res[:,1], c = 'green', label = 'UH 32')
plt.plot(res[:,2], c = 'cyan', label = 'UH 52')
plt.legend()
# plt.title("fin de remplissage de l'array", fontsize=20) 
#%%
# for j in range(0,n_pulse,5):
#     print(j)
#     test_exp.pulses[j].plot(traji+"spec_{}".format(j))
plt.close('all') 
deb = 7    
indice = 0
commande = Output2_IGT_mat[indice+1][deb:]
commande_0 = commande - np.min(commande)


plt.figure(figsize=(15,8))

plt.title("Commande IGT et F0", fontsize=20)
plt.plot(commande_0 * np.max(C1F0_mat[indice+1][deb:])/np.max(commande_0), c = 'black')
plt.plot(np.transpose(C1F0_mat[indice+1][deb:]), c = 'blue')    
plt.figure(figsize=(15,8))
plt.title("UH 3/2f0 et sa limite ", fontsize=20)
plt.plot(commande_0 /np.max(commande_0), c = 'blue')
plt.plot(np.transpose(UH32_mat[indice+1][deb:]/np.max(UH32_mat[indice+1][deb:])) , c = 'green')
plt.plot(np.transpose(CUT_32[indice+1][deb:]) , c = 'black')
# plt.plot(np.transpose(mean_UH32[indice+1][deb:])/np.max(UH32_mat[indice+1][deb:]) , c = 'green')
# plt.plot(np.transpose(std_UH32[indice+1][deb:])/np.max(UH32_mat[indice+1][deb:]) , c = 'green')
plt.plot(np.transpose(lim_UH32[indice+1][deb:])/np.max(UH32_mat[indice+1][deb:]), c = 'red')
plt.figure(figsize=(15,8))
plt.title("UH 5/2f0 et sa limite ", fontsize=20)
plt.plot(commande_0 /np.max(commande_0), c = 'blue')
plt.plot(np.transpose(UH52_mat[indice+1][deb:]/np.max(UH52_mat[indice+1][deb:])) , c = 'green')
plt.plot(np.transpose(CUT_52[indice+1][deb:]) , c = 'black')
plt.plot(np.transpose(lim_UH52[indice+1][deb:]/np.max(mean4_UH52[indice+1][deb:])), c = 'red')
plt.show()
plt.figure(figsize=(15,8))
plt.title("BB et sa limite ", fontsize=20)
plt.plot(commande_0 /np.max(commande_0), c = 'blue')
plt.plot(np.transpose(BB[indice+1][deb:]/np.max(BB[indice+1][deb:])) , c = 'red')
plt.plot(np.transpose(CUT_BB[indice+1][deb:]) , c = 'black')
plt.plot(np.transpose(lim_BB[indice+1][deb:]/np.max(mean4_BB[indice+1][deb:])), c = 'red')
plt.show()
plt.figure(figsize=(15,8))
plt.title("Différentes limites dépassées",fontsize=30, fontweight = 'bold')
plt.plot(commande_0 /np.max(commande_0), c = 'blue')
plt.plot(np.transpose(CUT_BB[indice+1][deb:]) , c = 'red', label = 'BB')
plt.plot(np.transpose(CUT_52[indice+1][deb:]) , c = 'green', label = 'UH 52')
plt.plot(np.transpose(CUT_32[indice+1][deb:]) , c = 'cyan', label = 'UH 32')
plt.plot(np.transpose(UHorBB[indice+1][deb:]) , c = 'black', label = 'OR')
plt.legend()
plt.show()
plt.figure(figsize=(15,8))
plt.title("Différents indices dans un pulse dans avec des MBs 1:100",fontsize=30, fontweight = 'bold')
plt.plot(commande_0 /np.max(commande_0), c = 'black')
UH_52_norm = (np.transpose(UH52_mat[indice+1][deb:]/np.max(UH52_mat[indice+1][deb:])))
UH_32_norm = (np.transpose(UH32_mat[indice+1][deb:]/np.max(UH52_mat[indice+1][deb:])))
BB_norm = np.transpose(BB[indice+1][deb:]/np.max(BB[indice+1][deb:]))
plt.plot(UH_52_norm - np.mean(UH_52_norm[:50]) , c = 'green', label = 'UH 5f0/2')
plt.plot(UH_32_norm - np.mean(UH_32_norm[:50]), c = 'darkcyan', label = 'UH 3f0/2')
plt.plot(BB_norm - np.mean(BB_norm[:50]), c = 'red', label = 'Broadband')
plt.legend(fontsize=20)
plt.ylim([-0.25,1.03])
plt.show()
plt.figure(figsize=(15,8))
plt.title("Différents indices dans un pulse dans avec des MBs 1:100  ",fontsize=30, fontweight = 'bold')
plt.plot(commande_0 /np.max(commande_0), c = 'black')
plt.plot(UH_52_norm - np.mean(UH_52_norm[:50]) , c = 'green', label = 'UH 5f0/2')
plt.plot(lim_UH52[indice+1][deb:]/np.max(UH52_mat[indice+1][deb:]) - np.mean(UH_52_norm[:50]) , c = 'green')
plt.plot(UH_32_norm - np.mean(UH_32_norm[:50]), c = 'darkcyan', label = 'UH 3f0/2')
plt.plot(lim_UH32[indice+1][deb:]/np.max(UH52_mat[indice+1][deb:]) - np.mean(UH_32_norm[:50]) , c = 'darkcyan')
plt.plot(BB_norm - np.mean(BB_norm[:50]), c = 'red', label = 'Broadband')
plt.plot(lim_BB[indice+1][deb:]/np.max(BB[indice+1][deb:]) - np.mean(BB_norm[:50]) , c = 'red')
plt.legend(fontsize=20)
plt.ylim([-0.25,1.03])
plt.show()
# plt.figure(figsize=(15,8))
# plt.title("ratio de STD ", fontsize=20)
# plt.plot(np.transpose((BB[indice+1][deb:-7]-np.mean(BB[indice+1][0:40]))/np.std(BB[indice+1][0:40])) , c = 'red')
# plt.plot(np.transpose((UH32_mat[indice+1][deb:-7]-np.mean(UH32_mat[indice+1][0:40]))/np.std(UH32_mat[indice+1][0:40])) , c = 'green')
# plt.plot(np.transpose((UH52_mat[indice+1][deb:-7]-np.mean(UH52_mat[indice+1][0:40]))/np.std(UH52_mat[indice+1][0:40])) , c = 'darkcyan')
# plt.show()
#%%  
indice = 20
a = []
b = []
c = []
d , e , f =[] , [], []
for k in range(100):
    arr = np.transpose(lim_UH32[indice+k][deb:])
    arr_bis = np.transpose(lim_UH52[indice+k][deb:])
    #arr_val = np.transpose(lim_BB[indice+k][deb:])*5
    arr_d = np.transpose(UH32_mat[indice+k][deb:])
    arr_e = np.transpose(UH52_mat[indice+k][deb:])
    arr_f = np.transpose(BB[indice+k][deb:])
    for j in range(len(arr))  :
        a += [arr[j]]
        #b += [arr_val[j]]
        c += [arr_bis[j]]
        d += [arr_d[j]]
        e += [arr_e[j]]
        f += [arr_f[j]]
plt.plot(a,label = 'UH32')
plt.plot(f,label = 'BB')
plt.plot(c,label = 'UH52')
plt.plot(d,label = 'UH32_val')
plt.plot(e,label = 'UH52_val')
plt.plot([0, len(a)],[np.mean(d+e)+2*np.std(d+e),np.mean(d+e)+2*np.std(d+e)], label = "VAL moy all")
print("std UH 32= {}".format(np.std(d)))
print("std UH 52= {}".format(np.std(e)))
print("std bb= {}".format(np.std(f)))

print("\n\nratio UH 32= {}".format(np.std(d)/np.std(f)))
print("ratio UH 52= {}".format(np.std(e)/np.std(f)))

#plt.plot(f,label = 'BB_val')
plt.grid()
#plt.ylim(4e-5,6e-5)
plt.legend()
plt.show()
#plt.plot(np.concatenate(np.transpose(lim_UH32[indice+1][deb:]),np.transpose(lim_UH32[indice+2][deb:])))



#%%
sys.exit()


print("\nPLot des différents indices")
print(legend[i])
#test_exp.plot_indice_RAMP("treatment",traj,x_press, vivo= True, all_plot= True)
for std in range(5,10,2):
        print("STD = ", std)
        test_exp.plot_indice_RAMP_std(legend[i],traj+"MOY_std_{}\\".format(std),still_wind = 25, std_tresh = std, true_harm = False, plot_true = False)
#%%   

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

