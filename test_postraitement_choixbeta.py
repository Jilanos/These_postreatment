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
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.stats import norm, expon, lognorm, gamma, poisson, binom



gc.collect(generation=2)

start, start_detec, end = 55, 62, 150
end_detec = 160


def gaussian_fit(array,nom, xlim_m=None, xlim_p=None):
    mean, std = norm.fit(array)
    x = np.linspace(min(array), max(array), 1000)
    pdf = norm.pdf(x, mean, std)
    plt.figure(figsize=(10,7))
    plt.hist(array, bins=400, density=True, alpha=0.8, color='k', label='Données')
    plt.plot(x, pdf, 'r', linewidth=2, label='Gaussienne ajustée')
    plt.title(nom, fontsize=17, fontweight='bold')
    plt.xlabel('Valeurs', fontsize=15)
    plt.ylabel('Densité de probabilité', fontsize=15)
    if xlim_m is not None and xlim_p is not None:
        plt.xlim(xlim_m,xlim_p)
    plt.legend(fontsize=15)
    plt.show()
    print("\nMoyenne de la gaussienne ajustée:", mean)
    print("Écart type de la gaussienne ajustée:", std)
    return mean, std

def lognormal_fit(array,seuil,nom , xlim_m=None, xlim_p=None):
    params = lognorm.fit(array)
    x = np.linspace(min(array), max(array)+4, 10000)
    pdf = lognorm.pdf(x, *params)
    plt.figure(figsize=(10,7))
    plt.hist(array, bins=400, density=True, alpha=0.75, color='k', label='Données')
    plt.plot(x, pdf, 'r', linewidth=1, label='Loi log-normale ajustée')
    #print("intégrale = ", np.sum(pdf) * (x[1]-x[0]))
    tot, i = 0, 1
    while tot<seuil:
        tot += pdf[-i] * (x[1]-x[0])
        i += 1
    print(nom, np.round(x[-i],decimals=2)) # , " avec ", tot, "à l'indice ", -i)
    plt.fill_between(x[-i:], pdf[-i:], color='g', alpha=1, label='Faux positifs avec un seuil de 10/27000')
    plt.plot([x[-i],x[-i]],[0,0.15], 'g--', linewidth=1, label='Seuil déterminé numériquement')
    plt.title(nom, fontsize=17, fontweight='bold')
    plt.xlabel('Valeurs', fontsize=15)
    plt.ylabel('Densité de probabilité', fontsize=15)
    plt.legend(fontsize=15)
    plt.show()
    if xlim_m is not None and xlim_p is not None:
        plt.xlim(xlim_m,xlim_p)
    #print("\nMoyenne de la loi lognormale ajustée:", mean)
    #print("Écart type de la loi lognormale ajustée:", std)
    return params

def gaussian_fit_double_fit(array,nom, xlim_m=None, xlim_p=None):
    array_plus = array[np.where(array>=0)]
    array_minus = array[np.where(array<0)]
    fake_plus =  np.concatenate((array_plus,-array_plus))
    fake_minus =  np.concatenate((array_minus,-array_minus))
    mean_plus, std_plus = norm.fit(fake_plus)
    mean_minus, std_minus = norm.fit(fake_minus)
    x_minus = np.linspace(xlim_m,0,500)
    x_plus = np.linspace(0,xlim_p,500)
    plt.figure(figsize=(10,7))
    plt.hist(array, bins=400, density=True, alpha=0.8, color='k', label='Données')
    pdf_minus = norm.pdf(x_minus, mean_minus, std_minus)
    pdf_plus = norm.pdf(x_plus, mean_plus,std_plus)  
    plt.plot(x_minus, pdf_minus, 'b', linewidth=2, label='Loi normale ajustée partie négative')
    plt.plot(x_plus, pdf_plus, 'r', linewidth=2, label='Loi normale ajustée partie positive')
    plt.plot([0,0],[0,0.5], 'r--', linewidth=2, label='0')
    plt.title(nom, fontsize=17, fontweight='bold')
    plt.xlabel('Valeurs', fontsize=15)
    plt.ylabel('Densité de probabilité', fontsize=15)
    if xlim_m is not None and xlim_p is not None:
        plt.xlim(xlim_m,xlim_p)
    plt.legend(fontsize=15)
    plt.show()

    print("\nMoyenne de la gaussienne ajustée:", mean_plus)
    print("Écart type de la gaussienne ajustée:", std_plus)
    return mean_plus, std_plus

def Scatter_simple_Paul(ax,BN_norm,UH2_norm,list_born,legend=['',''],list_thr=[]):  
  
    Val0=np.array(BN_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts=Val

    Val0=np.array(UH2_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts=Val

    Val0=np.array(BN_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_2=Val

    Val0=np.array(UH2_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_2=Val
         
    # start with a square Figure
    print("coucou")
    ax.plot(BN_pts,iUH2_pts ,linestyle='None',marker='o',markersize=5,color='steelblue', alpha = 0.2, label = 'Différentes fenêtres')
    #ax.plot(BN_pts_2,iUH2_pts_2,label=legend[1],linestyle='None',marker='o',markersize=5,color='indianred', alpha = 0.2)    
    ax.grid('on')  
    
    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    # ax.set_xlim([bxmin,bxmax])
    # ax.set_ylim([bymin,bymax])
    if len(list_thr)!=0:
        thr_y=list_thr[1]
        thr_x=list_thr[0]
        ax.plot([bxmin,bxmax],[thr_y,thr_y],color='r', label = 'Seuil', linewidth = 3)
        ax.plot([thr_x,thr_x],[bymin,bymax],color='r', linewidth = 3)
        
plt.close('all')
tra = r"F:\data_vitro\EXPE_aout_2023\20230810__PHD_PAUL_long baseline"
trajet = tra + "\\test_new_control"
path(trajet)
name = ["ControlArr.npy"] 
legend = ["test vivo"] 
i = 0


print("\nLoading pulses......")
data_raw = np.load(tra+"\\"+name[i]) #Bubbles #Control
print(np.shape(data_raw)) # = (312, 32, 203)
Y = data_raw #(np.transpose(data_raw[:,:]))

for j in range(10):
    debut, fin = j*100+33, (j+1)*100+33
    X = Y[:,:,debut : fin+1 ]
    print("\nPartie de la baseline nr ", j+1, " sur 10")    

    Clock=np.transpose(X[:,0,:])
    i_max=np.argmax(np.max(Clock,axis=1))
    BN1_mat=np.transpose(X[:,7,:i_max])
    UH1_mat=np.transpose(X[:,4,:i_max])
    UH2_mat=np.transpose(X[:,5,:i_max])

    BN1_mat_std = np.subtract(np.transpose(BN1_mat) , np.mean(BN1_mat[:,0:40],axis = 1))/np.std(BN1_mat[:,0:40],axis = 1)
    UH1_mat_std = np.subtract(np.transpose(UH1_mat) , np.mean(UH1_mat[:,0:40],axis = 1))/np.std(UH1_mat[:,0:40],axis = 1)
    UH2_mat_std = np.subtract(np.transpose(UH2_mat) , np.mean(UH2_mat[:,0:40],axis = 1))/np.std(UH2_mat[:,0:40],axis = 1)
    BN1_mat_std_c = np.array(BN1_mat_std)[2:-7,:]
    UH1_mat_std_c = np.array(UH1_mat_std)[2:-7,:]
    UH2_mat_std_c = np.array(UH2_mat_std)[2:-7,:]

    pas = 0.002
    nombre = 10


    UH_STD_MAt = np.array(UH1_mat_std_c)#[UH1_mat_std_c,UH2_mat_std_c]
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
        if std_y_uh[-(j+1)]>nombre:
            break
        std_opt_UH = val_std_uh[-(j+1)]
        if j == len(val_std_uh)-2:
            print("error bout de liste atteint pour les UH")
    for j in range(len(val_std_BB)-1):
        if std_y_BB[-(j+1)]>nombre:
            break
        std_opt_BB = val_std_BB[-(j+1)]
        if j == len(val_std_BB)-2:
            print("error bout de liste atteint pour le BB")
    print("UH val = {}".format(np.round(std_opt_UH, decimals=2)))
    print("BB val = {}".format(np.round(std_opt_BB, decimals=2)))


    params_UH1 = lognormal_fit(UH1_mat_std_c.flatten(),10./len(UH1_mat_std_c.flatten()),"UH mod = " , -3, 8)#"Indice des ultra-harmoniques normalisé et sa modélisation"
    params_UBB = lognormal_fit(BN1_mat_std_c.flatten(),10./len(BN1_mat_std_c.flatten()), "BB mod = ", -3, 8)#"Indice du BB normalisé et sa modélisation"





