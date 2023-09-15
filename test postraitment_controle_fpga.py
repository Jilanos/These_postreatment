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
    x = np.linspace(min(array), max(array)+4, 10000)
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

def lognormal_fit(array,nom, xlim_m=None, xlim_p=None):
    params = lognorm.fit(array)
    x = np.linspace(min(array), max(array)+4, 10000)
    pdf = lognorm.pdf(x, *params)
    plt.figure(figsize=(10,7))
    plt.hist(array, bins=400, density=True, alpha=0.75, color='k', label='Données')
    plt.plot(x, pdf, 'r', linewidth=1, label='Loi log-normale ajustée')
    print("intégrale = ", np.sum(pdf) * (x[1]-x[0]))
    seuil = 10/25000.
    tot, i = 0, 1
    while tot<seuil:
        tot += pdf[-i] * (x[1]-x[0])
        i += 1
    print("seuil = ", x[-i], " avec ", tot, "à l'indice ", -i)
    plt.fill_between(x[-i:], pdf[-i:], color='g', alpha=1, label='Faux positifs avec un seuil de 10/25000')
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

plt.close('all')
tra = r"C:\Users\PM263553\Desktop\These\big_projects\in_vitro\data_vivo"
trajet = tra + "\\test_new_control"
path(trajet)
name = ["baseline_souris588_arr.npy"] 
legend = ["test vivo"] 
i = 0


print("\nLoading pulses......")
data_raw = np.load(tra+"\\"+name[i]) #Bubbles #Control
print(np.shape(data_raw)) # = (312, 32, 203)
X = data_raw #(np.transpose(data_raw[:,:]))
print("Loading done!")    

Clock=np.transpose(X[:,0,:])
i_max=np.argmax(np.max(Clock,axis=1))
BN1_mat=np.transpose(X[:,7,:i_max])
UH1_mat=np.transpose(X[:,4,:i_max])
UH2_mat=np.transpose(X[:,5,:i_max])

BB_mat_no_0 = BN1_mat.flatten()[np.where(BN1_mat.flatten()>0)]
UH1_mat_no_0 = UH1_mat.flatten()[np.where(UH1_mat.flatten()>0)]
UH2_mat_no_0 = UH2_mat.flatten()[np.where(UH2_mat.flatten()>0)]



BN1_mat_std = np.subtract(np.transpose(BN1_mat) , np.mean(BN1_mat[:,0:40],axis = 1))/np.std(BN1_mat[:,0:40],axis = 1)
UH1_mat_std = np.subtract(np.transpose(UH1_mat) , np.mean(UH1_mat[:,0:40],axis = 1))/np.std(UH1_mat[:,0:40],axis = 1)
UH2_mat_std = np.subtract(np.transpose(UH2_mat) , np.mean(UH2_mat[:,0:40],axis = 1))/np.std(UH2_mat[:,0:40],axis = 1)
BN1_mat_std_c = np.array(BN1_mat_std)[2:-7,:]
UH1_mat_std_c = np.array(UH1_mat_std)[2:-7,:]
UH2_mat_std_c = np.array(UH2_mat_std)[2:-7,:]


# je veux modéliser mes données par des lois lognormales
# je veux trouver les paramètres de ces gaussiennes du coté positif uniquement
BB_plus = BN1_mat_std_c.flatten()[np.where(BN1_mat_std_c.flatten()>=0)]
UU1_plus = UH1_mat_std_c.flatten()[np.where(UH1_mat_std_c.flatten()>=0)]
UU2_plus = UH2_mat_std_c.flatten()[np.where(UH2_mat_std_c.flatten()>=0)]
BB_minus = -BB_plus
UU1_minus = -UU1_plus
UU2_minus = -UU2_plus


BB = np.concatenate((BB_plus,BB_minus))
UU1 = np.concatenate((UU1_plus,UU1_minus))
UU2 = np.concatenate((UU2_plus,UU2_minus))

data = UH2_mat_std_c.flatten()


# params_BB = gaussian_fit_double_fit(BN1_mat_std_c.flatten(), "double loi normale du BB")
# params_UH1 = gaussian_fit_double_fit(UH1_mat_std_c.flatten(), "double loi normale du UH1")
#params_UH2 = gaussian_fit_double_fit(UH2_mat_std_c.flatten(), "Indice des Ultra-harmoniques normalisé et sa modélisation", -3, 6)

# params_BB = lognormal_fit(BN1_mat_std_c.flatten(), "Log-normale du BB", -5, 5)
# params_UH1 = lognormal_fit(UH1_mat_std_c.flatten(), "Log-normale du UH1", -3, 6)
params_UH2 = lognormal_fit(UH2_mat_std_c.flatten(), "Indice des ultra-harmoniques normalisé et sa modélisation", -3, 6)



# UH1_mat_no_0 = UH1_mat.flatten()[np.where(UH1_mat.flatten()>0)]
# params_UH1_raw = lognormal_fit(UH1_mat_no_0.flatten(), "Log-normale du UH1raw")

sys.exit()


BB_mean, BB_std = gaussian_fit(BB.flatten(), "Gaussienne du BB +")
UU1_mean, UU1_std = gaussian_fit(UU1.flatten(), "Gaussienne du UH1 +")
UU2_mean, UU2_std = gaussian_fit(UU2.flatten(), "Gaussienne du UH2 +")






#on veut afficher les histogrammes des ces valeurs
plt.figure(figsize=(20,10))
plt.subplot(1,3,1)
plt.hist(BN1_mat_std_c.flatten(),bins=100)
plt.title("BN1")
plt.subplot(1,3,2)
plt.hist(UH1_mat_std_c.flatten(),bins=100)
plt.title("UH1")
plt.subplot(1,3,3)
plt.hist(UH2_mat_std_c.flatten(),bins=100)
plt.title("UH2")
plt.show()
plt.savefig(trajet+"\\histo.png", bbox_inches='tight')
