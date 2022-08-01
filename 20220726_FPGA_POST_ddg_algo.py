# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 09:04:44 2021

@author: MIDAS
"""

import time
import numpy as np
import os
import matplotlib.pyplot as pl
import random
from classes import *




doss = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_8\\2022_07_26_13_22_01_p_0.5_rep_4_seuil_100_pas_7\\'
folder = doss + "plot_bis\\"
doss_data = doss + "data.npy"
doss_amp = doss + "amp.npy"
data = np.load(doss_data)
amp=np.load(doss_amp)

#%%

path(folder)

Delay_Acq_trigger = 100000 ### Nombre de ticks d'horloge (1 tick equivaut à 10 ns pour la Data Clock)
amp_stock = []
Signal_size=32 ## 32 bien pour un pulse de 10 ms
Elts_par_pulse=Signal_size*8188

PL=10e-3
bit_max=196
average=30
freq=1.5e6
fe = 25e6
bubbles=True
pause=100000
p= 0.5
range_bit = [0,200]



        

val_dic = {}
for press in range(range_bit[0],range_bit[1]+1):
    val_dic[press] = []

def algo(amp,val,rep,seuil,pas,p):
    ptitpas = max([int(pas/2),1])
    if amp.count(amp[-1])<rep:
        return amp[-1], False
    elif max(amp)<seuil:
        if seuil*0.8>max(amp):
            return max(amp)+pas, False
        else :
            return max(amp)+ptitpas, False
    else:
        courbe_Y = []
        courbe_X = []
        for i in val_dic:
            if len(val_dic[i])>0:
                courbe_Y.append(moyenne_histo(val_dic[i],p))
                courbe_X.append(i)
        if abs(courbe_X[courbe_Y.index(max(courbe_Y))]-amp[-1])/float(courbe_X[courbe_Y.index(max(courbe_Y))])> 0.15:
            #print("ouuuuut {}, ecart : {}".format(courbe_X[courbe_Y.index(max(courbe_Y))], abs(courbe_X[courbe_Y.index(max(courbe_Y))]-amp[-1])/float(courbe_X[courbe_Y.index(max(courbe_Y))])))
            return courbe_X[courbe_Y.index(max(courbe_Y))], True
        else:
            amp_1 = amp[-1]
            val_1 = moyenne_histo(val_dic[amp_1],p)
            for j in range(1,len(amp)):
                if amp[-j] != amp_1:
                    amp_2 = amp[-j]
                    break
            
            val_2 = moyenne_histo(val_dic[amp_2],p)
            grad = (val_2-val_1)/(amp_2-amp_1)
            if grad>=0:
                #print("pluuus : last : amp : {}  val{}  // avant : amp : {}  val {}  ".format(amp_1,val_1,amp_2,val_2))
                return amp_1+ptitpas, False
            else :
                #print("moins")
                return amp_1-ptitpas, False


temps=30
Signal_size=32 ## 32 bien pour un pulse de 10 ms
Elts_par_pulse=Signal_size*8188



test_exp=experiment(fe,freq,start=0,end=-1)##
test_exp=experiment(fe,freq,start=15000,end=200000)##
    
Data_stock = np.zeros((temps*11,Elts_par_pulse))

  

indice = 0
pression = 0
amp_stock = []
val_stock = []
err_stock = []
maxi = []

UH_stock = []
BB_stock = []
SC_stock = []

shap = np.shape(data)

for j in range(shap[0]):
    if j%100 ==0:
        print(j)
    Data_el= data[j]
    pression = amp[j]
    test_exp.add_pulse(Data_el, spacer =100e3)
    amp_stock.append(pression)
    UH = np.mean(test_exp.pulses[-1].indice_Uharm[1:4])
    UH_stock.append(UH)
    BB = np.mean(test_exp.pulses[-1].indice_BB)
    BB_stock.append(BB)
    SC = np.mean(test_exp.pulses[-1].indice_harm[1:4])
    SC_stock.append(SC)
    val = 3* UH - 1*20* BB
    val_dic[pression].append(val)
    indice += 1
    val_stock.append(val)



courbe_Y = []
courbe_Y_moy = []
courbe_X = []
for i in val_dic:
    if len(val_dic[i])>0:
        courbe_Y.append(moyenne_histo(val_dic[i],p))
        courbe_Y_moy.append(np.mean(val_dic[i]))
        courbe_X.append(i)
        
        
 #%%
pl.figure(figsize=(20,10))
pl.plot(val_stock)
pl.title("values obtained")
pl.ylabel("amplitude")
pl.xlabel("PULSES)")
pl.tight_layout()
pl.savefig(folder+"\\plot_VALUES_bis.png")

pl.clf()
pl.plot(courbe_X,courbe_Y,label = "courbe regardée algo")
pl.plot(courbe_X,filtr(courbe_X,courbe_Y),label = "courbe gauss")
pl.title("Courbe lissée des valeurs pour chaque bit")
plt.legend()
pl.ylabel("amlitude")
pl.xlabel("bit)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_UH_BB_par_bit_bis.png")

pl.clf()
pl.plot(UH_stock,c='green')
pl.title("Ultraharmonics")
pl.ylabel("Magnitude")
pl.xlabel("pulse)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_courbe_UH_bis.png")

pl.clf()
pl.plot(SC_stock,c='blue')
pl.title("Harmonics")
pl.ylabel("Magnitude")
pl.xlabel("pulse)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_courbe_H_bis.png")

pl.clf()
pl.plot(BB_stock,c='red')
pl.title("Inertiel")
pl.ylabel("Magnitude")
pl.xlabel("pulse)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_courbe_BB_bis.png")


pl.clf()
# pl.figure(figsize=(20,10))
pl.plot(amp_stock,c='black', ls = ":", marker = "o",ms=4)
pl.plot(maxi,c='blue', ls = "--", marker = "+",ms=4)
pl.title("amplitudes shot ")
for i in err_stock:
    pl.scatter(i[0],i[1],s=150,c = 'red', marker = 'x',lw = 1,)
pl.ylabel("amplitude")
pl.xlabel("PULSES)")
pl.grid()
pl.tight_layout()
pl.savefig(folder+"\\plot_AMPLITUDES.png")


#%%
pl.figure(figsize=(20,10))


            
pl.plot(courbe_X,courbe_Y,label = "courbe regardée algo")
pl.plot(courbe_X,filtr(courbe_X,courbe_Y),label = "courbe gauss")
pl.title("Courbe lissée des valeurs pour chaque bit")
plt.legend()
pl.ylabel("amlitude")
pl.xlabel("bit)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_UH_BB_par_bit_ter.png")
