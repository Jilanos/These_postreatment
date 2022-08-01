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
import random as rd
import sys
import gc
pl.close('all')


dossier="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_6\\bubbles_25_100"
data = np.load(dossier+'\\data_order.npy')
shape = np.shape(data)
data = np.reshape(data, (shape[0]//30,30,shape[1]))
amp = np.load(dossier+'\\amp_order.npy')


#%%
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

p = 0.4
rep = 5
start = 5
pas_0 = 3
seuil = 40
range_bit = [0,99]


        
val_dic = {}
for press in range(range_bit[0],range_bit[1]+1):
    val_dic[press] = []

mini_doss="SIMU_"+time.strftime("%Y_%m_%d_%H_%M_%S")+"_p_{}_rep_{}_seuil_{}_pas_{}\\".format(p,rep,seuil,pas_0)
folder='C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_8\\'+mini_doss
folder="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_8\\"+mini_doss
if not os.path.exists(folder):
    os.makedirs(folder)
folder_plot = folder + "\\plot"
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)


temps=60
Signal_size=32 ## 32 bien pour un pulse de 10 ms
Elts_par_pulse=Signal_size*8188



test_exp=experiment(fe,freq,start=0,end=-1)
    
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

Data_stock = np.zeros(((bit_max)*average,Elts_par_pulse))



t1= time.time()
t0=t1
while time.time()-t0<temps:
    if time.time()-t0 <temps/10.:
        pas = max(1,int(pas_0/1.))
    elif time.time()-t0 <temps*3/10.:
        pas = max(1,int(pas_0/3.))
    else :
        pas = max(1,int(pas_0/10.))

        
    if len(amp_stock)==0:
            pression, err= start, False
    else :
        pression, err= algo(amp_stock, val_stock, rep, seuil, pas, p, val_dic)
    #print(pression)
    if err:
        err_stock.append([indice,pression])
    courbe_Y = []
    courbe_X = []
    for i in val_dic:
        if len(val_dic[i])>0:
            courbe_Y.append(moyenne_histo(val_dic[i],p))
            courbe_X.append(i)
    if len(courbe_Y)==0:
        maxi.append(0)
    else:
        courbe_Y = list(filtr(courbe_X,courbe_Y))
        
        maxi.append(courbe_X[courbe_Y.index(max(courbe_Y))])
    if pression>range_bit[1]:
        pression=range_bit[1]-10
    if pression<range_bit[0]:
        pression=range_bit[0]+10
    print("\nAmplitude : {}".format(int(pression)))
    Data=np.array(data[pression,rd.randint(0,29)])
    
    
    test_exp.add_pulse(Data, spacer =100e3)
    amp_stock.append(pression)
    UH = np.mean(test_exp.pulses[-1].indice_Uharm[1:4])
    UH_stock.append(UH)
    BB = np.mean(test_exp.pulses[-1].indice_BB)
    BB_stock.append(BB)
    SC = np.mean(test_exp.pulses[-1].indice_harm[1:4])
    SC_stock.append(SC)
    val = 3* UH - 1.5*20* BB
    val_dic[pression].append(val)
    #Data_stock[indice] = Data
    indice += 1
    val_stock.append(val)
    
    
    #print(str(Offset_OUTPUT.read()))
    #print(str(Ampl_Max_OUTPUT.read()))
    t2 = time.time()
    print("Durée  : {}ms".format(int((t2-t1)*1000)))
    t1=t2


        
print("nombre de tirs : {}".format(indice+1))
print('saving data')
#np.save((folder+'\\data.npy'), Data_stock)
np.save((folder+'\\amp.npy'), amp_stock)
text_file = open(folder+"\\"+"param.txt",'w')
text_file.write("PARAMETRES :\n----------\n\nFPGA\n\nSampling frequency : {}\nNombre de samples : {}\nNombres de pulses : {} \nDuree de pulses (us): {} \nDuree de pause (us): {} \nbits igt de tir : {} - {} \nFrequence : {} \nMicrobulles : {}\np : {}\npas : {}\nSeuil : {}\nrep : {}\n".format(fe,Elts_par_pulse,average,PL,pause,2,bit_max,freq,bubbles,p,pas_0,seuil,rep))
text_file.close()

courbe_Y = []
courbe_Y_moy = []
courbe_X = []
for i in val_dic:
    if len(val_dic[i])>0:
        courbe_Y.append(moyenne_histo(val_dic[i],p))
        courbe_Y_moy.append(np.mean(val_dic[i]))
        courbe_X.append(i)
        
        
#%
pl.figure(figsize=(20,10))
pl.plot(val_stock)
pl.title("values obtained")
pl.ylabel("amplitude")
pl.xlabel("PULSES)")
pl.tight_layout()
pl.savefig(folder+"\\plot_VALUES.png")

pl.clf()
pl.plot(courbe_X,filtr(courbe_X,courbe_Y),label = "courbe regardée algo")
pl.plot(courbe_X,courbe_Y_moy,label = "courbe moyennée total")
pl.title("Courbe lissée des valeurs pour chaque bit")
plt.legend()
pl.ylabel("amlitude")
pl.xlabel("bit)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_UH_BB_par_bit.png")

pl.clf()
pl.plot(UH_stock,c='green')
pl.title("Ultraharmonics")
pl.ylabel("Magnitude")
pl.xlabel("pulse)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_courbe_UH.png")

pl.clf()
pl.plot(SC_stock,c='blue')
pl.title("Harmonics")
pl.ylabel("Magnitude")
pl.xlabel("pulse)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_courbe_H.png")

pl.clf()
pl.plot(BB_stock,c='red')
pl.title("Inertiel")
pl.ylabel("Magnitude")
pl.xlabel("pulse)")
pl.tight_layout()
pl.grid()
pl.savefig(folder+"\\plot_courbe_BB.png")


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

pl.close("all")
del Data_stock
gc.collect(generation=2)
