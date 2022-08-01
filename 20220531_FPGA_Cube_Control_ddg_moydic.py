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
import optuna
from optuna import Trial
import random as rd
import sys
import gc

gc.collect(generation=2)

dossier="C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter5\\big_data\\\\bubbles_ok_99_181"
dossier="C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_6\\bubbles_50_100"
data = np.load(dossier+'\\data_order.npy')
shape = np.shape(data)
data = np.reshape(data, (shape[0]//30,30,shape[1]))
#amp = np.load(dossier+'\\amp_of.npy')

def descente2(b_n_1, b_n, v_n_1, v_n, lr=2):
    grad = (v_n-v_n_1)/(b_n_1-b_n) #approximation numérique de la dérivée
    if grad>=0:
        return b_n+lr
    else :

        return b_n-lr
    
     
#%%


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


pl.close('all')


for p in [0.1,0.3,0.6]:
    for pas in [3,5,7]:
                
        range_bit = [0,100]
        # p = 0.5
        rep = 5
        start = 5
        # pas = 3
        seuil = 40
        
        val_dic = {}
        for press in range(range_bit[0],range_bit[1]+1):
            val_dic[press] = []
        
        print("p : {}, rep : {}, seuil : {}, pas : {}".format(p,rep,seuil,pas))
        
        doss='resultat_controle_descente_gradiant\\'
        mini_doss=time.strftime("%Y_%m_%d_%H_%M_%S")+"_p_{}_rep_{}_seuil_{}_99\\".format(p,rep,seuil)
        folder='C:\\Users\\PM263553\\Desktop\\'+doss+mini_doss
        folder='C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\EXPE_RESULTS\\'+doss+mini_doss
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        temps=120
        iteration = 1200
        PL=10e-3
        bit_max=180
        average=1
        freq=1.5e6
        fe = 25e6
        bubbles=False
        pause=100000
        Signal_size=32 ## 32 bien pour un pulse de 10 ms
        Elts_par_pulse=Signal_size*8188
        
        
        
        test_exp=experiment(fe,freq,start=0,end=-1)
            
        Data_stock = np.zeros((temps*110,Elts_par_pulse))
        
          
        t1= time.time()
        t0=t1
        indice = 0
        pression = 0
        amp_stock = []
        val_stock = []
        err_stock = []
        maxi = []
        
        while indice<iteration-1:
            
            if len(amp_stock)==0:
                pression, err= start, False
            else :
                pression, err= algo(amp_stock, val_stock, rep, seuil, pas, p)
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
                maxi.append(courbe_X[courbe_Y.index(max(courbe_Y))])
            Data=np.array(data[pression,rd.randint(0,29)])
            test_exp.add_pulse(Data, spacer =100e3)
            amp_stock.append(pression)
            UH = np.mean(test_exp.pulses[-1].indice_Uharm[1:4])
            BB = np.mean(test_exp.pulses[-1].indice_BB)
            val = 3* UH - 20* BB
            val_dic[pression].append(val)
            #Data_stock[indice] = Data
            indice += 1
            val_stock.append(val)
            #print("\nAmplitude : {}, Valeur : {}".format(amp_stock[-1],val_stock[-1]))
            # repos = (pause*1e-6-time.time()+t1)
            # if repos*1000 > 1 :
            #     time.sleep(repos)
            # t2 = time.time()
            # #print("Durée  : {}ms".format(int((t2-t1)*1000)))
            # t1=t2
        
        print("nombre de tirs : {}".format(indice+1))
        print('saving data')
        #np.save((folder+'\\data.npy'), Data_stock)
        np.save((folder+'\\amp.npy'), amp_stock)
            
        text_file = open(folder+"\\"+"param.txt",'w')
        text_file.write("PARAMETRES :\n----------\n\nFPGA\n\nSampling frequency : {}\nNombre de samples : {}\nNombres de pulses : {} \nDuree de pulses (us): {} \nDuree de pause (us): {} \nbits igt de tir : {} - {} \nFrequence : {} \nMicrobulles : {}\n".format(fe,Elts_par_pulse,average,PL,pause,2,bit_max,freq,bubbles))
        text_file.close()
        
        courbe_Y = []
        courbe_Y_moy = []
        courbe_X = []
        for i in val_dic:
            if len(val_dic[i])>0:
                courbe_Y.append(moyenne_histo(val_dic[i],p))
                courbe_Y_moy.append(np.mean(val_dic[i]))
                courbe_X.append(i)
                
        
        
        pl.figure(figsize=(20,10))
        pl.plot(val_stock)
        pl.title("values obtained")
        pl.ylabel("amplitude")
        pl.xlabel("PULSES)")
        pl.tight_layout()
        pl.savefig(folder+"\\plot_VALUES.png")
        
        pl.clf()
        pl.plot(courbe_X,courbe_Y)
        pl.title("Courbe lissée des valeurs pour chaque bit")
        pl.ylabel("amlitude")
        pl.xlabel("bit)")
        pl.tight_layout()
        pl.grid()
        pl.savefig(folder+"\\plot_courbe_lisse.png")
        
        pl.clf()
        pl.plot(courbe_X,courbe_Y_moy)
        pl.title("Courbe lissée des valeurs pour chaque bit")
        pl.ylabel("amlitude")
        pl.xlabel("bit)")
        pl.tight_layout()
        pl.grid()
        pl.savefig(folder+"\\plot_courbe_moy.png")
        
        
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



