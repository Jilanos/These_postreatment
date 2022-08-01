# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 09:04:44 2021

@author: MIDAS
"""

from nifpga import Session
import time
import numpy as np
import os
import pga
import matplotlib.pyplot as pl
import random
from classes import *
gen = pga.Generator()

pl.close('all')
### DECLARER PARAMETRES
#Duree_Tir = 15 ### Durée de la séquence de tir (en s)

Adjust_Offset = -90 ### Nombre de pas de quantification (maxi 32767 mini -32768)
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

p = 0.4
rep = 5
start = 60
pas_0 = 9
seuil = 180
range_bit = [0,220]


        
val_dic = {}
for press in range(range_bit[0],range_bit[1]+1):
    val_dic[press] = []


mini_doss=time.strftime("%Y_%m_%d_%H_%M_%S")+"_p_{}_rep_{}_seuil_{}_pas_{}\\".format(p,rep,seuil,pas)
folder='C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\iter_8\\'+mini_doss
if not os.path.exists(folder):
    os.makedirs(folder)
folder_plot = folder + "\\plot"
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)


temps=30
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


  

def USshot(sequence):
    #print('triying TO SEND SEQ')
    gen.sendSequence (sequence)
    #print('SEQ. SENT')
    exec_flags = pga.ExecFlag.ASYNC_EXECUTION_RESULT #| pga.ExecFlag.TRIM_RESULTS_10
    gen.executeSequence (1, 0, exec_flags)
    #print('SEQ. EXECUTED')
    gen.waitExecution()

### PRECALCULER PARAMETRES


gen.loadConfig("generator.json")
# Try to auto connect to the serial port
if not gen.autoConnect():
    raise Exception("Can not connect to the PGA board.")
# Print initial state of the generator (temperature, voltage)
print (gen.readExecutionStatus())
# Start the amplifier (required before execution)
gen.enableAmplifier(True) 
gen.selectOutput(pga.Output.EXTERNAL)      



### OUVRIR BITFILE
path='C:\\Users\\MIDAS\\Documents\\LabVIEW Projects\\20210316_Monitoring_Control\\FPGA Bitfiles'
traj = "C:\\Users\\MIDAS\\Desktop\\code_UH_long\\GENE_MOD\\plot\\"
bitfile='p20210316monitor_FPGATarget_20211122FPGAPyth_TKtrFqpJHlY.lvbitx'
#os.chdir(path)

deb_ramp_up = 10 #300 #10
dur = 16200 #653 #16200
deb_ramp_down = 200000 #9300 #235000

with Session(bitfile,"PXI1Slot2") as session:

    ### Reset le bitfile
    print('RESETTING THE BITFILE')
    session.reset()
    time.sleep(5)
    ### Entrer/Ecrire les valeurs de registres
    print('CREATING REGISTERS')

    RUN_INPUT = session.registers['RUN']
    RUN_INPUT.write(False)

    Nb_values_per_Pulse_INPUT = session.registers['Values per Pulse'] #A DEFINIR PLUS TARD
    Nb_values_per_Pulse_INPUT.write(Elts_par_pulse+1)
    
    Adjust_Offset_INPUT = session.registers['Adjust']#A DEFINIR PLUS TARD
    Adjust_Offset_INPUT.write(Adjust_Offset)

    Amplitude_Pulse_INPUT = session.registers['Amplitude Pulse'] #A DEFINIR PLUS TARD
    Amplitude_Pulse_INPUT.write(32700)
    
    Initial_mod_point_INPUT = session.registers['Initial Mod Point']#A DEFINIR PLUS TARD     
    Initial_mod_point_INPUT.write(16000)
    
    pas_ramp_INPUT = session.registers['Pas'] 
    pas_ramp_INPUT.write(1)
    
    dt_Ramp_Up_INPUT = session.registers['Debut ramp up'] #A DEFINIR PLUS TARD
    dt_Ramp_Up_INPUT.write(deb_ramp_up)
    
    fin_Ramp_Up_INPUT = session.registers['Fin ramp up'] #A DEFINIR PLUS TARD
    fin_Ramp_Up_INPUT.write(deb_ramp_up+dur)
    
    dt_Ramp_Down_INPUT = session.registers['Debut ramp down'] #A DEFINIR PLUS TARD
    dt_Ramp_Down_INPUT.write(deb_ramp_down)
    
    fin_Ramp_Down_INPUT = session.registers['Fin ramp down'] #A DEFINIR PLUS TARD
    fin_Ramp_Down_INPUT.write(deb_ramp_down+dur)

    ### OUTPUT ###
    Offset_OUTPUT = session.registers['IO Module mean']
    Ampl_Max_OUTPUT = session.registers['IO Module max']  
    
    time.sleep(1)

    print('CONFIGURING THE DMA')
    DMA= session.fifos['20210527_DMA_Python_Scope_8188els_I16']
    DMA.configure(Elts_par_pulse*10)
    (Data0,N_elms)=DMA.read(0,timeout_ms=-1)
    Data_Clear=DMA.read(N_elms,timeout_ms=-1)
    time.sleep(1)

    ### Lancer le programme - Attendre le démarage de la clock   
    session.run()
    DMA.start()  
    print('INITIALIZING')

    time.sleep(2)
    ### Commencer le compteur
    RUN_INPUT.write(True) 
    time.sleep(1)  
    t1= time.time()
    t0=t1
    while time.time()-t0<temps:
        if time.time()-t0 <temps/10.:
            pas = max(1,int(pas_0/1.))
        elif time.time()-t0 <temps*3/10.:
            pas = max(1,int(pas_0/3.))
        else :
            pas = max(1,int(pas_0/10.))
        try:
            
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
                courbe_Y = filtr(courbe_X,courbe_Y)
                maxi.append(courbe_X[courbe_Y.index(max(courbe_Y))])
            if pression>range_bit[1]:
                pression=range_bit[1]-10
            if pression<range_bit[0]:
                pression=range_bit[0]+10
            print("\nAmplitude : {}".format(int(pression)))
            sequence = [pga.Pulse(int(PL*1e6), 0, int(pression*2), int(freq))]
            USshot(sequence)
            Data0=(DMA.read(Elts_par_pulse,timeout_ms=10000))
            Data=np.array(Data0[0])
            
            
            test_exp.add_pulse(Data, spacer =100e3)
            amp_stock.append(pression)
            UH = np.mean(test_exp.pulses[-1].indice_Uharm[1:4])
            UH_stock.append(UH)
            BB = np.mean(test_exp.pulses[-1].indice_BB)
            BB_stock.append(BB)
            SC = np.mean(test_exp.pulses[-1].indice_harm[1:4])
            SC_stock.append(SC)
            val = 3* UH - 1*20* BB
            val_dic[pression].append(val)
            Data_stock[indice] = Data
            indice += 1
            val_stock.append(val)
            
            
            #print(str(Offset_OUTPUT.read()))
            #print(str(Ampl_Max_OUTPUT.read()))
            repos = (pause*1e-6-time.time()+t1)
            if repos*1000 > 1 :
                time.sleep(repos)
            
            t2 = time.time()
            
            print("Durée  : {}ms".format(int((t2-t1)*1000)))
            t1=t2
        except Exception as why:
            print ("Exception: "+ str(why))

    try:   
        DMA.stop() 
        RUN_INPUT.write(False)
        time.sleep(1)
        session.abort()
        print('FIN DU TIR')
    except Exception as why:
        print ("Exception: "+ str(why))
        
print("nombre de tirs : {}".format(indice+1))
print('saving data')
np.save((folder+'\\data.npy'), Data_stock)
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
        
        
 #%%
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

