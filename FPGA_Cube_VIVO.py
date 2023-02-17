# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 09:04:44 2021

@author: MIDAS
"""

from nifpga import Session
import time
from classes import *
import numpy as np
import os
import pga
import matplotlib.pyplot as pl
import random
gen = pga.Generator()

pl.close('all')
### DECLARER PARAMETRES
#Duree_Tir = 15 ### Durée de la séquence de tir (en s)

Adjust_Offset = -90 ### Nombre de pas de quantification (maxi 32767 mini -32768)
Delay_Acq_trigger = 100000 ### Nombre de ticks d'horloge (1 tick equivaut à 10 ns pour la Data Clock)
amp_stock = []
Signal_size=48 ## 32 bien pour un pulse de 10 ms
Elts_par_pulse=Signal_size*8188

PL=15e-3
bit_max=45
repetition_treatment=1200
repetition_baseline = 10
freq=1.5e6
fe = 25e6
bubbles=False
pause=100000

doss='VIVO_Mouse587_221122_{}'.format(bit_max)
folder='D:\\code_UH_long\\GENE_MOD\\VIVO_1\\'+doss
if not os.path.exists(folder):
    os.makedirs(folder)
folder_plot = folder + "\\plot"
if not os.path.exists(folder_plot):
    os.makedirs(folder_plot)




Data_stock_baseline = np.zeros((repetition_baseline,Elts_par_pulse))

Data_stock = np.zeros((repetition_treatment,Elts_par_pulse))


  

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
#bitfile='p20210316monitor_FPGATarget_20211122FPGAPyth_TKtrFqpJHlY.lvbitx'#classique
bitfile='p20210316monitor_FPGATarget_20220831FPGAPyth_wIwDXgOSnXE.lvbitx'#slow ramp

ramp_divider=20
deb_ramp_up = 10 #300 #10
dur = 16200*ramp_divider #653 #16200
deb_ramp_down = 450000 #9300 #235000

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
        
    pas_ramp_INPUT = session.registers['Ramp Divider'] 
    pas_ramp_INPUT.write(ramp_divider)
    
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
    ok = 0
    z=0
    while  ok==0:
        z+=1
        Data_stock_baseline = np.zeros((repetition_baseline,Elts_par_pulse))
        for i in range(repetition_baseline):
            amp = bit_max*2
            try:
                if i%20==0:
                    print("tir : {}".format(i))
                sequence = [pga.Pulse(int(PL*1e6), 0, int(amp), int(freq))]
                USshot(sequence)
                Data0=(DMA.read(Elts_par_pulse,timeout_ms=10000))
                Data=np.array(Data0[0])
                Data_stock_baseline[i] = Data
                repos = (pause*1e-6-time.time()+t1)
                if repos*1000 > 1 :
                    time.sleep(repos)
                t2 = time.time()
                t1=t2
            except Exception as why:
                print ("Exception: "+ str(why))
            
        test_exp=experiment(25000000,1500000,start=0,end=-1)
        print("adding pulses exp")
        test_exp.add_pulses(Data_stock_baseline, spacer =100e3)
        x_press = [indi for indi in range(test_exp.pulses[0].n_window)]
        test_exp.plot_indice_RAMP("baseline_{}".format(z),folder_plot,x_press)
        test_exp.pulses[-1].plot(folder_plot+"baseline_{}".format(z))
        print("Pour refaire une baselin taper 0 sinon taper un autre chiffre pour Valider")
        ok = int(input())
            
    for i in range(repetition_treatment):
        amp = bit_max*2
        try:
            if i%20==0:
                print("tir : {}".format(i))
            sequence = [pga.Pulse(int(PL*1e6), 0, int(amp), int(freq))]
            USshot(sequence)
            Data0=(DMA.read(Elts_par_pulse,timeout_ms=10000))
            Data=np.array(Data0[0])
            Data_stock[i] = Data
            repos = (pause*1e-6-time.time()+t1)
            if repos*1000 > 1 :
                time.sleep(repos)
            t2 = time.time()
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
        
print('saving data')
np.save((folder+'\\data_baseline.npy'), Data_stock_baseline)
np.save((folder+'\\data_treatment.npy'), Data_stock)
          
print('data saved') 
text_file = open(folder+"\\"+"param.txt",'w')
text_file.write("PARAMETRES :\n----------\n\nFPGA\n\nSampling frequency : {}\nNombre de samples : {}\nRepetition treatment : {} \nDuree de pulses (us): {} \nDuree de pause (us): {} \nbits igt de tir : {} - {} \nFrequence : {} \nMicrobulles : {}\n".format(fe,Elts_par_pulse,repetition_treatment,PL,pause,2,bit_max,freq,bubbles))
text_file.close()


#%%

pl.figure(figsize=(20,10))
n=1
Data=Data_stock[n]
x_abs = [i * 1/ fe*1000 for i in range(len(Data))]
pl.plot(x_abs[:],Data[:])
pl.title("Pulse shot amp {} bits".format(bit_max))
pl.ylabel("amplitude")
pl.xlabel("time (ms)")
pl.tight_layout()
pl.savefig(folder_plot+"\\plot_{}.png".format(n))




