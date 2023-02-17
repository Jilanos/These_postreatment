# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 09:17:11 2023

@author: PM263553
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from classes import *
import sys 
import gc

gc.collect(generation=2)


tra11 = "D:\\data_vivo\\20221205_PCDsig_BBBop_M1\\"
tra21 = "D:\\data_vivo\\20221205_PCDsig_BBBop_M2\\"
tra41 = "D:\\data_vivo\\20221205_PCDsig_BBBop_M4\\"
tra51 = "D:\\data_vivo\\20221205_PCDsig_BBBop_M6\\"
tra61 = "D:\\data_vivo\\20221206_PCDsig_BBBop_M8\\"
tra71 = "D:\\data_vivo\\20221206_PCDsig_BBBop_M9\\"
tra81 = "D:\\data_vivo\\20221206_PCDsig_BBBop_M10\\"


tra1 = "D:\\data_vivo\\20230124_PCDsig_M1\\"
tra2 = "D:\\data_vivo\\20230124_PCDsig_M2\\"
tra3 = "D:\\data_vivo\\20230124_PCDsig_M3\\"
tra4 = "D:\\data_vivo\\20230124_PCDsig_M4\\"
tra5 = "D:\\data_vivo\\20230124_PCDsig_M5\\"
tra6 = "D:\\data_vivo\\20230124_PCDsig_M6\\"
tra7 = "D:\\data_vivo\\20230125_PCDsig_M7\\"
tra8 = "D:\\data_vivo\\20230125_PCDsig_M8\\"
tra9 = "D:\\data_vivo\\20230125_PCDsig_M9\\"

tra99 = "D:\\data_vivo\\20221125_PCDsig_CMUT2\\"
trajet = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Analyse_corentin_M1_10_start\\"
trajet = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\Analyse_corentin_+Cmut\\"
trajet = "C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\data_vivo\\TESSSTETTS\\"
path(trajet)
std = 5
# frequencies= np.arange(2048//2) * 25000000/2048
# freq_harm=np.array([1500000.*i for i in range(1,4+2)])
# freq_harm_n= list(map(lambda x: find_nearest(25000000., x), freq_harm))


legend = ["M{}".format(i) for i in [1,2,4,6,8,9,10]] + ["CMUT"] + ["M{}_bis".format(i) for i in [1,2,3,4,5,6,7,8,9,10]] 



#for i, tra in enumerate([tra1, tra2, tra4, tra5, tra6, tra7, tra8]):
for i, tra in enumerate([tra11, tra21, tra41, tra51, tra61, tra71, tra81, tra99, tra1, tra2, tra3, tra4, tra5, tra6, tra7, tra8, tra9]):
    traj = trajet + legend[i]+ "\\"
    path(traj)
    data = np.transpose(np.load(tra+"Bubbles.npy"))[:-8]
    
    
    
    start,end = 0,262000
       
    print("Loading done!")
    test_exp=experiment(25000000,1500000,start=start,end=end)
    print("\nAdding pulses exp")
    test_exp.add_pulses(data, spacer = 100e3)
    
    
    print("Adding done!")
    x_press = [indi for indi in range(test_exp.pulses[0].n_window)]
    
    print("\nPLot des différents indices")
    print(legend[i])
    for std in [5]: #7,8,
        print("STD : {}".format(std))
        test_exp.plot_indice_RAMP_std(legend[i],traj+"std_{}\\".format(std),still_wind = 10, std_tresh = std, plot_true = False)
    sys.exit()
    #test_exp.plot_indice_RAMP("all_pulses",traj,x_press, vivo= True, all_plot=True)


