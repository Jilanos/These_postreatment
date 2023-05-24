# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 16:05:17 2023

@author: PM263553
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys 


tra = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\'
traj='C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\Analyse_PULSE\\'


fe = 25000000.
f0 = 1500000.

data = np.load(tra + "PULSE_666_40\\data_o.npy")


#%%
pulse_n = 200
y = data[pulse_n]-np.mean(data[pulse_n])


y_s = y[100000:101001]
x = np.arange(len(y_s))/fe

plt.plot(y_s)


#%%
def J(q,m,N,Y,w0):
    W = 
    val = Y.dot(W).dot(Y.transpose())
    ret = N/2*np.log(val)+q*(m+1)*np.log(N)
     
     
     
N = len(y_s)
Y_t = y_s



