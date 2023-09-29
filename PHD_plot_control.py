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
general_path = 'F:\\data_vitro\\iter_21\\PULSE_0_75\\'



traj_plot = 'C:\\Users\\PM263553\\Desktop\\These\\report\\Manuscript\\img\\3_img\\plot_python\\'
path(traj_plot)

#%%
# start = 50000
window = 2048
decalage = 1280
f0 = 0.1
n = 4
signal = np.zeros((n,window))
echantillon = np.zeros((n,window))
sig_hanninged = np.zeros((n,window))
for j in range(n):
    echantillon[j] = np.array([i + decalage*j for i in range(window)])
    noise = np.random.normal(0, 0.03, window)
    signal[j] = np.sin( echantillon[j]*f0) + noise
for j in range(n):
    sig_hanninged[j] = signal[j]*np.hanning(window)



style = ['--', ':','--', ':','--', ':','--', ':','--', ':','--', ':']
color = ['firebrick', 'limegreen', 'royalblue', 'darkorange', 'lightyellow', 'lightgrey', 'lightpink', 'lightcyan', 'lightsteelblue', 'lightseagreen', 'lightsalmon', 'lightgoldenrodyellow']

plt.figure(figsize=(21,14))
plt.subplot(211)
for j in range(n):
    plt.plot(echantillon[j], signal[j], label='Fenêtre '+str(j+1), color=color[j], linestyle=style[j], linewidth=2)
plt.legend(fontsize=16,framealpha = 1, loc = 'upper right')
plt.title("Signaux temporels fenêtrés des tirs avant et après fenêtrage avec Hanning", fontsize=20, fontweight='bold')
plt.xlabel('Échantillons', fontsize=18)
plt.ylabel('Amplitude', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True)
plt.xlim(0, 5888)



plt.subplot(212)
for j in range(n):
    plt.plot(echantillon[j], sig_hanninged[j], label='Fenêtre '+str(j+1), color=color[j], linestyle=style[j], linewidth=2)
plt.legend(fontsize=16,framealpha = 1, loc = 'upper right')
plt.xlabel('Échantillons', fontsize=18)
plt.ylabel('Amplitude', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True)
plt.xlim(0, 5888)
plt.savefig(traj_plot+'signal_hanning.png', bbox_inches='tight')



data = np.load(general_path+'data_o.npy')




#%%
plt.close("all")

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx
def plot_span(f,ecart,array,c,leg,val):
    for i, f1 in enumerate(f):
        center = find_nearest(array, f1)
        if i == 0:
            plt.axvspan(array[center-ecart], array[center+ecart],color=c,alpha=0.2,lw=0,label=leg)
        else:
            plt.axvspan(array[center-ecart], array[center+ecart],color=c,alpha=0.2,lw=0)
        plt.scatter(array[center], np.mean(val[center-ecart:center+ecart+1]), color=c, s=100, marker='o')

start, end = 60, int(2048/4)
signal_1 = data[1000, 17000:17000+window]
signal_hanninged = signal_1*np.hanning(window)
fft = np.abs(np.fft.fft(signal_hanninged))
freq = np.fft.fftfreq(window, d=1/25000000.)/1e6

mean5 = [0,0]
for i in range(2, len(fft)):
    mean5.append(np.mean(fft[i-2:i+3]))
mean5.append(0)
mean5.append(0)

f0 = 1.5
mini, maxi = np.min(fft[start:end]), np.max(fft[start:end])

plt.figure(figsize=(17,7))

plt.plot(freq[start:end] ,fft[start:end] , color='k', linewidth=1, label='FFT du signal fenêtré')
#plt.plot(freq[start:end] ,mean5[start:end] , color='r', linewidth=2, label='Moyenne glissante sur 5 points')
plot_span([1.5, 3, 4.5, 6],2,freq[start:end],'b',"Harmoniques",fft[start:end])
# plt.plot([1.5,1.5], [mini, maxi], color='b', linewidth=1, linestyle='--')
# plt.plot([3,3], [mini, maxi], color='b', linewidth=1, label='Harmoniques', linestyle='--')
# plt.plot([6,6], [mini, maxi], color='b', linewidth=1, linestyle='--')
# plt.plot([4.5,4.5], [mini, maxi], color='b', linewidth=1, linestyle='--')
plot_span([2.25,3.75,5.25],2,freq[start:end],'g',"Ultra-harmoniques",fft[start:end])
#plt.plot([2.25,2.25], [mini, maxi], color='g', linewidth=1, label = "Ultra-harmoniques", linestyle='--')
#plt.plot([3.75,3.75], [mini, maxi], color='g', linewidth=1, linestyle='--')
#plt.plot([5.25,5.25], [mini, maxi], color='g', linewidth=1, linestyle='--')
plot_span([1.5 + 1.5/4 + i*1.5/2 for i in range(5)],8,freq[start:end],'r',"Bruit large bande",fft[start:end])
# plt.plot([f0 * 1.25, f0 * 1.25], [mini, maxi], color='r', linewidth=1, label = "Bruit large bande", linestyle='--')
# plt.plot([f0 * 1.75, f0 * 1.75], [mini, maxi], color='r', linewidth=1, linestyle='--')
# plt.plot([f0 * 2.25, f0 * 2.25], [mini, maxi], color='r', linewidth=1, linestyle='--')
# plt.plot([f0 * 2.75, f0 * 2.75], [mini, maxi], color='r', linewidth=1, linestyle='--')
# plt.plot([f0 * 3.25, f0 * 3.25], [mini, maxi], color='r', linewidth=1, linestyle='--')
# plt.plot([f0 * 3.75, f0 * 3.75], [mini, maxi], color='r', linewidth=1, linestyle='--')
plt.scatter(0.8, 1, color='k', s=100, marker='o', label="Valeur moyenne de l'indice")
plt.legend(fontsize=16,framealpha = 1, loc = 'upper right')
plt.yscale('log')
# je veux construire les ticks_X avec un enchevetrement en f0 et 1.5*f0
#je veux que les f0 soient ecrit en italique avec 0 en indice du f
foo = r'$f_0$'
add = ' = 1.5 MHz'
ticks_X = []
for i in range(3):
    ticks_X.append(str(i+1)+foo)
    ticks_X.append(str(i+1.5)+foo)
ticks_X[0] += add
plt.xticks([1.5 + 1.5/2*i for i in range(6)], ticks_X,fontsize=16)
plt.title("Différentes composantes fréquentielles calculées à partir du spectre", fontsize=20, fontweight='bold')
plt.xlabel('Fréquence', fontsize=18)
plt.ylabel('Amplitude', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True)
plt.xlim(1, 5)
plt.ylim(mini, maxi)
plt.savefig(traj_plot+'spectre_fft.png', bbox_inches='tight')