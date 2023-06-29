# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:01:38 2022

@author: PM263553
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sklearn as skl
from sklearn import linear_model
from sklearn import ensemble
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from classes import *
from scipy import signal
import sys 
import gc
import random
import time

gc.collect(generation=2)


def choix_model(mod):
    if mod =="LinearSVR" :
        model = skl.svm.LinearSVR(max_iter=50000)
        espacer="\t"
    elif mod =="Lasso"  :
        model = skl.linear_model.Lasso(max_iter=50000)
        espacer="\t\t"
    elif mod =="ElasticNet"  :
        model = skl.linear_model.ElasticNet(max_iter=50000)
        espacer="\t"
    elif mod =="ElasticNetCV"  :
        model = skl.linear_model.ElasticNetCV(max_iter=50000)
        espacer="\t"
    elif mod =="Ridge"  :
        model = skl.linear_model.Ridge()
        espacer="\t\t"
    elif mod =="SGDRegressor"  :
        model = skl.linear_model.SGDRegressor(max_iter=50000)
        espacer="\t"
    elif mod =="SVR"  :
        model = skl.svm.SVR(cache_size=2000)
        espacer="\t\t"
    elif mod =="RandomForestRegressor"  :
        model = skl.ensemble.RandomForestRegressor()
        espacer=""
    elif mod =="KNeighborsRegressor"  :
        model = skl.neighbors.KNeighborsRegressor(n_neighbors = 8)
        espacer=""
    elif mod =="LinearRegression"  :
        model = skl.linear_model.LinearRegression()
        espacer=""
    return model, espacer

def carte_concentration(data,chemin, titre,legend):
    fig=plt.figure(figsize=(20,11))
    color = ['black', 'lightcoral', 'red','darkred']
    for i in range(np.shape(data)[0]):
        plt.plot(data[i,:], label = legend[i], color = color[i])
    plt.legend(fontsize = 15)
    plt.ylabel("Concentration prédite",fontsize = 15)
    plt.xlabel("Pulses",fontsize = 15)
    plt.title("Variation de la concentration pour chaque concentration selon la pression", fontsize = 17)
    plt.savefig(chemin+titre+".png",bbox_inches='tight')
    plt.close('all')

def concent(i):
    if i==0:
        return 0
    else:
        return 3**(i-1) 
fit = np.array([7.92060316*2, 2.42161125])
sys.exit()

pc_fixe = True
if pc_fixe:
    tra = 'D:\\code_UH_long\\GENE_MOD\\iter_21\\'
    traj='C:\\Users\\MIDAS\\Desktop\\resultats_concentration\\iter_21\\'
    doss=["bubbles_0_75","bubbles_240_75","bubbles_80_75","bubbles_27_75"] 
    doss=["PULSE_0_75","PULSE_240_75","PULSE_80_75","PULSE_27_75"]
    legend=["Eau","Dilution 1:240","Dilution 1:80","Dilution 1:27"]
    nexp=4
    rep = 20
else:
    tra = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\'
    traj = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\resultats_concentration\\'
    doss=["PULSE_0_40","PULSE_666_40"]
    legend=["Eau","Dilution 1:666"]
    nexp=2
    rep = 6
path(traj)


start,end = 0,-1

test_m_exp=experiment_mult(25000000,1500000,start=start,end=end,size_decal = 1024)

fit = np.array([7.92060316*2, 2.42161125])
pression_max = 1100
pression_min = 1
bitmax=int(np.round(((pression_max-fit[1])/fit[0])))
bitmin=int(np.round(((pression_min-fit[1])/fit[0])))
bitpress = int(np.round(((400-fit[1])/fit[0])))
nbit=[1,bitmax]


order = True
for j in range(nexp):
    test_m_exp.creat_expe()

for i in range(0,nexp):
    print("\nloading data : "+doss[i])
    dossier=tra+doss[i]
    traj_comp = tra+doss[i]
    path(traj_comp)
    data = np.load(dossier+'\\data_o.npy') #
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:rep*(bitmax-1)], i, spacer =100e3)
    print('done')
    del data
    gc.collect(generation=2)


#%%

# on cree la carte de cavitation
# les dimensions de la carte sont : nbr d'indice, nbr d'exp, nbr de pulses, nbr de fenetres, nbr de composantes
# on supprime les n_supr dernieres fenetres du pulse car moins interessantes
# cela permet de les regrouper par k
carte = test_m_exp.extract_composante(nbit,fit,legend = legend)
shape = np.shape(carte)
#on regroupe les fenetres par k
#(exp, pulses, fenetres, indice, composantes)
#carte_reshape = np.reshape(carte_rognee,(shape[0],shape[1],shape[2]//k,k,shape[3],shape[4]))
#(exp, pulses, fenetres//k, k, indice, composantes)
#on donnera à l'algo de labelisation toutes les fenêtres, chaque point sera composé d'un lot de k fenetres et de toutes les composantes 
n_exp = shape[0]
n_pulse = shape[1]
n_fen = shape[2]
# il faut donc réaplatire les dimensions 3,4,5
X_plat = np.reshape(carte,(shape[0],shape[1],shape[2]*shape[3]*shape[4]))
# (exp, pulses, fenetres//k, k*indice*composantes)
# on veut maintenant rajouter la pression dans la dernière dimension
x_pre = np.zeros((n_exp,n_pulse,n_fen*shape[3]*shape[4]+1))
x_pre[:,:,:-1] = X_plat
# on va maintenant rajouter la pression dans la dernière dimension
for i_exp in range(n_exp):
    for i_pulse in range(n_pulse):
        x_pre[i_exp,i_pulse,-1] = i_pulse//rep *fit[0] + fit[1]
# on veut écrire une fct qui va nous donner les valeure suivantes : i = 0 : 0; i = 1 : 1, i = 2 : 3, i = 3 : 9, i = 4 : 27, ...


#on va créer l'array de sortie y :
y_pre = np.zeros((n_exp,n_pulse))
for i_exp in range(n_exp):
    print("exp : ",i_exp," concentration : ",concent(i_exp))
    for i_pulse in range(n_pulse):
        y_pre[i_exp,i_pulse] = concent(i_exp)
#%%
# on va fusionner les 3 première dimensions
X_final = np.reshape(x_pre,(n_exp*n_pulse,n_fen*shape[3]*shape[4]+1))
Y_final = np.reshape(y_pre,(n_exp*n_pulse))
shape_X = np.shape(X_final)
shape_Y = np.shape(Y_final)


# on va maintenant regarder la carte de concentration créé par la régression linéaire


prop_train = 0.05
X_train, X_test, y_train, y_test = train_test_split(X_final, Y_final, train_size=prop_train, test_size=1-prop_train, random_state=42)
models = ["Lasso","ElasticNet","ElasticNetCV","Ridge","SVR","LinearSVR","RandomForestRegressor","KNeighborsRegressor","LinearRegression"]
for mod in models:
    model, espacer = choix_model(mod)
    print(mod)
    t1=time.time()
    model.fit(X_train, y_train)
    t_train = int(time.time()-t1)
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2s = r2_score(y_test, y_pred)
    lr = model.score(X_test, y_test)
    print("MSE : ",mse)
    print("R2 : ",r2s)
    print(mod+" score : ",lr)
    t1=time.time()
    c_pred = model.predict(X_final)
    carte_concentration(np.reshape(c_pred,(n_exp,n_pulse)),traj,mod,legend)
    t_crea = int(time.time()-t1)
    file1 = open(traj+"log_record.txt", "a")
    file1.write(mod+espacer+" \ttrain size = {}% \tacc = {} \ttrain time = {}s \tplot time = {}s\n".format(int(np.round(prop_train*100,decimals=1)), lr,t_train,t_crea))
    file1.close()  

#%%
#Version pulse rampe
pc_fixe = True
if pc_fixe:
    tra = 'D:\\code_UH_long\\GENE_MOD\\iter_21\\'
    
    traj='C:\\Users\\MIDAS\\Desktop\\resultats_concentration\\iter_21_ramp_bis\\'
    doss=["bubbles_0_75","bubbles_240_75","bubbles_80_75","bubbles_27_75"] 
    doss=["RAMP_0_75","RAMP_240_75","RAMP_80_75","RAMP_27_75"]
    legend=["Eau","Dilution 1:240","Dilution 1:80","Dilution 1:27"]
    nexp=4
    rep = 20
else:
    tra = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\'
    traj = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\in_vitro\\iter_19\\resultats_concentration\\'
    doss=["RAMP_0_40","RAMP_666_40"]
    legend=["Eau","Dilution 1:666"]
    nexp=2
    rep = 6
path(traj)
start,end = 0,-1
test_file1 = "D:\\code_UH_long\\GENE_MOD\\VIVO_1\\VIVO_Mouse564_221122_65\\"
test_file2 = "D:\\code_UH_long\\GENE_MOD\\VIVO_1\\VIVO_Mouse587_221122_45\\"
test_file3 = "D:\\code_UH_long\\GENE_MOD\\VIVO_1\\VIVO_rat409_291122_86\\"
test_file4 = "D:\\code_UH_long\\GENE_MOD\\VIVO_1\\VIVO_rat_445_030123_75\\"
legend_mice = ["mouse 564","mouse 587","rat 409","rat 445"]

path_mice = [test_file1,test_file2,test_file3,test_file4]
fit = np.array([7.92060316*2, 2.42161125])
for i, elt in enumerate(path_mice):
    print("loading pulses exp souris test")
    mice_exp=experiment_mult(25000000,1500000,start=start,end=end,size_decal = 1024)
    mice_exp.creat_expe()
    data = np.load(elt+'data_treatment.npy') #
    print("adding pulses multi exp souris test ")
    mice_exp.add_pulses(data[:], 0, spacer =100e3)
    carte_mice = mice_exp.extract_composante(False,fit,legend = legend)
    print('done')
    del data
    gc.collect(generation=2)
    shape = np.shape(carte_mice)
    n_exp = shape[0]
    n_pulse = shape[1]
    n_pulse_mice = n_pulse
    n_fen = shape[2]
    mice_plat = np.reshape(carte_mice,(shape[0],shape[1],shape[2]*shape[3]*shape[4]))

    x_mice_pre = mice_plat
    X_mice_final = np.reshape(x_mice_pre,(n_exp*n_pulse,n_fen*shape[3]*shape[4]))
    shape_mice = np.shape(X_mice_final)
    if i==0:
        X_mice_data = np.zeros((len(path_mice),shape_mice[0],shape_mice[1]))
    X_mice_data[i,:,:] = X_mice_final
    print("map of concent for exp souris created ")
    del mice_exp

test_m_exp=experiment_mult(25000000,1500000,start=start,end=end,size_decal = 1024)



order = True
for j in range(nexp):
    test_m_exp.creat_expe()

for i in range(0,nexp):
    print("\nloading data : "+doss[i])
    dossier=tra+doss[i]
    traj_comp = tra+doss[i]
    path(traj_comp)
    data = np.load(dossier+'\\data.npy') #
    print("adding pulses multi exp")
    test_m_exp.add_pulses(data[:], i, spacer =100e3)
    print('done')
    del data
    gc.collect(generation=2)




# on cree la carte de cavitation
# les dimensions de la carte sont : nbr d'indice, nbr d'exp, nbr de pulses, nbr de fenetres, nbr de composantes
# on supprime les n_supr dernieres fenetres du pulse car moins interessantes
# cela permet de les regrouper par k
carte = test_m_exp.extract_composante(False,fit,legend = legend)
shape = np.shape(carte)
#on regroupe les fenetres par k
#(exp, pulses, fenetres, indice, composantes)
#carte_reshape = np.reshape(carte_rognee,(shape[0],shape[1],shape[2]//k,k,shape[3],shape[4]))
#(exp, pulses, fenetres//k, k, indice, composantes)
#on donnera à l'algo de labelisation toutes les fenêtres, chaque point sera composé d'un lot de k fenetres et de toutes les composantes 
n_exp = shape[0]
n_pulse = shape[1]
n_fen = shape[2]
# il faut donc réaplatire les dimensions 3,4,5
X_plat = np.reshape(carte,(shape[0],shape[1],shape[2]*shape[3]*shape[4]))

x_pre = X_plat
#on va créer l'array de sortie y :
y_pre = np.zeros((n_exp,n_pulse))
for i_exp in range(n_exp):
    print("exp : ",i_exp," concentration : ",concent(i_exp))
    for i_pulse in range(n_pulse):
        y_pre[i_exp,i_pulse] = concent(i_exp)

# on va fusionner les 3 première dimensions
X_final = np.reshape(x_pre,(n_exp*n_pulse,n_fen*shape[3]*shape[4]))
Y_final = np.reshape(y_pre,(n_exp*n_pulse))
shape_X = np.shape(X_final)
shape_Y = np.shape(Y_final)


# on va maintenant regarder la carte de concentration créé par la régression linéaire


prop_train = 0.1
X_train, X_test, y_train, y_test = train_test_split(X_final, Y_final, train_size=prop_train, test_size=1-prop_train, random_state=42)
models = ["Lasso","ElasticNet","ElasticNetCV","Ridge","SVR","LinearSVR","RandomForestRegressor","KNeighborsRegressor","LinearRegression"]
for mod in models:
    model, espacer = choix_model(mod)
    print(mod)
    t1=time.time()
    model.fit(X_train, y_train)
    t_train = int(time.time()-t1)
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2s = r2_score(y_test, y_pred)
    lr = model.score(X_test, y_test)
    print("MSE : ",mse)
    print("R2 : ",r2s)
    print(mod+" score : ",lr)
    t1=time.time()
    c_pred = model.predict(X_final)
    c_mice_pred = model.predict(np.reshape(X_mice_data,(shape_mice[0]*len(path_mice),shape_mice[1])))
    carte_concentration(np.reshape(c_mice_pred,(len(path_mice),n_pulse_mice)),traj+"mice_",mod,legend_mice)
    carte_concentration(np.reshape(c_pred,(n_exp,n_pulse)),traj,mod,legend)
    t_crea = int(time.time()-t1)
    file1 = open(traj+"log_record.txt", "a")
    file1.write(mod+espacer+" \ttrain size = {}% \tacc = {} \ttrain time = {}s \tplot time = {}s\n".format(int(np.round(prop_train*100,decimals=1)), lr,t_train,t_crea))
    file1.close()  

