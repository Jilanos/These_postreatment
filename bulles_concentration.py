# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:01:38 2022

@author: PM263553
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors, KNeighborsRegressor
from sklearn.metrics import r2_score
from classes import *
from scipy import signal
import sys 
import gc

gc.collect(generation=2)

pc_fixe = False
if pc_fixe:
    tra = 'D:\\code_UH_long\\GENE_MOD\\iter_13\\'
    traj='C:\\Users\\MIDAS\\Desktop\\resultats_concentration\\iter_13\\'
    doss=["PULSE_0_75","PULSE_240_75","PULSE_80_75","PULSE_27_75"]
    doss=["bubbles_0_75","bubbles_240_75","bubbles_80_75","bubbles_27_75"] 
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

#VITRO
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
k=5
# on cree la carte de cavitation
# les dimensions de la carte sont : nbr d'indice, nbr d'exp, nbr de pulses, nbr de fenetres, nbr de composantes
# on supprime les n_supr dernieres fenetres du pulse car moins interessantes
# cela permet de les regrouper par k
carte = test_m_exp.extract_composante(nbit,fit,legend = legend)
shape = np.shape(carte)
n_supr = shape[2]%k
carte_rognee = carte[:,:,:-n_supr,:,:]  #on supprime les fenetres en trop
shape = np.shape(carte_rognee)
#on regroupe les fenetres par k
#(exp, pulses, fenetres, indice, composantes)
carte_reshape = np.reshape(carte_rognee,(shape[0],shape[1],shape[2]//k,k,shape[3],shape[4]))
#(exp, pulses, fenetres//k, k, indice, composantes)
#on donnera à l'algo de labelisation toutes les fenêtres, chaque point sera composé d'un lot de k fenetres et de toutes les composantes 
shape = np.shape(carte_reshape)
n_exp = shape[0]
n_pulse = shape[1]
n_fen = shape[2]
last_dim = [shape[2], shape[3],shape[4],shape[5]]
print("shape de carte_reshape : ",shape)
# il faut donc réaplatire les dimensions 3,4,5
X_plat = np.reshape(carte_reshape,(shape[0],shape[1],shape[2],shape[3]*shape[4]*shape[5]))
# (exp, pulses, fenetres//k, k*indice*composantes)
# on veut maintenant rajouter la pression dans la dernière dimension
x_pre = np.zeros((n_exp,n_pulse,n_fen,k*shape[4]*shape[5]+1))
x_pre[:,:,:,:-1] = X_plat
# on va maintenant rajouter la pression dans la dernière dimension
for i_exp in range(n_exp):
    for i_pulse in range(n_pulse):
        for i_fen in range(n_fen):
            x_pre[i_exp,i_pulse,i_fen,-1] = i_pulse//rep *fit[0] + fit[1]
# on veut écrire une fct qui va nous donner les valeure suivantes : i = 0 : 0; i = 1 : 1, i = 2 : 3, i = 3 : 9, i = 4 : 27, ...
def puissance(i):
    if i==0:
        return 0
    else:
        return 3**(i-1) 

#on va créer l'array de sortie y :
y_pre = np.zeros((n_exp,n_pulse,n_fen))
for i_exp in range(n_exp):
    print("exp : ",i_exp," concentration : ",puissance(i_exp))
    for i_pulse in range(n_pulse):
        for i_fen in range(n_fen):
            y_pre[i_exp,i_pulse,i_fen] = puissance(i_exp)

# on va fusionner les 3 première dimensions
X_final = np.reshape(x_pre,(n_exp*n_pulse*n_fen,k*shape[4]*shape[5]+1))
Y_final = np.reshape(y_pre,(n_exp*n_pulse*n_fen))
shape_X = np.shape(X_final)
shape_Y = np.shape(Y_final)


# on va maintenant regarder la carte de concentration créé par la régression linéaire
def carte_concentration(data,chemin, titre,legend):
    fig=plt.figure(figsize=(20,11))
    fig, axes = plt.subplots(nrows=1, ncols=n_exp, figsize=(20,11))
    fig.suptitle(titre,fontsize=18, fontweight = 'bold',y=0.94)
    for i,ax in enumerate(axes.flat):
        im = ax.imshow(data[i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,11.,((nbit[1])*fit[0]+fit[1]),((nbit[0])*fit[0]+fit[1])])  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
        ax.title.set_text(legend[i].upper())
        ax.title.set_fontweight('bold')
        ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
        if i==0:
            ax.set_ylabel("Pression (kPa)")
            fig.colorbar(im,ax=axes.ravel().tolist())
            #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
    plt.savefig(chemin+titre+".png",bbox_inches='tight')
    plt.close('all')


X_train, X_test, y_train, y_test = train_test_split(X_final, Y_final, train_size=0.8, test_size=0.2, random_state=42)
print("\nModel : LinearRegression")
lr = LinearRegression()
lr.fit(X_train, y_train)
y_pred = lr.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
r2s = r2_score(y_test, y_pred)
lrs = lr.score(X_test, y_test)
print("MSE : ",mse)
print("R2 : ",r2s)
print("LR score : ",lrs)
c_pred = lr.predict(X_final)
carte_concentration(np.reshape(c_pred,(n_exp,n_pulse,n_fen)),traj,"LinearRegression",legend)

for k in [3,6,11,20]:
    print("\nModel : KNeighborsRegressor with k = ",k," neighbors")
    knn = KNeighborsRegressor(n_neighbors=k)
    knn.fit(X_train, y_train)
    y_pred = knn.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2s = r2_score(y_test, y_pred)
    knns = knn.score(X_test, y_test)
    print("MSE : ",mse)
    print("R2 : ",r2s)
    print("KNN score : ",knns)
    c_pred = knn.predict(X_final)
    carte_concentration(np.reshape(c_pred,(n_exp,n_pulse,n_fen)),traj,"KNeighborsRegressor_{}_neighbors".format(k),legend)  



sys.exit()




def carte_cluster(data,chemin, titre,legend):
    fig=plt.figure(figsize=(20,11))
    fig, axes = plt.subplots(nrows=1, ncols=n_exp, figsize=(20,11))
    fig.suptitle(titre,fontsize=18, fontweight = 'bold',y=0.94)
    for i,ax in enumerate(axes.flat):
        im = ax.imshow(data[i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,11.,((nbit[1])*fit[0]+fit[1]),((nbit[0])*fit[0]+fit[1])])  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
        ax.title.set_text(legend[i].upper())
        ax.title.set_fontweight('bold')
        ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
        if i==0:
            ax.set_ylabel("Pression (kPa)")
            fig.colorbar(im,ax=axes.ravel().tolist())
            #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
    plt.savefig(chemin+titre+".png",bbox_inches='tight')
    plt.close('all')


#on veut maintenant appliquer un algo de labelisation sur les points de la carte
#on va utiliser kmeans, dbscan ou gmm
for n_clusters in [2,3,4,5,6,7,8,9,10]:
    kmeans = KMeans(n_clusters=n_clusters, n_init=10).fit(X)
    y_kmeans = kmeans.predict(X)

    gmm = sklearn.mixture.GaussianMixture(n_components=n_clusters).fit(X)
    y_gmm = gmm.predict(X)

    # on va maintenant regarder le pourcentage de points de chaque cluster 

    pourcentage_kmeans = [np.count_nonzero(y_kmeans == i) for i in range(np.max(y_kmeans)+1)]
    pourcentage_gmm = [np.count_nonzero(y_gmm == i) for i in range(np.max(y_gmm)+1)]

    # plt.figure()
    # for i_meth, pour_res in enumerate([pourcentage_kmeans,pourcentage_dbscan,pourcentage_gmm]):
    #     #on veut afficher visuellement le résultat
    #     plt.subplot(3,1,i_meth+1)
    #     plt.bar(range(len(pour_res)),pour_res)
    #     plt.title("Methode "+str(i_meth))
    #     plt.xlabel("Cluster")
    #     plt.ylabel("Pourcentage de points")
    # plt.suptitle("Nombre de clusters : "+str(n_clusters))

    # on veut maintenant regarder la répartition des clusters en fonction de la pression et de la feêtre
    #pour cela on va créer un carte 2D (une pour chaque xp) avec en abscisse les fenêtres et en ordonnée les pressions
    #les différents clussters seront différenciés par des couleurs

    #on va d'abord reshape no Y :
    shape = np.shape(y_kmeans)
    y_kmeans_reshape = np.reshape(y_kmeans,(n_exp,n_pulse,n_fen))
    y_gmm_reshape = np.reshape(y_gmm,(n_exp,n_pulse,n_fen))

    #on va maintenant créer les cartes au travers d'une fonction
    carte_cluster(y_kmeans_reshape,traj,"kmeans_{}_clusters".format(n_clusters),legend)
    carte_cluster(y_gmm_reshape,traj,"gmm_{}_clusters".format(n_clusters),legend)

    # on va controler en refaisant apparaitre la carte de la composante inertielle
    control = True
    if control and n_clusters==2:
        control_reshape = np.reshape(X,(n_exp,n_pulse,n_fen,last_dim[1],last_dim[2], last_dim[3]))
        #(exp, pulses, fenetres//k, k, indice, composantes)
        control_inertielle = np.mean(control_reshape[:,:,:,:,2,:],axis=(3,4))  
        carte_cluster(control_inertielle,traj,"control_inertielle",legend)





