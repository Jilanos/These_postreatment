# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 14:25:52 2022

@author: PM263553
"""
import numpy as np
import matplotlib.pyplot as plt
import os

from scipy import signal

plt.close("all")


def sortir_bruit(signal, seuil_min, debut_detec):
    n = len(signal)
    res = np.zeros((n))
    for j in range(debut_detec, n):
        moy = np.mean(signal[:j])
        std = np.std(signal[:j])
        depassement = (signal[j]-moy)/std
        if depassement >= seuil_min:
            res[j] = depassement
    maxi, maxi_ind = max_pos(res)
    return maxi, maxi_ind

     
def max_pos(array):
    n=len(array)
    maxi,Imaxi = array[0],0
    for i in range(n):
        elt = array[i]
        if elt>maxi:
            maxi,Imaxi=elt,i
    return maxi,Imaxi 



def gradation(n):
    if n==0 :
        return 0
    else :
        return (n-1)*0.20

def plott(X,Y,chemin,color='black',titre = "titre",Xlabel="X axis",Ylabel="Y axis"):
    plt.figure(figsize = (20,11))
    plt.plot(X,Y,c=color)
    plt.title(titre,color=color,fontsize=30, fontweight = 'bold')
    plt.xlabel(Xlabel,fontsize=20)
    plt.ylabel(Ylabel, color=color,fontsize=20)
    Ymax = 1.01*np.amax(Y)
    Ymin = 0.99*np.amin(Y)
    plt.ylim([Ymin, Ymax])    
    plt.grid(True)
    plt.tight_layout()
    plt.yscale('linear')
    plt.savefig(chemin,bbox_inches='tight')
    plt.close("all")
    return 1

def algo(amp,val,rep,seuil,pas,p,val_dic):
    ptitpas = max([int(pas/2),1])
    if amp.count(amp[-1])<rep:
        return amp[-1], False
    elif max(amp)<seuil:
        return max(amp)+pas, False
    else:
        courbe_Y = []
        courbe_X = []
        for i in val_dic:
            if len(val_dic[i])>0:
                courbe_Y.append(moyenne_histo(val_dic[i],p))
                courbe_X.append(i)
        courbe_Y = list(filtr(courbe_X,courbe_Y))
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
                return amp_1+pas, False
            else :
                #print("moins")
                return amp_1-pas, False

def gauss(x,b,c):
    e = -(x-b)**2/2./c
    return np.exp(e)

def filtr(Xarr, Yarr, sigma = 4):
    n= len(Xarr)
    YarrF = np.copy(Yarr)
    for j in range(n):
        deb, fin = -2, 3 
        if j == 0 : deb = 0
        if j == 1 : deb = -1
        if j == n-1 : fin = 1
        if j == n-2 : fin = 2
        poids = np.array([gauss(Xarr[j+i],Xarr[j],sigma) for i in range(deb,fin)])
        somme = np.sum(poids)
        YarrF[j]=np.sum(np.multiply(Yarr[j+deb:j+fin],poids))/somme
    return YarrF

def valeurs(data, press_max):
    VAL = []
    m = np.mean(data, axis = 0)
    m = m-np.mean(m)
    size_decal = 1024
    size_window = 2048
    Nb_spectres=len(m) //size_decal-1
    
    for i_spectro in range(Nb_spectres):
        sig=m[i_spectro*size_decal:i_spectro*size_decal+size_window]
        #Sig_fen=np.hanning(size_window)*sig
        Sig_fen=signal.windows.flattop(size_window)*sig
        VAL.append(abs(np.min(Sig_fen)))
    max_val = np.mean(np.array(VAL[-4:]))
    conv = press_max / max_val
    VAL = conv * np.array(VAL)
    plt.plot(VAL)
    return VAL

def valeurs_rampe(data, rep, p_deb, size_window = 2048):
    m = np.mean(data[-rep::], axis = 0)
    m = m-np.mean(m)
    maxi = np.max(m)
    pos_max = 0
    for i in range(len(m)):
        if m[i]>= maxi :
            pos_max = i
    d = pos_max - p_deb
    n_fen = int(d/size_window)
    return n_fen, pos_max
    
    
 

def moyenne_histo(liste, percent = 0.3):
    if len(liste)==1:
        return liste[0]
    else :
        return liste[-1] * percent + (1 - percent) * moyenne_histo(liste[:-1], percent)

def path(path):
    if not os.path.isdir(path): # check if folder exists, otherwise create it
        os.mkdir(path)

def nommage(num,den):
    nom=""
    for j in num:
        nom += "H{}".format(j)
    nom += "_"
    for j in den:
        nom += "H{}".format(j)
    return nom
     
def moyenning(carte,n):
    dim=np.shape(carte)
    carte_bis = np.zeros((dim[0],dim[1],dim[2]-n+1))
    for i in range(dim[2]-n+1):
        carte_bis[:,:,i]=np.mean(carte[:,:,i:i+n],axis=2)
    return carte_bis  

def anteriorite(carte,concent,ant):
    dim=np.shape(carte)
    carte_bis = []
    concent_list = []
    for i in range(dim[0]):
        for j in range(dim[1]):#self.exp[i].pulses[j].indice_harm_w[1]+
            for k in range(dim[2]-ant+1):
                carte_bis.append(np.concatenate([carte[i,j,k-ij] for ij in range(ant-1,-1,-1)]))
                concent_list.append(concent[i,j,k])
    return carte_bis, concent_list


def reorder(val_ini,amp_ini):
    val , amp =  np.ndarray.tolist(np.copy(val_ini)),  np.ndarray.tolist(np.copy(amp_ini))
    new_val =[]
    new_amp =[]
    while len(amp)>0:
        ind = np.argmin(amp)
        new_val.append(val[ind])
        new_amp.append(amp[ind])
        del(val[ind])
        del(amp[ind])
    return new_val,new_amp



def peak_detect(array):
    n=len(array)
    mini,Imini,maxi,Imaxi = array[0],0,array[0],0
    for i in range(n):
        elt = array[i]
        if elt>maxi:
            maxi,Imaxi=elt,i
        if elt<mini:
            mini,Imini=elt,i
    moy= np.mean(array)
    if abs(mini-moy)>abs(maxi-moy):
        return mini,Imini
    else:
        return maxi,Imaxi
    
def max_detect(array):
    n=len(array)
    maxi,Imaxi = array[0],0
    for i in range(n):
        elt = array[i]
        if elt>maxi:
            maxi,Imaxi=elt,i
    moy= np.mean(array)
    std = np.std(array)
    if maxi> moy+1.5*std:
        return maxi,Imaxi
    else :
        return moy, n//2
    
def peak_detect_pos(array):
    n=len(array)
    maxi,Imaxi = array[0],0
    for i in range(n):
        elt = array[i]
        if elt>maxi:
            maxi,Imaxi=elt,i
    moy= np.mean(array)
    if np.std(array)>abs(maxi-moy):
        return moy,int(len(array)/2)
    else:
        return maxi,Imaxi

def fft3(Data,size_win):
    SIG_FFT=[]
    Nb_spectres=len(Data)//size_win
    for i_spectro in range(Nb_spectres-1):
        sig=Data[i_spectro*size_win:(i_spectro+1)*size_win]
        Sig_fen=np.hanning(2048)*sig
        SIG_FFT.append(abs(np.fft.fft(Sig_fen)))
    return np.array(SIG_FFT)

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

def p_to_b(fit,p):
    return int(np.round(((p-fit[1])/fit[0])))


def b_to_p(fit,b):
    return int(np.round(fit[1]+fit[0]*b))


class temp_sample():
    def __init__(self, fe, f0, data,start:int,end:int, size_window : int = 2048, n_harm : int = 4,n_BB : int =5, delta_harm =50e3,delta_Uharm =50e3, spacer = 50e3,peak_detection : bool =True, size_decal : int = 2048):
        self.data=np.copy(np.array(data[start:end]))
        self.data=self.data-np.mean(self.data)
        self.freq_ech = float(fe)
        self.freq_emiss = float(f0)
        self.num_samples = len(self.data)
        self.num_BB = n_BB
        self.n_harm=n_harm
        self.start=start
        self.end=end
        self.freq_harm=np.array([self.freq_emiss*i for i in range(1,n_harm+2)])
        self.freq_Uharm=np.array([self.freq_emiss*(i-0.5) for i in range(1,n_harm+2)])
        self.freq_Uharm_BBref1=np.array([self.freq_emiss*(i-0.75) for i in range(1,n_harm+2)])
        self.freq_Uharm_BBref2=np.array([self.freq_emiss*(i-0.25) for i in range(1,n_harm+2)])
        self.freq_BB=np.array([self.freq_emiss*(i/2.+1.75) for i in range(self.num_BB)])
        
        self.peak_detection=peak_detection
        
        
        self.size_window=size_window
        self.size_decal=size_decal
        self.spectrogram=self.spectogramage()
        self.n_window = len(self.spectrogram)
        self.fftsimple=np.mean(self.spectrogram,axis=0)
        #self.fftsimple=self.spectoUnik()
        self.frequencies= np.arange(self.size_window//2) * self.freq_ech/self.size_window
        #self.frequencies= np.arange(self.num_samples//2) * self.freq_ech/self.num_samples
        self.delta_harm=delta_harm
        self.delta_Uharm=delta_Uharm
        self.delta_BB=0.25*self.freq_emiss-max(self.delta_harm,self.delta_Uharm)-spacer
        self.freq_harm_n= list(map(lambda x: find_nearest(self.frequencies, x), self.freq_harm))
        self.freq_fonda_n= find_nearest(self.frequencies, self.freq_emiss)
        self.freq_Uharm_n= list(map(lambda x: find_nearest(self.frequencies, x), self.freq_Uharm))
        self.freq_Uharm_BBref1_n= list(map(lambda x: find_nearest(self.frequencies, x), self.freq_Uharm_BBref1))
        self.freq_Uharm_BBref2_n= list(map(lambda x: find_nearest(self.frequencies, x), self.freq_Uharm_BBref2))
        self.freq_BB_n= list(map(lambda x: find_nearest(self.frequencies, x), self.freq_BB))
        self.delta_Uharm_n = find_nearest(self.frequencies, self.delta_Uharm)
        self.delta_harm_n = find_nearest(self.frequencies, self.delta_harm)
        self.delta_BB_n = find_nearest(self.frequencies, self.delta_BB)
        
        self.broadband_fft_copy=np.mean(self.spectrogram,axis=0) # [0:12]
        self.indice_harm=np.zeros((self.n_harm+1))
        self.indice_Uharm=np.zeros((self.n_harm+1))
        self.indice_BB=0
        self.indice_BB_sliced=np.zeros((self.num_BB))

        
        self.indice_calc()
        self.BB_calc()
        
        self.indice_harm_w=np.zeros((self.n_harm+1,self.n_window))
        self.indice_Uharm_w=np.zeros((self.n_harm+1,self.n_window))
        self.indice_Uharm_norm_div_w=np.zeros((self.n_harm+1,self.n_window))
        self.indice_Uharm_norm_div=np.zeros((self.n_harm+1))
        self.indice_BB_w=np.zeros((self.n_window))
        self.indice_BB_sliced_w=np.zeros((self.num_BB,self.n_window))
        self.fondamental_w = np.zeros((self.n_window))
        self.indice_calc_windowed()
        self.BB_calc_windowed()
        del(self.data)
        
        
    def calcul_indice(self,arr,freq_0,freq_1):
        if self.peak_detection:
            val=np.max(arr[freq_0:freq_1+1])
        else :
            val=np.mean(arr[freq_0:freq_1+1])
        return val
    
    def indice_calc(self):
        for i in range(self.n_harm+1):
            self.indice_harm[i]=self.calcul_indice(self.fftsimple,self.freq_harm_n[i]-self.delta_harm_n,self.freq_harm_n[i]+self.delta_harm_n+1)
        for i in range(self.n_harm+1):
            self.indice_Uharm[i]=self.calcul_indice(self.fftsimple,self.freq_Uharm_n[i]-self.delta_Uharm_n,self.freq_Uharm_n[i]+self.delta_Uharm_n+1)
            
    def indice_calc_windowed(self):
        for n in range(self.n_window):
            for i in range(self.n_harm+1):
                self.indice_harm_w[i,n]=self.calcul_indice(self.spectrogram[n,:],self.freq_harm_n[i]-self.delta_harm_n,self.freq_harm_n[i]+self.delta_harm_n+1)
            for i in range(self.num_BB):
                self.peak_detection = False
                self.indice_BB_sliced_w[i,n]=self.calcul_indice(self.spectrogram[n,:],self.freq_BB_n[i]-self.delta_BB_n,self.freq_BB_n[i]+self.delta_BB_n+1)
                self.peak_detection = True
            for i in range(self.n_harm+1):
                self.indice_Uharm_w[i,n]=self.calcul_indice(self.spectrogram[n,:],self.freq_Uharm_n[i]-self.delta_Uharm_n,self.freq_Uharm_n[i]+self.delta_Uharm_n+1)
            self.fondamental_w[n]=self.calcul_indice(self.spectrogram[n,:],self.freq_fonda_n-self.delta_harm_n,self.freq_fonda_n+self.delta_harm_n+1)
            for i in range(self.n_harm+1):
                val=self.indice_Uharm_w[i,n]
                ref1=self.calcul_indice(self.spectrogram[n,:],self.freq_Uharm_BBref1_n[i]-self.delta_BB_n,self.freq_Uharm_BBref1_n[i]+self.delta_BB_n+1)
                ref2=self.calcul_indice(self.spectrogram[n,:],self.freq_Uharm_BBref2_n[i]-self.delta_BB_n,self.freq_Uharm_BBref2_n[i]+self.delta_BB_n+1)
                self.indice_Uharm_norm_div_w[i,n] = val*2 / (ref1 + ref2)
        self.indice_Uharm_norm_div = np.mean(self.indice_Uharm_norm_div_w, axis = 1)
            
    
    def BB_calc(self):
        self.broadband_value=np.array([])
        for i in range(self.num_BB):
            self.indice_BB_sliced[i]=np.mean(self.broadband_fft_copy[self.freq_BB_n[i]-self.delta_BB_n:self.freq_BB_n[i]+self.delta_BB_n])
            self.broadband_value=np.concatenate((self.broadband_value,self.broadband_fft_copy[self.freq_BB_n[i]-self.delta_BB_n:self.freq_BB_n[i]+self.delta_BB_n])) 
        self.indice_BB=np.mean(self.broadband_value[0:2])

    def BB_calc_windowed(self):
        for n in range(self.n_window):
            self.broadband_value=np.array([])
            for i in range(self.num_BB):
                self.broadband_value=np.concatenate((self.broadband_value,self.spectrogram[n,self.freq_BB_n[i]-self.delta_BB_n:self.freq_BB_n[i]+self.delta_BB_n])) 
            self.indice_BB_w[n]=np.mean(self.broadband_value)      

        
    def spectogramage(self):
        SIG_FFT=[]
        Nb_spectres=self.num_samples //self.size_decal-3
        for i_spectro in range(Nb_spectres):
            sig=self.data[i_spectro*self.size_decal:i_spectro*self.size_decal+self.size_window]
            #Sig_fen=np.hanning(self.size_window)*sig
            Sig_fen=signal.windows.flattop(self.size_window)*sig
            SIG_FFT.append(abs(np.fft.fft(Sig_fen))/self.size_window)
        return np.array(SIG_FFT)
    
    def spectoUnik(self):
        sig=self.data
        #Sig_fen=np.hanning(self.num_samples)*sig
        Sig_fen=signal.windows.flattop(self.num_samples)*sig
        SIG_FFT=abs(np.fft.fft(Sig_fen))/self.num_samples
        return np.array(SIG_FFT)
    
    def plot_window(self):
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(18,10))
        axes[0].plot(np.mean(self.indice_harm_w, axis =0))
        axes[1].plot(np.mean(self.indice_Uharm_w, axis =0))
        axes[2].plot(self.indice_BB_w)
        plt.figure(figsize=(18,10))
        print(np.shape(self.spectrogram))
        plt.imshow(10*np.log10(self.spectrogram[:,:512]))
        plt.show()
        
        
    def plot(self,chemin):
        plt.figure(figsize=(20,11))
        spec=self.fftsimple[:self.size_window//2]#size_window
        Y=spec
        Ymax = 1.05*np.amax(spec)
        Ymin = 0.99*np.amin(spec)
        
        plt.plot(self.frequencies/1e6, Y, 'k', label="FFT", linewidth=2.0)
        # Fundamental
        plt.plot((self.freq_emiss/1e6 , self.freq_emiss/1e6), (0, Ymax), 'r-')
        # Inertial bandwidth
        for n in range(self.num_BB):
            #plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='w',alpha=None,lw=0)
            plt.axvspan((self.freq_BB[n]-self.delta_BB)/1e6,(self.freq_BB[n]+self.delta_BB)/1e6, color='r',alpha=0.2,lw=0)
            plt.plot(((self.freq_BB[n]-self.delta_BB)/1e6,(self.freq_BB[n]+self.delta_BB)/1e6), (self.indice_BB_sliced[n], self.indice_BB_sliced[n]), 'r', linewidth=2.0)
            #plt.plot(((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6), (self.indice_harm[n], self.indice_harm[n]), 'b', linewidth=1.0)
        #plt.axvspan(self.frequencies[self.freq_Uharm_n[1]]/1e6,self.frequencies[self.freq_harm_n[2]]/1e6, color='r',alpha=0.2,lw=0)  # ultraharm_f[0]+harm_df
# =============================================================================
#         # Sub-harmonic
#         plt.plot(((subharm_f-harm_df)/1e6 , (subharm_f-harm_df)/1e6), (0, Ymax), 'y--', linewidth=2.0)
#         plt.plot(((subharm_f+harm_df)/1e6 , (subharm_f+harm_df)/1e6), (0, Ymax), 'y--', linewidth=2.0)
#         plt.axvspan((subharm_f-harm_df)/1e6,(subharm_f+harm_df)/1e6, color='y',alpha=0.2,lw=0)
# =============================================================================
        # Harmonics
        plt.plot(((self.freq_harm-self.delta_harm)/1e6 , (self.freq_harm-self.delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        plt.plot(((self.freq_harm+self.delta_harm)/1e6 , (self.freq_harm+self.delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        for n in range(self.n_harm+1):
            #plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='w',alpha=None,lw=0)
            plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='b',alpha=0.2,lw=0)
            #plt.scatter(self.freq_harm[n]/1e6, self.indice_harm[n], 'b', s=100, marker='o')
            plt.plot(((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6), (self.indice_harm[n], self.indice_harm[n]), 'b', linewidth=2.0)
        # Ultra-harmonics
        plt.plot(((self.freq_Uharm-self.delta_Uharm)/1e6 , (self.freq_Uharm-self.delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        plt.plot(((self.freq_Uharm+self.delta_Uharm)/1e6 , (self.freq_Uharm+self.delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        for n in range(self.n_harm+1):
            #plt.axvspan((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6, color='w',alpha=None,lw=0)
            plt.axvspan((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6, color='g',alpha=0.2,lw=0)
            #plt.scatter(self.freq_Uharm[n]/1e6, self.indice_Uharm[n], 'g', s=100, marker='o')
            plt.plot(((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6), (self.indice_Uharm[n], self.indice_Uharm[n]), 'g', linewidth=2.0)
        plt.grid(True, which='major')
        plt.ylabel('Magnitude',fontsize=20)
        plt.yscale('log')
        plt.xlabel('Frequency [MHz]',fontsize=20)
        plt.title('Spectre moyenné du pulse avec visualisation des différentes bandes fréquentielles',fontsize=20, fontweight = 'bold')
        plt.xlim([self.freq_emiss/1e6/5, (self.freq_harm[self.n_harm-1]+self.delta_harm)/1e6*1.01])#(harm_f[-2]-2*harm_df)/1e6])
        plt.ylim([Ymin, Ymax])
        #plt.xlim([self.freq_emiss/1e6/5, 9.5*1.01])#(harm_f[-2]-2*harm_df)/1e6])
        #plt.ylim([1e-3, Ymax])
        
        plt.savefig(chemin+'spec.png',bbox_inches='tight')
        plt.close("all")

    def plot_FR(self,chemin):
        plt.figure(figsize=(18,8))
        spec=self.fftsimple[:self.size_window//2]#size_window
        Y=spec
        Ymax =3.5*np.amax(spec)
        Ymin = 0.995*np.amin(spec)
        plt.scatter(15,0,c='k',s=100,marker='o',label='Valeurs des indices')
        plt.plot(self.frequencies/1e6, Y, 'k', label="FFT", linewidth=2.0)
        # Fundamental
        plt.plot((self.freq_emiss/1e6 , self.freq_emiss/1e6), (0, Ymax), 'r-', label='Fondamental')
        # Inertial bandwidth
        for n in range(self.num_BB):
            #plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='w',alpha=None,lw=0)
            if n==0:
                plt.axvspan((self.freq_BB[n]-self.delta_BB)/1e6,(self.freq_BB[n]+self.delta_BB)/1e6, color='r',alpha=0.2,lw=0, label='Bruit large bande')
            else:
                plt.axvspan((self.freq_BB[n]-self.delta_BB)/1e6,(self.freq_BB[n]+self.delta_BB)/1e6, color='r',alpha=0.2,lw=0)
            #plt.plot(((self.freq_BB[n]-self.delta_BB)/1e6,(self.freq_BB[n]+self.delta_BB)/1e6), (self.indice_BB_sliced[n], self.indice_BB_sliced[n]), 'r', linewidth=2.0)
        plt.scatter(self.freq_BB/1e6, self.indice_BB_sliced, c='r', s=100, marker='o')
            #plt.plot(((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6), (self.indice_harm[n], self.indice_harm[n]), 'b', linewidth=1.0)
        #plt.axvspan(self.frequencies[self.freq_Uharm_n[1]]/1e6,self.frequencies[self.freq_harm_n[2]]/1e6, color='r',alpha=0.2,lw=0)  # ultraharm_f[0]+harm_df
# =============================================================================
#         # Sub-harmonic
#         plt.plot(((subharm_f-harm_df)/1e6 , (subharm_f-harm_df)/1e6), (0, Ymax), 'y--', linewidth=2.0)
#         plt.plot(((subharm_f+harm_df)/1e6 , (subharm_f+harm_df)/1e6), (0, Ymax), 'y--', linewidth=2.0)
#         plt.axvspan((subharm_f-harm_df)/1e6,(subharm_f+harm_df)/1e6, color='y',alpha=0.2,lw=0)
# =============================================================================
        # Harmonics
        plt.plot(((self.freq_harm-self.delta_harm)/1e6 , (self.freq_harm-self.delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        plt.plot(((self.freq_harm+self.delta_harm)/1e6 , (self.freq_harm+self.delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
        for n in range(self.n_harm+1):
            #plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='w',alpha=None,lw=0)
            if n==0:
                plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='b',alpha=0.2,lw=0, label='Harmoniques')
            else:
                plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='b',alpha=0.2,lw=0)
            #plt.plot(((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6), (self.indice_harm[n], self.indice_harm[n]), 'b', linewidth=1.0)
        plt.scatter(self.freq_harm/1e6, self.indice_harm, c='b', s=100, marker='o')
            #plt.scatter(self.freq_harm[n]/1e6, self.indice_harm[n], c='b', s=100, marker='o')
        # Ultra-harmonics
        plt.plot(((self.freq_Uharm[1:-1]-self.delta_Uharm)/1e6 , (self.freq_Uharm[1:-1]-self.delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        plt.plot(((self.freq_Uharm[1:-1]+self.delta_Uharm)/1e6 , (self.freq_Uharm[1:-1]+self.delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
        for n in range(1,self.n_harm):
            #plt.axvspan((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6, color='w',alpha=None,lw=0)
            if n==0:
                plt.axvspan((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6, color='g',alpha=0.2,lw=0, label='Ultra-harmoniques')
            else:
                plt.axvspan((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6, color='g',alpha=0.2,lw=0)
            #plt.plot(((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6), (self.indice_Uharm[n], self.indice_Uharm[n]), 'g', linewidth=1.0)
            plt.scatter(self.freq_Uharm[n]/1e6, self.indice_Uharm[n], c='g', s=100, marker='o')

        plt.plot(((self.freq_Uharm[0]-self.delta_Uharm)/1e6 , (self.freq_Uharm[0]-self.delta_Uharm)/1e6), (0, Ymax), 'purple', linestyle = '--', linewidth=2.0)
        plt.plot(((self.freq_Uharm[0]+self.delta_Uharm)/1e6 , (self.freq_Uharm[0]+self.delta_Uharm)/1e6), (0, Ymax), 'purple', linestyle = '--', linewidth=2.0)
        plt.axvspan((self.freq_Uharm[0]-self.delta_Uharm)/1e6,(self.freq_Uharm[0]+self.delta_Uharm)/1e6, color='purple',alpha=0.2,lw=0, label='Sous-harmonique')
        plt.scatter(self.freq_Uharm[0]/1e6, self.indice_Uharm[0], c='purple', s=100, marker='o')
        plt.grid(True, which='major')
        plt.ylabel('Amplitude',fontsize=20)
        plt.yscale('log')
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel('Fréquence [MHz]',fontsize=20)
        legend = plt.legend(fontsize=20, framealpha=1.0, loc = 'upper right')
        legend.get_frame().set_edgecolor('black')  # Couleur de la bordure
        legend.get_frame().set_linewidth(1.5)  # Épaisseur de la bordure
        plt.title("Spectre d'un signal acoustique rétro-diffusé",fontsize=22, fontweight = 'bold')
        plt.xlim([self.freq_emiss/1e6/5, 7*1.01])#(harm_f[-2]-2*harm_df)/1e6])
        plt.ylim([Ymin, Ymax])
        #plt.xlim([self.freq_emiss/1e6/5, 9.5*1.01])#(harm_f[-2]-2*harm_df)/1e6])
        #plt.ylim([1e-3, Ymax])
        
        plt.savefig(chemin+'spec.png',bbox_inches='tight')
        plt.close("all")

        
    def plot_W(self,chemin,n_max=10000):
        for win in range(min(self.n_window,n_max)):
            plt.clf()
            spec=self.spectrogram[win,:self.size_window//2]#size_window
            Y=spec
            Ymax = 1.01*np.amax(spec)
            Ymin = 0.99*np.amin(spec)
            Ymax = 500
            Ymin = 1e-2
            
            plt.plot(self.frequencies/1e6, Y, 'k', label="Waveform", linewidth=2.0)
            # Fundamental
            plt.plot((self.freq_emiss/1e6 , self.freq_emiss/1e6), (0, Ymax), 'r-')
            # Inertial bandwidth
            for n in range(self.num_BB):
                #plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='w',alpha=None,lw=0)
                plt.axvspan((self.freq_BB[n]-self.delta_BB)/1e6,(self.freq_BB[n]+self.delta_BB)/1e6, color='r',alpha=0.2,lw=0)
                plt.plot(((self.freq_BB[n]-self.delta_BB)/1e6,(self.freq_BB[n]+self.delta_BB)/1e6), (self.indice_BB_w[win], self.indice_BB_w[win]), 'r', linewidth=2.0)
                #plt.plot(((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6), (self.indice_harm[n], self.indice_harm[n]), 'b', linewidth=1.0)
            #plt.axvspan(self.frequencies[self.freq_Uharm_n[1]]/1e6,self.frequencies[self.freq_harm_n[2]]/1e6, color='r',alpha=0.2,lw=0)  # ultraharm_f[0]+harm_df
    # =============================================================================
    #         # Sub-harmonic
    #         plt.plot(((subharm_f-harm_df)/1e6 , (subharm_f-harm_df)/1e6), (0, Ymax), 'y--', linewidth=2.0)
    #         plt.plot(((subharm_f+harm_df)/1e6 , (subharm_f+harm_df)/1e6), (0, Ymax), 'y--', linewidth=2.0)
    #         plt.axvspan((subharm_f-harm_df)/1e6,(subharm_f+harm_df)/1e6, color='y',alpha=0.2,lw=0)
    # =============================================================================
            # Harmonics
            plt.plot(((self.freq_harm-self.delta_harm)/1e6 , (self.freq_harm-self.delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
            plt.plot(((self.freq_harm+self.delta_harm)/1e6 , (self.freq_harm+self.delta_harm)/1e6), (0, Ymax), 'b--', linewidth=2.0)
            for n in range(self.n_harm+1):
                #plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='w',alpha=None,lw=0)
                plt.axvspan((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6, color='b',alpha=0.2,lw=0)
                plt.plot(((self.freq_harm[n]-self.delta_harm)/1e6,(self.freq_harm[n]+self.delta_harm)/1e6), (self.indice_harm_w[n,win], self.indice_harm_w[n,win]), 'b', linewidth=1.0)
            # Ultra-harmonics
            plt.plot(((self.freq_Uharm-self.delta_Uharm)/1e6 , (self.freq_Uharm-self.delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
            plt.plot(((self.freq_Uharm+self.delta_Uharm)/1e6 , (self.freq_Uharm+self.delta_Uharm)/1e6), (0, Ymax), 'g--', linewidth=2.0)
            for n in range(self.n_harm+1):
                #plt.axvspan((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6, color='w',alpha=None,lw=0)
                plt.axvspan((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6, color='g',alpha=0.2,lw=0)
                plt.plot(((self.freq_Uharm[n]-self.delta_Uharm)/1e6,(self.freq_Uharm[n]+self.delta_Uharm)/1e6), (self.indice_Uharm_w[n,win], self.indice_Uharm_w[n,win]), 'g', linewidth=1.0)
            plt.grid(True, which='major')
            plt.ylabel('Magnitude',fontsize=20)
            plt.yscale('log')
            plt.xlabel('Frequency [MHz]',fontsize=20)
            plt.title('Spectre moyenne du pulse avec visualisation des differentes bandes frequentielle',fontsize=20, fontweight = 'bold')
            plt.xlim([self.freq_emiss/1e6/5, (self.freq_harm[self.n_harm-1]+self.delta_harm)/1e6*1.01])#(harm_f[-2]-2*harm_df)/1e6])
            plt.ylim([ Ymin, Ymax])
            #plt.xlim([self.freq_emiss/1e6/5, 9.5*1.01])#(harm_f[-2]-2*harm_df)/1e6])
            #plt.ylim([1e-3, Ymax])
            
            plt.savefig(chemin+'spec_{}.png'.format(win),bbox_inches='tight')

class experiment():
    
    def __init__(self, fe, f0,start:int,end:int, size_window : int = 2048, n_harm : int = 4, n_BB : int =5, delta_harm =50e3,delta_Uharm =50e3, spacer = 50e3, peak_detection : bool =True, size_decal : int = 2048):
        #valeurs initiales de l'experience
        self.freq_ech = float(fe)
        self.freq_emiss = float(f0)
        self.num_BB = n_BB
        self.n_harm=n_harm
        self.start=start
        self.end=end
        self.size_decal=size_decal
        self.freq_harm=np.array([self.freq_emiss*i for i in range(1,n_harm+3)])
        self.freq_Uharm=np.array([self.freq_emiss*(i-0.5) for i in range(1,n_harm+3)])
        self.freq_BB=np.array([self.freq_emiss*(i/2.+1.75) for i in range(self.num_BB)])
        self.peak_detection=peak_detection
        self.size_window=size_window
        self.frequencies= np.arange(self.size_window//2) * self.freq_ech/self.size_window
        self.delta_harm=delta_harm
        self.delta_Uharm=delta_Uharm
        self.delta_BB=0.25*self.freq_emiss-max(self.delta_harm,self.delta_Uharm)-spacer
        self.pulses=[]
        self.n_pulse=0
        self.ratio=0
    
    def add_pulse(self,data,spacer=50e3):
        pulse=temp_sample(self.freq_ech,self.freq_emiss,data,spacer = spacer, size_window= self.size_window, n_harm= self.n_harm,delta_harm =self.delta_harm,delta_Uharm =self.delta_Uharm,peak_detection =self.peak_detection,start=self.start,end=self.end,size_decal=self.size_decal)
        self.pulses.append(pulse)
        self.n_pulse+=1
    
    def add_pulses(self,array_data,spacer=50e3):
        for data in array_data:
            pulse=temp_sample(self.freq_ech,self.freq_emiss,data,spacer = spacer, size_window= self.size_window, n_harm= self.n_harm,delta_harm =self.delta_harm,delta_Uharm =self.delta_Uharm,peak_detection =self.peak_detection,start=self.start,end=self.end,size_decal=self.size_decal)
            self.pulses.append(pulse)
            self.n_pulse+=1
    
    def plot_windowed(self,chemin,nbit,fit,mini_value=0,maxi_value=100,n=1,legend = ['data '+str(int(i))for i in range(1,11)], ramp = False):
        y_axis = 'Pressure (KPa)'
        if ramp:
            y_axis = 'Différents tirs'
            
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##recuperation des indices
        n_plot = 4
        plot=np.zeros((n_plot,self.n_pulse//n,self.pulses[0].n_window))
        
        for j in range(self.n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
            val = np.mean(np.mean([self.pulses[j*n+a].indice_harm_w[1:3] for a in range(n)] ,axis=0),axis=0)
            plot[0,j,:]=20*np.log10(val) #20*np.log10
            
            val = np.mean(np.mean([self.pulses[j*n+a].indice_Uharm_w[1:3] for a in range(n)] ,axis=0),axis=0)
            plot[1,j,:]=20*np.log10(val) #20*np.log10
            
            val = np.mean(np.mean([self.pulses[j*n+a].indice_Uharm_norm_div_w[1:4] for a in range(n)] ,axis=0),axis=0)
            plot[2,j,:]=20*np.log10(val) #20*np.log10

            val = np.mean(np.mean([self.pulses[j*n+a].indice_BB_sliced_w[0:2] for a in range(n)] ,axis=0),axis=0)
            #val = np.mean([self.pulses[j*n+a].indice_BB_w for a in range(n)] ,axis=0)
            plot[3,j,:]=20*np.log10(val) #20*np.log10

        fig=plt.figure(figsize=(20,11))
        amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
        pression=[fit[0]*j+fit[1] for j in amp]
        self.n_window = self.pulses[0].n_window
        nom_str = ["\\H_n{}_perLow{}_perHigh{}_ramp.png".format(n,mini_value,maxi_value),"\\UH_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\UH_norm_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\BB_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value)]
        title_str = ["Composante harmonique ","Composante ultra-harmonique ","Composante ultra-harmonique émergente ","Composante inertielle "]
        for j in range(n_plot):
            if j == 2:
                continue
            plt.clf()
            plt.title(title_str[j]+ legend,fontsize=22, fontweight = 'bold')
            mini=np.percentile(plot[j],mini_value)
            maxi=np.percentile(plot[j],maxi_value)
            #print(np.shape(plot[j])[1]*self.size_window/self.freq_ech*1000.)
            #print((nbit[1])*fit[0]+fit[1])
            #print((nbit[0])*fit[0]+fit[1])
            #plt.imshow(plot[j,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[j])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
            plt.imshow(plot[j,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[j])[1]*self.size_decal/self.freq_ech*1000.,self.n_pulse,nbit[0]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
            #plt.xlabel('Intra-Pulse Time'.upper()+' (ms)',fontsize=20)
            plt.xlabel('Temps au sein du tir (ms)',fontsize=20)
            plt.ylabel(y_axis,fontsize=20)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            
            plt.colorbar()
            #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
            plt.savefig(chemin+nom_str[j],bbox_inches='tight')
        plt.close("all")
        

        # path(chemin)
        # plot=np.zeros((self.n_pulse,self.pulses[0].n_window))
        
        # for j in range(self.n_pulse):#self.exp[i].pulses[j].indice_harm_w[1]+
        #     val = self.pulses[j].indice_BB_w
        #     plot[j,:]=20*np.log10(val) #20*np.log10

        # self.ratio=plot
        # plt.figure(figsize=(20,11))
        # self.n_window = self.pulses[0].n_window

        # nom = "\\BB_moy_windowed_perLow{}_perHigh{}.png".format(mini_value,maxi_value)
        # mini=np.percentile(plot,mini_value)
        # maxi=np.percentile(plot,maxi_value)
        # plt.imshow(plot[:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot)[1]*self.size_window/self.freq_ech*1000.,self.n_pulse,0], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
        # plt.title("INERTIAL cavitation components ".upper() + ' logarithmic',fontsize=30, fontweight = 'bold')
        # plt.xlabel('Intra-Pulse Time'.upper()+' (ms)',fontsize=20)
        # plt.ylabel('Pressure (KPa)',fontsize=20)
        
        # plt.savefig(chemin+nom,bbox_inches='tight')
        # plt.close("all")

    def plot_comp_windowed(self,chemin,mini_value=0,maxi_value=100):
        y_axis = '#tirs'
        n=1  
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##recuperation des indices
        n_plot = 3
        n_BB = self.num_BB
        n_harm=self.n_harm+1
        long=[n_harm,n_harm,n_BB,n_harm,n_harm]
        n_comp = max(self.pulses[0].n_harm,self.pulses[0].num_BB)
        plot=np.zeros((n_plot,self.n_pulse,n_comp,self.pulses[0].n_window))
        
        for j in range(self.n_pulse):#self.exp[i].pulses[j].indice_harm_w[1]+
            val = self.pulses[j].indice_harm_w[:]
            plot[0,j,:val.shape[0]]=20*np.log10(val) #20*np.log10

            val = self.pulses[j].indice_Uharm_w[:]
            plot[1,j,:val.shape[0]]=20*np.log10(val) #20*np.log10
            
            val = self.pulses[j].indice_BB_sliced_w[:]
            plot[2,j,:val.shape[0]]=20*np.log10(val) #20*np.log10

        self.n_window = self.pulses[0].n_window
        nom_str = ["\\separates_H_n{}_perLow{}_perHigh{}_ramp.png".format(n,mini_value,maxi_value),"\\separates_UH_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\separates_UH_norm_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\BB_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value)]
        title_str = ["Harmonic components","Ultra-harmonic components","Composante ultra-harmonique émergente","Inertial components"]
        title_str = ["Indices harmoniques","Indices ultra-harmoniques","Indices inertiels","Indices inertiels"]
        legend = [["H{}".format(i) for i in range(1,n_harm+1)], ["UH{}".format(i) for i in range(1,n_harm+1)], ["BB{}".format(i) for i in range(1,n_BB+1)], ["H{}".format(i) for i in range(1,n_harm+1)], ["BB{}".format(i) for i in range(1,n_BB+1)]]
        for j in range(n_plot):
            fig, axes = plt.subplots(nrows=1, ncols=long[j], figsize=(20,11))
            fig.suptitle(title_str[j],fontsize=25, fontweight = 'bold',y=0.94)
            for i,ax in enumerate(axes.flat):
                mini=np.percentile(plot[j,:,i,:],mini_value)
                maxi=np.percentile(plot[j,:,i,:],maxi_value)
                im = ax.imshow(plot[j,:,i,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,self.n_window*self.size_decal/self.freq_ech*1000.,self.n_pulse,1], vmin=mini,vmax=maxi) 
                ax.title.set_text(legend[j][i].upper())
                ax.title.set_fontweight('bold')
                ax.title.set_fontsize(18)
                ax.set_xlabel('Temps intra'+' (ms)',fontsize=20)
                if i==0:
                    ax.set_ylabel(y_axis,fontsize=30)
                ax.tick_params(axis = 'both',which = 'both',  labelsize=16)
                
            #je veux grossir la taille de la police des ticks de la colorbar
            colorbar = fig.colorbar(im, ax=axes.ravel().tolist())
            colorbar.ax.yaxis.set_tick_params(labelsize=20)  # Taille des ticks
            colorbar.ax.set_ylabel('dB', fontsize=23, rotation=90, labelpad=15, fontweight='bold')  # Label
            plt.savefig(chemin+nom_str[j],bbox_inches='tight')
            plt.close("all")


        fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(20,11))
        fig.suptitle("Composante proche du fondamental",fontsize=25, fontweight = 'bold',y=0.94)
        different_carte = [plot[0,:,0,:], plot[0,:,1,:], plot[1,:,0,:], plot[1,:,1,:], plot[2,:,0,:]]
        titres = ["H1","H2","UH1","UH2","BB1"]
        for i,ax in enumerate(axes.flat):
            mini=np.percentile(different_carte[i] ,mini_value)
            maxi=np.percentile(different_carte[i] ,maxi_value)
            im = ax.imshow(different_carte[i] ,aspect='auto',interpolation='none',cmap='turbo',extent=[0,self.n_window*self.size_decal/self.freq_ech*1000.,self.n_pulse,1], vmin=mini,vmax=maxi) 
            ax.title.set_text(titres[i].upper())
            ax.title.set_fontweight('bold')
            ax.title.set_fontsize(18)
            ax.set_xlabel('Temps intra'+' (ms)',fontsize=20)
            if i==0:
                ax.set_ylabel(y_axis,fontsize=30)
            ax.tick_params(axis = 'both',which = 'both',  labelsize=16)
            
        #je veux grossir la taille de la police des ticks de la colorbar
        colorbar = fig.colorbar(im, ax=axes.ravel().tolist())
        colorbar.ax.yaxis.set_tick_params(labelsize=20)  # Taille des ticks
        colorbar.ax.set_ylabel('dB', fontsize=23, rotation=90, labelpad=15, fontweight='bold')  # Label
        plt.savefig(chemin+"\\Composantes_proche_fonda.png",bbox_inches='tight')
        plt.close("all")


    def plot_indice(self,nom,chemin,amp = False):
        path(chemin)
        chemin_log=chemin+'\\log\\'
        chemin_lin=chemin+'\\linear\\'
        path(chemin_log)
        path(chemin_lin)
        plot=np.zeros((4,self.n_pulse))
        n_plot = 4
        temp = np.zeros((2,n_plot))
        for j in range(self.n_pulse):
            plot[0,j] = np.mean(self.pulses[j].indice_harm[1:-1]) #np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
            plot[1,j] = np.mean(self.pulses[j].indice_Uharm[1:4]) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
            plot[2,j] = np.mean(self.pulses[j].indice_BB) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
            plot[3,j] = 3 * plot[1,j] - 20 * plot[2,j] 
        
        plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit de bande","Indices représentant les 3xUH-30xBB","Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        nom_img=["harm","U_harm","BB","ratio_"+nom+"_all","ratio_"+nom+"_start","ratio_"+nom+"_end"]
        colors=['black','red','peru','forestgreen','dodgerblue','gold']
        pression = np.arange(self.n_pulse)
        fit = np.array([1, 0])
        redblue=['b','g','r','teal','black','black','black']
        # Stable cavitation dose
        
        for i in range(0,4):
            fig, host = plt.subplots(figsize=(20,11)) # (width, height) in inches
            host.plot(pression,plot[i,:],  c=redblue[i])
            if amp != False:
                par1 = host.twinx()
                par1.plot(pression,amp,color='gray')
            
            plt.title(plot_legend[i],color=redblue[i],fontsize=30, fontweight = 'bold')
            host.set_xlabel('Pression (KPa)',fontsize=20)
            host.set_ylabel('Magnitude [a.u.]', color=redblue[i],fontsize=20)
            par1.set_ylabel('Pulse amplitude [bit IGT]', color=redblue[i],fontsize=20)
            Ymax = 1.01*np.amax(plot[i,:])
            Ymin = 0.99*np.amin(plot[i,:])
            host.set_ylim([Ymin, Ymax])  
            par1.set_ylim([0.99*np.amin(amp),1.01*np.amax(amp)])    
            host.grid(True)
            fig.tight_layout()
            host.set_yscale('linear')
            plt.savefig(chemin_lin+nom_img[i]+'_linear.png',bbox_inches='tight')  
            host.set_yscale('log')
            plt.savefig(chemin_log+nom_img[i]+'_linear.png',bbox_inches='tight') 
            plt.close("all")
         
    def plot_indice_bis(self,nom,title,chemin,n,nbit,legend,fit = False,fit_indic = False):
        path(chemin)
        n_plot=4
        plot=np.zeros((n_plot,self.n_pulse//n,3))
        plot_raw=np.zeros((n_plot,self.n_pulse//n,3))

        temp = np.zeros((n_plot,n))
        for j in range(self.n_pulse//n):
            temp[0] =[np.mean(self.pulses[j*n+a].indice_harm[1:3]) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
            temp[1] =[np.mean(self.pulses[j*n+a].indice_Uharm[1:3]) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
            temp[2] =[np.mean(self.pulses[j*n+a].indice_BB) for a in range(n)] #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
            temp[3] =[np.mean(self.pulses[j*n+a].indice_Uharm_norm_div[1:4]) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
            for k in range(n_plot):
                plot_raw[k,j,0]=np.mean(temp[k])
                plot_raw[k,j,1]=plot_raw[k,j,0]+np.std(temp[k])
                plot_raw[k,j,2]=plot_raw[k,j,0]-np.std(temp[k])

        for k in range(n_plot):
            plot[k]=plot_raw[k]-np.min(plot_raw[k,:,0])
        for k in range(n_plot):
            plot[k]=plot[k]/np.max(plot[k,:,0])
            
        plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit large bande","Indices représentant les 3xUH-30xBB","Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        label = ["Harmonics", "Ultraharmonics", "Broadband noise"]
        label = ["Harmoniques", "Ultra-harmoniques", "Bruit large bande"]
        nom_img=["harm","U_harm","BB","ratio_"+nom+"_all","ratio_"+nom+"_start","ratio_"+nom+"_end"]
        colors=['blue','green','red','forestgreen','dodgerblue','gold']
        pression = np.arange(self.n_pulse)
        if fit_indic==False : 
            pression = np.arange(self.n_pulse//n)
            fit = np.array([1, 0])
        else :
            amp=np.arange(nbit[0],nbit[1])
            pression=[fit[0]*j+fit[1] for j in amp]
        redblue=['b','g','r','teal','black','black','black']
        # Stable cavitation dose
        
        plt.figure(figsize=(20,11)) # (width, height) in inches
        for i in range(0,3):
            plt.plot(pression,plot[i,:,0], c=redblue[i],label=label[i])
            plt.fill_between(pression,plot[i,:,1],plot[i,:,2], color=redblue[i],alpha=0.2)
            Ymax = 1.01*np.amax(plot[i,:])
            Ymin = 1.01*np.amin(plot[i,:])
            plt.ylim([Ymin, Ymax]) 
            plt.ylim([0, 1.08])
            plt.grid(True) #Normalized indexes
        plt.title("Indices normalisés       "+title,color='black', fontsize=30, fontweight = 'bold')
        plt.legend(fontsize=20, loc = 'upper left')
        plt.xlim([0, np.max(pression)]) 
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        if fit_indic == False :
            plt.xlabel('#tirs',fontsize=25)
        else :
            plt.xlabel('Pression (kPa)',fontsize=25)
        plt.ylabel('Amplitude normalisée', color=redblue[i],fontsize=25)
        plt.savefig(chemin+nom+'.png',bbox_inches='tight') 
        plt.close("all")    
        
    def plot_indice_RAMP(self,nom,title,chemin,pression):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        chemin_lin=chemin
        n_window = self.pulses[0].n_window
        ##HARMONIQUES
        n_plot = 3
        plot_raw=np.zeros((n_plot,n_window,3))
        plot=np.zeros((n_plot,n_window,3))
        temp = np.zeros((n_plot,self.n_pulse,n_window))
        for j in range(self.n_pulse):
            temp[0,j] =np.mean(self.pulses[j].indice_harm_w[1:-1],axis=0)#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
            temp[1,j] =np.mean(self.pulses[j].indice_Uharm_w[1:4],axis=0) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
            temp[2,j] =self.pulses[j].indice_BB_w  #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
            #temp[3,j] =np.mean(self.pulses[j].indice_Uharm_norm_div_w[1:4],axis=0) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
        for k in range(n_plot):
            plot_raw[k,:,0]=np.mean(temp[k], axis=0)
            plot_raw[k,:,1]=plot_raw[k,:,0]+np.std(temp[k], axis=0)
            plot_raw[k,:,2]=plot_raw[k,:,0]-np.std(temp[k], axis=0)
        for k in range(n_plot):
            plot[k]=plot_raw[k]-np.min(plot_raw[k,:,0])
        for k in range(n_plot):
            plot[k]=plot[k]/np.max(plot[k,:,0])
            
            
        plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit large bande","Indices représentant les ultra-harmoniques-BB","Indice representant UH-{}BB".format(7),"Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        nom_img=["harm","U_harm","BB","U_harm_norm","UH_BB","ratio_"+title+"_all","ratio_"+title+"_start","ratio_"+title+"_end"]
        colors=['black','red','peru','forestgreen','dodgerblue','gold']
        label = ["Harmoniques", "Ultra-harmoniques", "Bruit large bande"]
        redblue=['b','g','r','teal','m','black','black','black']
        # Stable cavitation dose
        fig=plt.figure(figsize=(20,11))
        for i in range(n_plot):
            plt.plot(pression,plot[i,:,0], c=redblue[i],label=label[i])
            plt.fill_between(pression,plot[i,:,1],plot[i,:,2], color=redblue[i],alpha=0.2)
        plt.legend(fontsize=20, loc = 'upper left')
        plt.title("Indices normalisés       "+title,color="black",fontsize=22, fontweight = 'bold')
        plt.xlabel('Pression (kPa)',fontsize=20)
        plt.ylabel('Amplitude normalisée',fontsize=20)
        #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')
        Ymax = 1.01*np.amax(plot[:,:,1])
        Ymin = 0.99*np.amin(plot[:,:,2])
        plt.ylim([-0.1, 1.1])    
        plt.grid(True)
        plt.tight_layout()    
        plt.xlim([0, np.max(pression)]) 
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.yscale('linear')
        plt.savefig(chemin_lin+nom+'_ramp.png',bbox_inches='tight')           
        plt.close("all")  
        
    def plot_indice_windowed(self,chemin,mapval,bits):
# =============================================================================
#         if not os.path.isdir(chemin): # check if folder exists, otherwise create it
#             os.mkdir(chemin)
# =============================================================================
        fig=plt.figure(figsize=(20,11))
        fit=np.array([29.95777559,  6.81690626])
        amp=np.arange(bits[0],(bits[1]-1)*30+1)
        pression=[fit[0]*j+fit[1] for j in amp]
        self.n_window = self.pulses[0].n_window
        ##HARMONIQUES
        plot=np.zeros((3,self.n_pulse,self.n_window))
        for i,pulse in enumerate(self.pulses):
            plot[0,i] = 20*np.log10(np.mean(pulse.indice_harm_w, axis = 0))
            plot[1,i] = 20*np.log10(np.mean(pulse.indice_Uharm_w, axis = 0 ))
            plot[2,i] = 20*np.log10(pulse.indice_BB_w)
        plot_legend=["Harmoniques","Ultra-harmoniques","Bruit large bande"]
        print("min : {}".format(np.min(plot[1])))
        print("max : {}".format(np.max(plot[1])))
        redblue=['b','g','r']
        
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20,11))
        for i,ax in enumerate(axes.flat):
            print(np.shape(plot[i]))
            im = ax.imshow(plot[i][:23*30+1,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[i])[1]*self.size_window/self.freq_ech*1000.,(23-1)*fit[0]+fit[1],0], vmin=mapval[i*2],vmax=mapval[i*2+1])  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
            ax.title.set_text(plot_legend[i].upper() + ' (logarithmic scale)')
            ax.title.set_color(redblue[i])
            ax.title.set_fontweight('bold')
# =============================================================================
#             axx=ax.get_xticks()
#             axy=ax.get_yticks()
#             x_ax=[i*self.size_window/self.freq_ech*1000 for i in axx]
#             y_ax=[i//30*fit[0]+fit[1] for i in axy]
#             ax.set_xticks(x_ax[1:-1])
#             ax.set_yticks(y_ax[1:-1])
# =============================================================================
            ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
            if i==0:
                ax.set_ylabel('Pressure (KPa)')
        
        fig.colorbar(im,ax=axes.ravel().tolist())
        plt.savefig(chemin,bbox_inches='tight')
        plt.close("all")
        return [np.min(plot[0,:]),np.max(plot[0,:]),np.min(plot[1,:]),np.max(plot[1,:]),np.min(plot[2,:]),np.max(plot[2,:])]
        
    def plot_indice_component(self,chemin,nbit,fit,n=1,fit_indic=False):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##HARMONIQUES
        if fit_indic==False : 
            pression = np.arange(self.n_pulse//n)
            fit = np.array([1, 0])
        else :
            amp=np.arange(nbit[0],nbit[1])
            pression=[fit[0]*j+fit[1] for j in amp]
        n_plot = 4
        n_BB = self.num_BB
        n_harm=self.n_harm+1
        plot = np.zeros((n_plot,self.n_pulse//n,n_harm,3))
        temp = np.zeros((n_plot,n,n_harm))
        for i in range(self.n_pulse//n):
            temp[1] =[self.pulses[i*n+a].indice_harm for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
            temp[2] =[self.pulses[i*n+a].indice_Uharm for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
            temp[3] =[self.pulses[i*n+a].indice_Uharm_norm_div for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+s
            temp[0] =[self.pulses[i*n+a].indice_BB_sliced for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
            for k in range(0,n_plot):
                for j in range(n_harm):
                    plot[k,i,j,0]=np.mean(temp[k,:,j])
                    plot[k,i,j,1]=plot[k,i,j,0]+np.std(temp[k,:,j])
                    plot[k,i,j,2]=plot[k,i,j,0]-np.std(temp[k,:,j])
            # plot[0,i,0,0]=np.mean(temp[0,:,0])
            # plot[0,i,0,1]=plot[0,i,0,0]+np.std(temp[0,:,0])
            # plot[0,i,0,2]=plot[0,i,0,0]-np.std(temp[0,:,0])
        plot_legend=["Indices représentant le bruit large bande","Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant les ultra-harmoniques émergents","Indice representant 3UH-20BB","Indices représentant le fondamental & les harmoniques",]
        legend=[["Bruit large bande {}".format(idbb+1) for idbb in range(n_BB)],["error","2f0","3f0","4f0","5f0","6f0","7f0"],["f0/2","3f0/2","5f0/2","7f0/2","9f0/2","11f0/2"],["f0/2","3f0/2","5f0/2","7f0/2","9f0/2","11f0/2"],["f0","2f0","3f0","4f0","5f0","6f0","7f0"],["f0/2","3f0/2","5f0/2","7f0/2","9f0/2","11f0/2"]]
        colors=['black','red','peru','forestgreen','dodgerblue','indigo']
        long=[n_BB,n_harm,n_harm,n_harm,n_harm]
        deb_ind=[0,1,1,1,1]
        redblue=['b','g','r','teal','m','black','black','black']
        # Stable cavitation dose
        fig=plt.figure(figsize=(20,11))
        for i in range(3):
            fig.clf()
            Ymax = -9999999
            Ymin = 999999999
            Ymax_zoom = -9999999
            Ymin_zoom = 999999999
            zomm_min = 50
            zoom_max = 800
            for j in range(deb_ind[i],long[i]):
                # if j==0 and i==1:
                #      continue
                plt.plot(pression,plot[i,:,j,0], c=colors[j],label=legend[i][j])
                plt.fill_between(pression,plot[i,:,j,1],plot[i,:,j,2], color=colors[j],alpha=0.2)
                
                if Ymax<1.01*np.amax(plot[i,:,j,1]):
                    Ymax = 1.01*np.amax(plot[i,:,j,1])
                if Ymin>0.99*np.amin(plot[i,:,j,2]):
                    Ymin = 0.99*np.amin(plot[i,:,j,2])
                if Ymax_zoom<1.01*np.amax(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,1]):
                    Ymax_zoom = 1.01*np.amax(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,1])
                if Ymin_zoom>0.99*np.amin(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,2]):
                    Ymin_zoom = 0.99*np.amin(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,2])
            plt.title(plot_legend[i],color=redblue[i],fontsize=30, fontweight = 'bold')
            if fit_indic == False :
                plt.xlabel('#tirs',fontsize=25)
            else :
                plt.xlabel('Pression (kPa)',fontsize=25)
            plt.legend(fontsize=24, loc = 'upper left')
            plt.ylabel('Amplitude [a.u.]', color=redblue[i],fontsize=25)
            #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')

            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.xlim([0, np.max(pression)])  
            plt.ylim([Ymin, Ymax])    
            plt.grid(True)
            plt.tight_layout()
            plt.yscale('linear')
            plt.savefig(chemin+"\\"+plot_legend[i]+'.png',bbox_inches='tight')
            #plt.xlim([zomm_min, zoom_max])
            #plt.ylim([Ymin_zoom, Ymax_zoom])
            #plt.savefig(chemin+plot_legend[i]+'_linear_ZOOM.png',bbox_inches='tight')
        plt.close("all")
        
        
        # val_harm = np.zeros((self.n_pulse//n,n,n_harm))
        # val_Uharm = np.zeros((self.n_pulse//n,n,n_harm))
        # val_BB = np.zeros((self.n_pulse//n,n))
        # for j in range(self.n_pulse//n):
        #     val_harm[j] =[self.pulses[j*n+a].indice_harm for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
        #     val_Uharm[j] =[self.pulses[j*n+a].indice_Uharm for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
        #     val_BB[j] =[self.pulses[j*n+a].indice_BB for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
        
        # for i in range(self.n_pulse//n):
        #     for j in range(n_harm):
        #         plot[0,i,j,0]=np.mean(val_harm[i,:,j])
        #         plot[0,i,j,1]=plot[0,i,j,0]+np.std(val_harm[i,:,j])
        #         plot[0,i,j,2]=plot[0,i,j,0]-np.std(val_harm[i,:,j])
        #     for j in range(n_harm):
        #         plot[1,i,j,0]=np.mean(val_Uharm[i,:,j])
        #         plot[1,i,j,1]=plot[1,i,j,0]+np.std(val_Uharm[i,:,j])
        #         plot[1,i,j,2]=plot[1,i,j,0]-np.std(val_Uharm[i,:,j])
        #     plot[2,i,0,0]=np.mean(val_BB[i,:])
        #     plot[2,i,0,1]=plot[2,i,0,0]+np.std(val_BB[i,:])
        #     plot[2,i,0,2]=plot[2,i,0,0]-np.std(val_BB[i,:])
        # plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit de bande","Indices représentant le fondamental & les harmoniques",]
        # legend=[["error","2f0","3f0","4f0","5f0","6f0","7f0"],["f0/2","3f0/2","5f0/2","7f0/2","9f0/2","11f0/2"],["Broadband"],["f0","2f0","3f0","4f0","5f0","6f0","7f0"]]
        # colors=['black','red','peru','forestgreen','dodgerblue','indigo']
        # #fit = np.array([4.66745471*2, 5.80567673])
        # amp=np.arange(nbit[0],nbit[1])
        # pression=[fit[0]*j+fit[1] for j in amp]
        # long=[n_harm,n_harm,1]
        # redblue=['b','g','r']
        # x_abs=np.arange(1,71)
        # # Stable cavitation dose
        # fig=plt.figure(figsize=(20,11))
        # for i_bis in range(4):
        #     i=i_bis
        #     if i==3:
        #         i=0
        #     fig.clf()
        #     Ymax = -9999999
        #     Ymin = 999999999
        #     Ymax_zoom = -9999999
        #     Ymin_zoom = 999999999
        #     zomm_min = 50
        #     zoom_max = 900
        #     for j in range(long[i]):
        #         if j==0 and i_bis==0:
        #             continue
        #         #print(i_bis,i,j)
        #         plt.plot(pression,plot[i,:,j,0], c=colors[j],label=legend[i_bis][j])
        #         plt.fill_between(pression,plot[i,:,j,1],plot[i,:,j,2], color=colors[j],alpha=0.2)
                
        #         #plt.plot(pression,plot[i,j], c=colors[j],label = legend[i][j])
        #         if Ymax<1.01*np.amax(plot[i,:,j,1]):
        #             Ymax = 1.01*np.amax(plot[i,:,j,1])
        #         if Ymin>0.99*np.amin(plot[i,:,j,2]):
        #             Ymin = 0.99*np.amin(plot[i,:,j,2])
        #         if Ymax_zoom<1.01*np.amax(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,1]):
        #             Ymax_zoom = 1.01*np.amax(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,1])
        #         if Ymin_zoom>0.99*np.amin(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,2]):
        #             Ymin_zoom = 0.99*np.amin(plot[i,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),j,2])
        #     plt.title(plot_legend[i_bis],color=redblue[i],fontsize=30, fontweight = 'bold')
        #     plt.xlabel('Pression KPa',fontsize=20)
        #     plt.legend()
        #     plt.ylabel(plot_legend[i_bis]+' [a.u.]', color=redblue[i],fontsize=20, fontweight = 'bold')
        #     #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')

        #     plt.ylim([Ymin, Ymax])    
        #     plt.grid(True)
        #     plt.tight_layout()
        #     plt.yscale('log')
        #     plt.savefig(chemin_log+plot_legend[i_bis]+'_log.png',bbox_inches='tight')
        #     plt.yscale('linear')
        #     plt.savefig(chemin_lin+plot_legend[i_bis]+'_linear.png',bbox_inches='tight')
        #     plt.xlim([zomm_min, zoom_max])
        #     plt.ylim([Ymin_zoom, Ymax_zoom])
        #     plt.savefig(chemin_lin+plot_legend[i_bis]+'_linear_ZOOM.png',bbox_inches='tight')
        #     plt.yscale('log')
        #     plt.savefig(chemin_log+plot_legend[i_bis]+'_log_ZOOM.png',bbox_inches='tight')
        # plt.close("all")
           
    def plot_indice_RAMP_std(self,nom,chemin,still_wind = 10, std_tresh = 3, true_harm = False, plot_true =True, val_std = 0.12, val_mean = 0.0, BB_mult = 1, UH_mult = 1):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        pression = [i for i in range(self.pulses[0].n_window)]
        chemin_lin=chemin
        n_window = self.pulses[0].n_window
        ##HARMONIQUES
        n_plot = 3
        plot_raw=np.zeros((n_plot,n_window,3))
        plot=np.zeros((n_plot,n_window,3))
        temp = np.zeros((n_plot,self.n_pulse,n_window))
        temp_norm = np.zeros((n_plot,self.n_pulse,n_window))
        for j in range(self.n_pulse):
            #temp[0,j] =np.mean(self.pulses[j].indice_harm_w[1:-1],axis=0)#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
            temp[0,j] =self.pulses[j].fondamental_w#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
            temp[1,j] =np.mean(self.pulses[j].indice_Uharm_w[1:4],axis=0) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
            temp[2,j] =self.pulses[j].indice_BB_w  #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
            
        for k in range(n_plot):
            plot_raw[k,:,0]=np.mean(temp[k], axis=0)
            plot_raw[k,:,1]=plot_raw[k,:,0]+np.std(temp[k], axis=0)
            plot_raw[k,:,2]=plot_raw[k,:,0]-np.std(temp[k], axis=0)
        for k in range(n_plot):
            plot[k]=plot_raw[k]-np.min(plot_raw[k,:,0])
        for k in range(n_plot):
            plot[k]=plot[k]/np.max(plot[k,:,0])
        

        for k in range(n_plot):
            for j in range(self.n_pulse):
                temp_norm[k,j]=temp[k,j]-np.min(temp[k,j])
                
                
        if true_harm :
            maxima = np.max(temp_norm[0,:,:], axis = 1)
            #print(np.shape(maxima))
        else :
            maxima = np.ones((self.n_pulse))
            #(np.shape(maxima))
        
        for j in range(self.n_pulse):
            temp_norm[0,j]=temp_norm[0,j]/np.max(temp_norm[0,j]) * maxima[j]
        for k in range(1,n_plot):
            for j in range(self.n_pulse):
                eta_0 = temp_norm[k,j]*val_std/np.std(temp_norm[k,j,:still_wind])/std_tresh
                eta_1 = eta_0-(np.mean(eta_0[:still_wind])-val_mean)
                temp_norm[k,j]=eta_1
        
        time = [i*self.size_decal*1000/self.freq_ech for i in range(n_window)]
        plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit large bande","Indices représentant les ultra-harmoniques-BB","Indice representant UH-{}BB".format(7),"Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        nom_img=["harm","U_harm","BB","U_harm_norm","UH_BB","ratio_"+nom+"_all","ratio_"+nom+"_start","ratio_"+nom+"_end"]
        colors=['black','red','peru','forestgreen','dodgerblue','gold']
        label = ["Fondamental", "Ultra-harmoniques", "Bruit large bande"]
        redblue=['b','g','r','teal','m','black','black','black']
        # Stable cavitation dose
        fig=plt.figure(figsize=(20,11))
        for i in range(n_plot):
            plt.plot(pression,plot[i,:,0], c=redblue[i],label=label[i])
            plt.fill_between(pression,plot[i,:,1],plot[i,:,2], color=redblue[i],alpha=0.2)
        plt.legend(fontsize=20)
        plt.title("Différentes composantes : "+nom,color="black",fontsize=22, fontweight = 'bold')
        plt.xlabel('Fenêtres temporelles',fontsize=20)
        plt.ylabel('Amplitude normalisée',fontsize=20)
        #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')
        Ymax = 1.01*np.amax(plot[:,:,1])
        Ymin = 0.99*np.amin(plot[:,:,2])
        plt.ylim([-0.01, 1.01])    
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('linear')
        plt.savefig(chemin_lin+'00'+nom+'.png',bbox_inches='tight')     
        plt.clf()
        uh,bb,both = 0, 0, 0
        win_BB = []
        win_UH = []
        amplit_rela_UH = []
        UH_abs = []
        UH_abs_j = []
        amplit_rela_BB = []
        BB_abs = []
        BB_abs_j = []
        harm_amp_UH = []
        harm_amp_BB = []
        for k in range(self.n_pulse):#
            signal =temp_norm[1:,k,:]
            signal_UH = signal[0]
            signal_BB = signal[1]
            max_UH, w_ind_UH = sortir_bruit(signal_UH, std_tresh, still_wind)
            max_BB, w_ind_BB = sortir_bruit(signal_BB, std_tresh, still_wind)
            #print(np.shape(signal))
            # print(" BB : max : {}   , pos : {}".format(max_BB,w_ind_BB))
            # print(" UH : max : {}   , pos : {}".format(max_UH,w_ind_UH))
            add = ""
            if max_UH>0 or max_BB>0:
                add = "event_"
                if max_UH>0:
                    moy_UH = np.mean(signal[0,:w_ind_UH])
                    std_UH = np.std(signal[0,:w_ind_UH])
                    win_UH.append(win_UH)
                    harm_amp_UH.append(temp[0,k,w_ind_UH])
                if max_BB>0:
                    moy_BB = np.mean(signal[1,:w_ind_BB])
                    std_BB = np.std(signal[1,:w_ind_BB])
                    win_BB.append(win_BB)
                    harm_amp_BB.append(temp[0,k,w_ind_BB])
                if max_UH>0 and max_BB>0:
                    colo = 'black'
                    both +=1
                    amplit_rela_UH.append(max_UH)
                    UH_abs.append(k)
                    UH_abs_j.append(w_ind_BB)
                    amplit_rela_BB.append(max_BB)
                    BB_abs.append(k)
                    BB_abs_j.append(w_ind_BB)
                    
                elif max_BB>0:
                    colo = 'r'
                    bb +=1
                    amplit_rela_BB.append(max_BB)
                    BB_abs.append(k)
                    BB_abs_j.append(w_ind_BB)
                
                else :
                    colo =  'g'
                    uh +=1
                    amplit_rela_UH.append(max_UH)
                    UH_abs.append(k)
                    UH_abs_j.append(w_ind_UH)
                    
                    
                # moy = np.mean(signal[0,:j],axis=1)
                # std = np.std(signal[:,:j],axis=1)
                # if signal[0,j]>moy[0]+std_tresh*std[0] or signal[1,j]>moy[1]+std_tresh*std[1]:
                #     resultat.append(j)
                #     ratio_UH = (signal[0,j]-moy[0])/std[0]
                #     ratio_BB = (signal[1,j]-moy[1])/std[1]
                #     harm_amp.append(temp[0,k,j])
                #     if signal[0,j]>moy[0]+std_tresh*std[0] and signal[1,j]>moy[1]+std_tresh*std[1]:
                #         both +=1
                #         amplit_rela_UH.append(ratio_UH)
                #         UH_abs.append(k)
                #         UH_abs_j.append(j)
                #         amplit_rela_BB.append(ratio_BB)
                #         BB_abs.append(k)
                #         BB_abs_j.append(j)
                        
                #     elif signal[1,j]>moy[1]+std_tresh*std[1] : 
                #         bb +=1
                #         amplit_rela_BB.append(ratio_BB)
                #         BB_abs.append(k)
                #         BB_abs_j.append(j)
                #     else :
                #         uh +=1
                #         amplit_rela_UH.append(ratio_UH)
                #         UH_abs.append(k)
                #         UH_abs_j.append(j) 
            multiplicateur = [1, UH_mult, BB_mult]
            if plot_true and (k==616 ):#or (max_UH>0 or max_BB>0)

                moy_UH = np.mean(signal[0,:still_wind])
                moy_BB = np.mean(signal[1,:still_wind])
                std_UH = np.std(signal[0,:still_wind])
                std_BB = np.std(signal[1,:still_wind])
                for i in range(n_plot):
                    plt.plot(time,temp_norm[i,k,:] *multiplicateur[i] , c=redblue[i],label=label[i], marker='.', markersize=10)
                
                #if max_UH>0:
                plt.fill_between([0,time[-1]], [moy_UH+std_tresh*std_UH,moy_UH+std_tresh*std_UH], [moy_UH-std_tresh*std_UH,moy_UH-std_tresh*std_UH],color=redblue[1], alpha=0.1)
                    #plt.plot([0,pression[-1]],[moy_UH+std_tresh*moy_UH,moy_UH+std_tresh*moy_UH],c=redblue[1])
                    #plt.plot([w_ind_UH,w_ind_UH],[0,1],c='g')
                #if max_BB>0:
                plt.fill_between([0,time[-1]], [moy_BB+std_tresh*std_BB,moy_BB+std_tresh*std_BB], [moy_BB-std_tresh*std_BB,moy_BB-std_tresh*std_BB],color=redblue[2], alpha=0.1)
                    #plt.plot([0,pression[-1]],[moy_BB+std_tresh*std_BB,moy_BB+std_tresh*std_BB],c=redblue[2])
                    #plt.plot([w_ind_BB,w_ind_BB],[0,1],c='r')
                plt.legend(fontsize=20)
                if (max_UH>0 or max_BB>0):
                    plt.title("Différentes composantes, tir avec événement inertiel",color="black",fontsize=35, fontweight = 'bold')
                else :
                    plt.title("Différentes composantes, tir sans événement inertiel",color="black",fontsize=35, fontweight = 'bold')
                plt.xlabel('Pression (kPa)',fontsize=30)
                plt.xlabel('Temps (ms)',fontsize=30)
                plt.ylabel('Amplitude normalisée',fontsize=30)
                #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')
                # Ymax = 1.01*np.amax(plot[:,:,1])
                # Ymin = 0.99*np.amin(plot[:,:,2])
                plt.ylim([0-val_std*1.1, 1.04]) 
                plt.xlim([0, time[-1]])
                plt.xticks(fontsize=22)
                plt.yticks(fontsize=22)
                plt.xlim([3,14])
                plt.savefig(chemin_lin+add+nom+'_{}_all.png'.format(k),bbox_inches='tight')
                plt.xlim([8,11])
                plt.savefig(chemin_lin+add+nom+'_{}_end.png'.format(k),bbox_inches='tight')
                
                #plt.ylim([0, np.max(temp_norm[:,k,:])])
                #plt.xlim([0, 40])    
                plt.grid(True)
                plt.tight_layout()
                plt.yscale('linear')
                plt.xlim([0,15])
                plt.savefig(chemin_lin+add+nom+'_{}.png'.format(k),bbox_inches='tight')
                plt.clf()
                    #break
        plt.scatter(UH_abs,amplit_rela_UH,c='green',marker = "x",s=100, linewidths = 2)
        plt.scatter(BB_abs,amplit_rela_BB,c='red',marker = "x",s=100, linewidths = 2)
        plt.title("Différents événements dépassant le seuil de bruit (pulses)", color='black',fontsize=30, fontweight = 'bold')
        plt.xlabel("Pulses",fontsize=20)
        plt.ylabel("ratio de STD supérieur à la moyenne du bruit", color="blue",fontsize=20) 
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('linear')
        plt.savefig(chemin_lin+"00"+nom+'recap_std_n_pulse.png',bbox_inches='tight')
        plt.clf()
        plt.scatter(UH_abs_j,amplit_rela_UH,c='green',marker = "x",s=100, linewidths = 2)
        plt.scatter(BB_abs_j,amplit_rela_BB,c='red',marker = "x",s=100, linewidths = 2)
        plt.title("Différents événements dépassant le seuil de bruit (windows)",color='black',fontsize=30, fontweight = 'bold')
        plt.xlabel("Fenêtre de déclenchement dans le pulse",fontsize=20)
        plt.ylabel("ratio de STD supérieur à la moyenne du bruit", color="blue",fontsize=20) 
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('linear')
        plt.savefig(chemin_lin+"00"+nom+'recap_std_window.png',bbox_inches='tight')
        plt.clf()
        plt.plot(harm_amp_BB,c='darkviolet', label = 'Dépassement BB')
        plt.plot(harm_amp_UH,c='darkcyan', label = 'Dépassement UH')
        plt.legend()
        plt.title("Valeurs des harmoniques au moment du depassement de seuil",color='black',fontsize=30, fontweight = 'bold')
        plt.ylabel("Valeur non normalisée de la composante harmonique", color="blue",fontsize=20) 
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('linear')
        plt.savefig(chemin_lin+"00"+nom+'harm_value.png',bbox_inches='tight')
        plt.clf()
        plt.scatter(UH_abs_j,UH_abs,c='green',marker = "x",s=100, linewidths = 2)
        plt.scatter(BB_abs_j,BB_abs,c='red',marker = "x",s=100, linewidths = 2)
        plt.title("Position des différents événements (dans le pulse et dans le traitement)",color='black',fontsize=30, fontweight = 'bold')
        plt.xlabel("Fenêtre de déclenchement dans le pulse",fontsize=20)
        plt.ylabel("Pulses de déclenchement", color="blue",fontsize=20) 
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('linear')
        plt.savefig(chemin_lin+"00"+nom+'recap_pulse_window.png',bbox_inches='tight')
        plt.close("all") 
        fig=plt.figure(figsize=(20,11))
        plt.subplot(1,2,1)
        plt.scatter(['UH' for j in amplit_rela_UH],amplit_rela_UH,c='green',marker = "o",s=100, linewidths = 2, alpha=0.15)
        plt.scatter(['BB' for j in amplit_rela_BB],amplit_rela_BB,c='red',marker = "o",s=100, linewidths = 2, alpha=0.15)
        plt.title("Différents événements durant le pulse (Linear)",color='black',fontsize=30, fontweight = 'bold')
        plt.xlabel("Valeur de ration de STD",fontsize=20)
        plt.ylabel("Type de déclenchement", color="blue",fontsize=20) 
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('linear')
        plt.subplot(1,2,1)
        plt.scatter(['UH' for j in amplit_rela_UH],amplit_rela_UH,c='green',marker = "o",s=100, linewidths = 2, alpha=0.15)
        plt.scatter(['BB' for j in amplit_rela_BB],amplit_rela_BB,c='red',marker = "o",s=100, linewidths = 2, alpha=0.15)
        plt.title("Différents événements durant le pulse (LOG)",color='black',fontsize=30, fontweight = 'bold')
        plt.xlabel("Valeur de ration de STD",fontsize=20)
        plt.ylabel("Type de déclenchement", color="blue",fontsize=20) 
        plt.grid(True)
        plt.tight_layout()
        plt.yscale('log')
        plt.savefig(chemin_lin+"00"+nom+'recap_pulse_ratio.png',bbox_inches='tight')
        plt.close("all")     
        print("déclenchement : UH {}, BB {}, UH&BB {}".format(uh,bb,both)) ,
        if len(amplit_rela_UH)>0:
            print("UH ratio de valeur de déclenchement : {} +- ".format(np.round(np.mean(amplit_rela_UH),decimals=2)),np.round(np.std(amplit_rela_UH),decimals=2))
        if len(amplit_rela_BB)>0:
            print("BB ratio de valeur de déclenchement : {} +- ".format(np.round(np.mean(amplit_rela_BB),decimals=2)),np.round(np.std(amplit_rela_BB),decimals=2))
        print("fenêtre de déclenchement moyen UH : {}".format(np.mean(UH_abs_j)))
        print("fenêtre de déclenchement median UH : {}".format(np.median(UH_abs_j))) 
        print("fenêtre de déclenchement moyen BB : {}".format(np.mean(BB_abs_j)))
        print("fenêtre de déclenchement median BB : {}".format(np.median(BB_abs_j)))
             
class experiment_mult():
    
    def __init__(self, fe, f0,start:int,end:int, size_window : int = 2048, n_harm : int = 4,delta_harm =50e3,delta_Uharm =50e3,peak_detection : bool =True, size_decal : int = 2048):
        #valeurs initiales de l'experience
        self.freq_ech = float(fe)
        self.freq_emiss = float(f0)
        self.n_harm=n_harm
        self.start=start
        self.end=end
        self.size_decal=size_decal
        self.freq_harm=np.array([self.freq_emiss*i for i in range(1,n_harm+3)])
        self.freq_Uharm=np.array([self.freq_emiss*(i-0.5) for i in range(1,n_harm+3)])
        self.peak_detection=peak_detection
        self.size_window=size_window
        self.frequencies= np.arange(self.size_window//2) * self.freq_ech/self.size_window
        self.delta_harm=delta_harm
        self.delta_Uharm=delta_Uharm
        self.exp=[]
        self.n_exp=0
        
    def creat_expe(self,freq_ech=None,start=None,end=None):
        if freq_ech==None:
            freq_ech=self.freq_ech
        if start==None:
            start=self.start
        if end==None:
            end=self.end
        self.exp.append(experiment(freq_ech,self.freq_emiss,start=start,end=end))
        self.n_exp+=1
        
    def add_pulse(self,data,exp,freq_ech=None,start=None,end=None):
        if freq_ech==None:
            freq_ech=self.freq_ech
        if start==None:
            start=self.start
        if end==None:
            end=self.end
        pulse=temp_sample(freq_ech,self.freq_emiss,data,size_window= self.size_window, n_harm= self.n_harm,delta_harm =self.delta_harm,delta_Uharm =self.delta_Uharm,peak_detection =self.peak_detection,start=start,end=end,size_decal=self.size_decal)
        self.exp[exp].pulses.append(pulse)
        self.exp[exp].n_pulse+=1
    
    def add_pulses(self,array_data,exp,spacer,freq_ech=None,start=None,end=None):
        if freq_ech==None:
            freq_ech=self.freq_ech
        if start==None:
            start=self.start
        if end==None:
            end=self.end
        for data in array_data:
            pulse=temp_sample(freq_ech,self.freq_emiss,data,spacer = spacer, size_window= self.size_window, n_harm= self.n_harm,delta_harm =self.delta_harm,delta_Uharm =self.delta_Uharm,peak_detection =self.peak_detection,start=start,end=end,size_decal=self.size_decal)
            self.exp[exp].pulses.append(pulse)
            self.exp[exp].n_pulse+=1
    
    def plot_indice_together(self,chemin):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##HARMONIQUES
        plot=np.zeros((3,self.n_exp,self.exp[0].n_pulse))
        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse):
                plot[0,i,j]=np.mean(self.exp[i].pulses[j].indice_harm)
                plot[1,i,j]=np.mean(self.exp[i].pulses[j].indice_Uharm )
                plot[2,i,j]=np.mean(self.exp[i].pulses[j].indice_BB)
        plot_legend=["Indice repr2sentant les harmoniques","Indice représentant les ultra-harmoniques","Indice representant le bruit de bande"]
        legend=["no Mbs","Mbs 0.25","Mbs 0.50","Mbs 0.75","Mbs 100","Mbs 100/2"]
        colors=['black','red','peru','forestgreen','dodgerblue','indigo']
        redblue=['b','g','r']
        # Stable cavitation dose
        fig=plt.figure(figsize=(20,11))
        for i in range(3):
            fig.clf()
            for j in range(self.n_exp):
                plt.plot(plot[i,j], c=colors[j],label=legend[j])
            plt.legend()
            plt.title(plot_legend[i],color=redblue[i],fontsize=30, fontweight = 'bold')
            plt.xlabel('Pulses',fontsize=25)
            plt.ylabel(plot_legend[i]+' [a.u.]', color=redblue[i],fontsize=25, fontweight = 'bold')
            if plot_legend[i]!="Ancienne SCD lineaire":
                plt.yscale('log')
            #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')
            Ymax = 1.01*np.amax(plot[i])
            Ymin = 0.99*np.amin(plot[i])
            plt.ylim([Ymin, Ymax])    
            plt.grid(True)
            plt.tight_layout()
            plt.yscale('log')
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.savefig(chemin+'\\'+plot_legend[i]+'_log.png',bbox_inches='tight')
            plt.yscale('linear')
            plt.savefig(chemin+'\\'+plot_legend[i]+'_linear.png',bbox_inches='tight')
        plt.close("all")  

    def extract_composante(self,nbit,fit,legend = ['data '+str(int(i))for i in range(1,11)]):
        ##creation des indices
        n_plot = 3
        n_harm_chosen = 3
        plot=np.zeros((self.n_exp,self.exp[0].n_pulse,n_plot,self.exp[0].pulses[0].n_window,n_harm_chosen))
        
        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse):
                plot[i,j,0,:] = np.transpose(20*np.log10(self.exp[i].pulses[j].indice_harm_w[1:n_harm_chosen+1]))
                plot[i,j,1,:] = np.transpose(20*np.log10(self.exp[i].pulses[j].indice_Uharm_w[1:n_harm_chosen+1]))
                plot[i,j,2,:] = np.transpose(20*np.log10(self.exp[i].pulses[j].indice_BB_sliced_w[1:n_harm_chosen+1]))
        #je souhaiterai avant de renvoyer le resultat
        #inverser les axes 2 et 3 pour plutot avoir dans plot : (exp, pulses, fenetres, indice, composantes) au lieu de (exp, pulses, indice, fenetres, composantes)
        plot = np.swapaxes(plot,2,3)
        return plot

    def plot_windowed(self,chemin,nbit,fit,mini_value=0,maxi_value=100,n=1,legend = ['data '+str(int(i))for i in range(1,11)], ramp = False):
        y_axis = 'Pression (kPa)'
        if ramp:
            y_axis = 'Différents tirs'
            
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##recuperation des indices
        n_plot = 4
        plot=np.zeros((n_plot,self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        
        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[1:-1] for a in range(n)] ,axis=0),axis=0)
                plot[0,i,j,:]=20*np.log10(val) #20*np.log10
                
                val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_w[1:4] for a in range(n)] ,axis=0),axis=0)
                plot[1,i,j,:]=20*np.log10(val) #20*np.log10
                
                val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_norm_div_w[1:4] for a in range(n)] ,axis=0),axis=0)
                plot[2,i,j,:]=20*np.log10(val) #20*np.log10
                
                val = np.mean([self.exp[i].pulses[j*n+a].indice_BB_w for a in range(n)] ,axis=0)
                plot[3,i,j,:]=20*np.log10(val) #20*np.log10

        fig=plt.figure(figsize=(20,11))
        amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
        pression=[fit[0]*j+fit[1] for j in amp]
        self.n_window = self.exp[0].pulses[0].n_window
        nom_str = ["\\H_n{}_perLow{}_perHigh{}_ramp.png".format(n,mini_value,maxi_value),"\\UH_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\UH_norm_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\BB_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value)]
        title_str = ["Harmonic components","Ultra-harmonic components","Composante ultra-harmonique émergente","Inertial components"]
        title_str = ["Indices harmoniques","Indices ultra-harmoniques","Indices ultra-harmoniques émergents","Indices inertiels"]
        for j in range(n_plot):
            plt.clf()
            fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
            fig.suptitle(title_str[j],fontsize=25, fontweight = 'bold',y=0.94)
            mini=np.percentile(plot[j],mini_value)
            maxi=np.percentile(plot[j],maxi_value)
            for i,ax in enumerate(axes.flat):
                im = ax.imshow(plot[j,i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,1*np.shape(plot[i])[1]*self.size_decal/self.freq_ech*1000.,((nbit[1])*fit[0]+fit[1]),((nbit[0])*fit[0]+fit[1])], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
                ax.title.set_text(legend[i].upper())
                ax.title.set_fontweight('bold')
                ax.title.set_fontsize(18)
                ax.set_xlabel('Temps'.upper()+' (ms)',fontsize=20)
                if i==0:
                    ax.set_ylabel(y_axis,fontsize=30)
                ax.tick_params(axis = 'both',which = 'both',  labelsize=16)
                
            #je veux grossir la taille de la police des ticks de la colorbar
            fig.colorbar(im,ax=axes.ravel().tolist())
            #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
            plt.savefig(chemin+nom_str[j],bbox_inches='tight')
        plt.close("all")
        
        
    def plot_windowed_RAMP(self,chemin,nbit,fit,mini_value=0,maxi_value=100,n=1,legend = ['data '+str(int(i))for i in range(1,11)], ramp = False):
        y_axis = 'Pressure (KPa)'
        if ramp:
            y_axis = 'Different pulses (ramping pressure)'
            
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##recuperation des indices
        n_fen = self.exp[0].pulses[0].n_window//2
        n_plot = 4
        plot=np.zeros((n_plot,self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        plot_ramp = np.zeros((n_plot,self.n_exp,self.exp[0].n_pulse//n,n_fen))
        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[1:-1] for a in range(n)] ,axis=0),axis=0)
                plot[0,i,j,:]=20*np.log10(val) #20*np.log10
                
                val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_w[1:4] for a in range(n)] ,axis=0),axis=0)
                plot[1,i,j,:]=20*np.log10(val) #20*np.log10aazazeazaazeazeaazaa
                
                val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_norm_div_w[1:4] for a in range(n)] ,axis=0),axis=0)
                plot[2,i,j,:]=20*np.log10(val) #20*np.log10
                
                val = np.mean([self.exp[i].pulses[j*n+a].indice_BB_w for a in range(n)] ,axis=0)
                plot[3,i,j,:]=20*np.log10(val) #20*np.log10
        for i in range(1,n_fen+1):
            plot_ramp[:,:,:,i-1] = plot[:,:,:,n_fen-i]-plot[:,:,:,n_fen+i-1]
            
            
            
        fig=plt.figure(figsize=(20,11))
        amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
        pression=[fit[0]*j+fit[1] for j in amp]
        self.n_window = self.exp[0].pulses[0].n_window
        nom_str = ["\\H_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\UH_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\UH_norm_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value),"\\BB_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value)]
        title_str = ["Harmonic components","Ultra-harmonic components","Composante ultra-harmonique émergente","Inertial components"]
        for j in range(n_plot):
            plt.clf()
            fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
            fig.suptitle(title_str[j],fontsize=18, fontweight = 'bold',y=0.94)
            mini=np.percentile(plot_ramp[j],mini_value)
            maxi=np.percentile(plot_ramp[j],maxi_value)
            for i,ax in enumerate(axes.flat):
                im = ax.imshow(plot_ramp[j,i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[i])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
                ax.title.set_text(legend[i].upper())
                ax.title.set_fontweight('bold')
                ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
                if i==0:
                    ax.set_ylabel(y_axis)
            
            fig.colorbar(im,ax=axes.ravel().tolist())
            #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
            plt.savefig(chemin+nom_str[j],bbox_inches='tight')
        plt.close("all")
        
    # def plot_UH_windowed(self,chemin,nbit,fit,mini_value=0,maxi_value=100,n=1,legend = ['data '+str(int(i))for i in range(1,11)], ramp = False):
    #     y_axis = 'Pressure (KPa)'
    #     if ramp:
    #         y_axis = 'Different pulses (ramping pressure)'
            
    #     if not os.path.isdir(chemin): # check if folder exists, otherwise create it
    #         os.mkdir(chemin)
    #     ##recuperation des indices
    #     plot=np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        
    #     for i in range(self.n_exp):
    #         for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
    #             val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_w[1:4] for a in range(n)] ,axis=0),axis=0)
    #             plot[i,j,:]=20*np.log10(val) #20*np.log10
                
    #     self.ratio=plot
    #     fig=plt.figure(figsize=(20,11))
    #     #fit = np.array([4.66745471*2, 5.80567673])
    #     amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
    #     pression=[fit[0]*j+fit[1] for j in amp]
    #     self.n_window = self.exp[0].pulses[0].n_window

    #     nom = "\\UH_windowed_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value)
    #     fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
    #     fig.suptitle("Composante ultra-harmonique",fontsize=18, fontweight = 'bold',y=0.94)
    #     mini=np.percentile(plot,mini_value)
    #     maxi=np.percentile(plot,maxi_value)
    #     for i,ax in enumerate(axes.flat):
    #         im = ax.imshow(plot[i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[i])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
    #         ax.title.set_text(legend[i].upper())
    #         ax.title.set_fontweight('bold')
    #         ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
    #         if i==0:
    #             ax.set_ylabel(y_axis)
        
    #     fig.colorbar(im,ax=axes.ravel().tolist())
    #     #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
    #     plt.savefig(chemin+nom,bbox_inches='tight')
    #     plt.clf()
    #     plt.close("all")

    # def plot_UH_norm_windowed(self,chemin,nbit,fit,mini_value=0,maxi_value=100,n=1,legend = ['data '+str(int(i))for i in range(1,11)],ramp = False):
    #     y_axis = 'Pressure (KPa)'
    #     if ramp:
    #         y_axis = 'Different pulses (ramping pressure)'
    #     if not os.path.isdir(chemin): # check if folder exists, otherwise create it
    #         os.mkdir(chemin)
    #     ##recuperation des indices
    #     plot=np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
    #     for i in range(self.n_exp):
    #         for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
    #             val = np.mean(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_norm_div_w[1:4] for a in range(n)] ,axis=0),axis=0)
    #             plot[i,j,:]=20*np.log10(val) #20*np.log10

    #     self.ratio=plot
    #     fig=plt.figure(figsize=(20,11))
    #     #fit = np.array([4.66745471*2, 5.80567673])
    #     amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
    #     pression=[fit[0]*j+fit[1] for j in amp]
    #     self.n_window = self.exp[0].pulses[0].n_window

    #     nom = "\\UH_norm_windowed_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value)
    #     fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
    #     fig.suptitle("Composante ultra-harmonique normalisée",fontsize=18, fontweight = 'bold',y=0.94)
    #     mini=np.percentile(plot,mini_value)
    #     maxi=np.percentile(plot,maxi_value)
    #     for i,ax in enumerate(axes.flat):
    #         im = ax.imshow(plot[i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[i])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
    #         ax.title.set_text(legend[i].upper())
    #         ax.title.set_fontweight('bold')
    #         ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
    #         if i==0:
    #             ax.set_ylabel(y_axis)
        
    #     fig.colorbar(im,ax=axes.ravel().tolist())
    #     #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
    #     plt.savefig(chemin+nom,bbox_inches='tight')
    #     plt.clf()
    #     plt.close("all")
        
        
        
    # def plot_H_windowed(self,chemin,nbit,fit,mini_value=0,maxi_value=100,n=1,legend = ['data '+str(int(i))for i in range(1,11)],ramp = False):
    #     y_axis = 'Pressure (KPa)'
    #     if ramp:
    #         y_axis = 'Different pulses (ramping pressure)'
    #     if not os.path.isdir(chemin): # check if folder exists, otherwise create it
    #         os.mkdir(chemin)
    #     ##recuperation des indices
    #     plot=np.zeros((self.n_harm , self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        
    #     for i in range(self.n_exp):
    #         for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
    #             for k in range(self.n_harm):
    #                 val = np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[k] for a in range(n)] ,axis=0)
    #                 plot[k,i,j,:]=20*np.log10(val) #20*np.log10

    #     self.ratio=plot
    #     fig=plt.figure(figsize=(20,11))
        
    #     #fit = np.array([4.66745471*2, 5.80567673])
    #     amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
    #     pression=[fit[0]*j+fit[1] for j in amp]
    #     self.n_window = self.exp[0].pulses[0].n_window

        
    #     for l in range(self.n_harm):
    #         nom = "\\H_windowed_n{}_f{}_perLow{}_perHigh{}.png".format(n,l,mini_value,maxi_value)
    #         fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
            
    #         fig.suptitle("Composante harmonique f{} normalisée".format(l),fontsize=18, fontweight = 'bold',y=0.94)
    #         mini=np.percentile(plot[l],mini_value)
    #         maxi=np.percentile(plot[l],maxi_value)
    #         for i,ax in enumerate(axes.flat):
                
    #             im = ax.imshow(plot[l,i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[l,i])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
    #             ax.title.set_text(legend[i].upper())
    #             ax.title.set_fontweight('bold')
    #             ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
    #             if i==0:
    #                 ax.set_ylabel(y_axis)
            
    #         fig.colorbar(im,ax=axes.ravel().tolist())
    #         #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
    #         plt.savefig(chemin+nom,bbox_inches='tight')
    #         plt.clf()
    #     plt.close("all")
        
    # def plot_BB_windowed(self,chemin,nbit,fit,mini_value=0,maxi_value=100,n=1,legend = ['data '+str(int(i))for i in range(1,11)],ramp = False):
    #     y_axis = 'Pressure (KPa)'
    #     if ramp:
    #         y_axis = 'Different pulses (ramping pressure)'
    #     if not os.path.isdir(chemin): # check if folder exists, otherwise create it
    #         os.mkdir(chemin)
    #     ##recuperation des indices
    #     plot=np.zeros(( self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        
    #     for i in range(self.n_exp):
    #         for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
    #             val = np.mean([self.exp[i].pulses[j*n+a].indice_BB_w for a in range(n)] ,axis=0)
    #             plot[i,j,:]=20*np.log10(val) #20*np.log10

    #     self.ratio=plot
    #     fig=plt.figure(figsize=(20,11))
        
    #     #fit = np.array([4.66745471*2, 5.80567673])
    #     amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
    #     pression=[fit[0]*j+fit[1] for j in amp]
    #     self.n_window = self.exp[0].pulses[0].n_window

        
    #     nom = "\\BB_windowed_n{}_perLow{}_perHigh{}.png".format(n,mini_value,maxi_value)

    #     fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
    #     fig.suptitle("Composante inertielle",fontsize=18, fontweight = 'bold',y=0.94)
    #     mini=np.percentile(plot,mini_value)
    #     maxi=np.percentile(plot,maxi_value)
    #     for i,ax in enumerate(axes.flat):
    #         im = ax.imshow(plot[i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[i])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
    #         ax.title.set_text(legend[i].upper())
    #         ax.title.set_fontweight('bold')
    #         ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
    #         if i==0:
    #             ax.set_ylabel(y_axis)
        
    #     fig.colorbar(im,ax=axes.ravel().tolist())
    #     #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
    #     plt.savefig(chemin+nom,bbox_inches='tight')
    #     plt.clf()
    #     plt.close("all")
        
    def plot_indice_together_grp(self,nom,chemin,nbit,n=1,legend = ['data '+str(int(i))for i in range(1,11)],fit = np.array([4.66745471*2, 5.80567673]),beta = 7):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        chemin_lin=chemin
        ##HARMONIQUES
        n_plot = 5
        plot=np.zeros((n_plot,self.n_exp,self.exp[0].n_pulse//n,3))
        #fenetre regardees :
        n0,n1, nf = 1,9,2
        temp = np.zeros((2,n_plot,n))
        #temp_bis = np.zeros((2,n_plot,n,len(harm_num)))
        for i in range(self.n_exp):
            temp = np.zeros((2,n_plot,n))
            temp_bis = np.zeros((2,n_plot,n,nf))
            for j in range(self.exp[0].n_pulse//n):
                temp[0,0] =[np.mean(self.exp[i].pulses[j*n+a].indice_harm[1:3]) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
                temp[0,1] =[np.mean(self.exp[i].pulses[j*n+a].indice_Uharm[1:4]) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
                temp[0,2] =[np.mean(self.exp[i].pulses[j*n+a].indice_BB) for a in range(n)] #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
                temp[0,3] =[np.mean(self.exp[i].pulses[j*n+a].indice_Uharm_norm_div[1:4]) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
                temp[0,4] = temp[0,1] - beta * temp[0,2]
                for k in range(n_plot):
                    plot[k,i,j,0]=np.mean(temp[0,k])
                    plot[k,i,j,1]=plot[k,i,j,0]+np.std(temp[0,k])
                    plot[k,i,j,2]=plot[k,i,j,0]-np.std(temp[0,k])
                    
        plot_legend=["Harmonic indexes","Ultraharmonic indexes","Broadband noise indexes","Indices représentant les ultra-harmoniques-BB","Indice representant UH-{}BB".format(beta),"Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit large bande","Indices représentant les ultra-harmoniques émergents","Indice representant UH-{}BB".format(beta),"Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        nom_img=["harm","U_harm","BB","U_harm_norm","UH_BB","ratio_"+nom+"_all","ratio_"+nom+"_start","ratio_"+nom+"_end"]
        colors=['black','red','peru','forestgreen','dodgerblue','gold']
        colors=['blue','green','red','forestgreen','dodgerblue','gold']
        colors_harm = ['black', 'mediumturquoise', 'dodgerblue','navy']
        colors_Uharm = ['black', 'olivedrab', 'limegreen','darkgreen']
        colors_BB = ['black', 'lightcoral', 'red','darkred']
        colors_grad = []
        for j in range(n_plot):
            colors_grad.append([])
            for k in range(self.n_exp):
                if k ==0 :
                    colors_grad[j].append('black')
                else :
                    colors_grad[j].append(colors[j])
        colors_grad[0] = colors_harm
        colors_grad[1] = colors_Uharm
        colors_grad[2] = colors_BB
        line= ['-','-','-','-','-','--','-']
        #fit = np.array([4.66745471*2, 5.80567673])
        if fit==False : 
            pression = np.arange(self.exp[0].n_pulse//n)
            fit = np.array([1, 0])
        else :
            amp=np.arange(nbit[0],nbit[1])
            pression=[fit[0]*j+fit[1] for j in amp]
        redblue=['b','g','r','teal','m','black','black','black']
        # Stable cavitation dose
        fig=plt.figure(figsize=(20,11))
        for i in range(n_plot):
            fig.clf()
            for j in range(self.n_exp):
                plt.plot(pression,plot[i,j,:,0], c=colors_grad[i][j],label=legend[j],alpha=1)
                plt.fill_between(pression,plot[i,j,:,1],plot[i,j,:,2], color=colors_grad[i][j],alpha=0.20)
            plt.legend(fontsize=24, loc = 'upper left')
            plt.title(plot_legend[i],color=redblue[i],fontsize=30, fontweight = 'bold')
            #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')
            Ymax = 1.01*np.amax(plot[i,:,:,1])
            Ymin = 0.99*np.amin(plot[i,:,:,2])
            zomm_min = 50
            zoom_max = 900
            Ymax_zoom = 1.01*np.amax(plot[i,:,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),1])
            Ymin_zoom = 0.99*np.amin(plot[i,:,p_to_b(fit,zomm_min):p_to_b(fit,zoom_max),2])
            #plt.ylim([Ymin_zoom, Ymax_zoom])        
            plt.grid(True)
            plt.tight_layout()
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.xlim([0, np.max(pression)]) 
            plt.xlabel('Pression (kPa)',fontsize=25)
            plt.ylabel('Amplitude [a.u.]', color=redblue[i],fontsize=25)
            plt.yscale('linear')
            plt.savefig(chemin_lin+nom_img[i]+'_linear.png',bbox_inches='tight')
            plt.xlim([zomm_min, zoom_max])
            plt.ylim([Ymin_zoom, Ymax_zoom])
            plt.savefig(chemin_lin+nom_img[i]+'_linear_ZOOM.png',bbox_inches='tight')            
        plt.close("all")  

    def plot_indice_RAMP(self,nom,title,chemin,pression,legend = ['data '+str(int(i))for i in range(1,11)]):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        chemin_lin=chemin
        n_window = self.exp[0].pulses[0].n_window
        ##HARMONIQUES
        n_plot = 4
        plot=np.zeros((n_plot,self.n_exp,n_window,3))
        for i in range(self.n_exp):
            temp = np.zeros((n_plot,self.exp[i].n_pulse,n_window))
            for j in range(self.exp[0].n_pulse):
                temp[0,j] =np.mean(self.exp[i].pulses[j].indice_harm_w[1:-1],axis=0)#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
                temp[1,j] =np.mean(self.exp[i].pulses[j].indice_Uharm_w[1:4],axis=0) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
                temp[2,j] =self.exp[i].pulses[j].indice_BB_w  #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
                temp[3,j] =np.mean(self.exp[i].pulses[j].indice_Uharm_norm_div_w[1:4],axis=0) #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+

            for k in range(n_plot):
                plot[k,i,:,0]=np.mean(temp[k], axis=0)
                plot[k,i,:,1]=plot[k,i,:,0]+np.std(temp[k], axis=0)
                plot[k,i,:,2]=plot[k,i,:,0]-np.std(temp[k], axis=0)
                
        plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit large bande","Indices représentant les ultra-harmoniques émergents","Indice representant UH-{}BB".format(5),"Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        nom_img=["harm","U_harm","BB","U_harm_norm","UH_BB","ratio_"+nom+"_all","ratio_"+nom+"_start","ratio_"+nom+"_end"]
        colors=['blue','green','red','forestgreen','dodgerblue','gold']
        colors_harm = ['black', 'mediumturquoise', 'dodgerblue','navy']
        colors_Uharm = ['black', 'olivedrab', 'limegreen','darkgreen']
        colors_BB = ['black', 'lightcoral', 'red','darkred']
        colors_grad = []
        for j in range(n_plot):
            colors_grad.append([])
            for k in range(self.n_exp):
                if k ==0 :
                    colors_grad[j].append('black')
                else :
                    colors_grad[j].append(colors[j])
        colors_grad[0] = colors_harm
        colors_grad[1] = colors_Uharm
        colors_grad[2] = colors_BB
        colors_grad[3] = colors_Uharm
        line= ['-','-','-','-','-','--','-']
        #fit = np.array([4.66745471*2, 5.80567673])
        redblue=['b','g','r','teal','m','black','black','black']
        # Stable cavitation dose
        fig=plt.figure(figsize=(20,11))
        for i in range(n_plot):
            fig.clf()
            for j in range(self.n_exp):
                plt.plot(pression,plot[i,j,:,0], c=colors_grad[i][j],label=title[j])
                plt.fill_between(pression,plot[i,j,:,1],plot[i,j,:,2], color=colors_grad[i][j],alpha=0.2)
            plt.legend(fontsize=24, loc = 'upper left')
            plt.title(plot_legend[i],color=redblue[i],fontsize=30, fontweight = 'bold')
            plt.xlabel('Pression (kPa)',fontsize=25)
            plt.ylabel('Amplitude [a.u.]', color=redblue[i],fontsize=25)
            #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')
            Ymax = 1.01*np.amax(plot[i,:,:,1])
            Ymin = 0.99*np.amin(plot[i,:,:,2])
            plt.ylim([Ymin, Ymax])         
            plt.grid(True)
            plt.tight_layout()
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.xlim([0, np.max(pression)]) 
            plt.yscale('linear')
            plt.savefig(chemin_lin+nom_img[i]+'_linear_ramp.png',bbox_inches='tight')           
        plt.close("all")  

    def plot_ratio_bulles_map(self,harm_num,harm_den,nom,nbit,n=1,legend = ['data '+str(int(i))for i in range(1,11)],):

        ##recuperation des indices
        plot=np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        
        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                num = np.zeros(self.exp[0].pulses[0].n_window)
                for elt in harm_num:
                    num += np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[elt] for a in range(n)] ,axis=0)
                den = np.zeros(self.exp[0].pulses[0].n_window)
                for elt in harm_den:
                    den +=  np.mean([ self.exp[i].pulses[j*n+a].indice_harm_w[elt] for a in range(n)] ,axis=0)
                ratio = num / den
                plot[i,j,:]=20*np.log10(ratio) #20*np.log10

        return plot

    def recup_harm_map(self,nbit,n=1,rep=30):
        n_harm=self.exp[0].pulses[0].n_harm
        ##recuperation des indices
        data = []
        c_MB = np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        data = np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window,9))#ayo n_harm*2+2
        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                val=np.zeros((self.exp[0].pulses[0].n_window,9))#ayo n_harm*2+2
                for k in range(3):#ayo
                    #print((i,j,k))
                    val[:,k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[k+1] for a in range(n)] ,axis=0))
                for k in range(2):#ayo
                    val[:,3+k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_w[k] for a in range(n)] ,axis=0))
                mb=[]
                for po,elt in enumerate(val):
                    elt[-1]=j//rep + nbit[0]
                    mb.append(i)
                data[i,j]=val
                c_MB[i,j]=mb
        return data,c_MB
    
    # def recup_harm_map(self,nbit,n=1,rep=30):
    #     n_harm=self.exp[0].pulses[0].n_harm
    #     ##recuperation des indices
    #     plot=np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window,n_harm*2+1))
    #     data = []
    #     c_MB = []
    #     for i in range(self.n_exp):
    #         for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                
    #             val=np.zeros((self.exp[0].pulses[0].n_window,n_harm*2+2))
    #             for k in range(n_harm):
    #                 val[:,k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[k+1] for a in range(n)] ,axis=0))
    #             for k in range(n_harm+1):
    #                 val[:,n_harm+k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_w[k] for a in range(n)] ,axis=0))
    #             for elt in val:
    #                 elt[-1]=j//rep + nbit[0]
    #                 data.append(elt)
    #                 c_MB.append(i)
    #     return data,c_MB
    
    
    def creat_concentration_map(self,model,legend,chemin,nbit,n=1,rep=30,ante=1,vmin=0,vmax=4,save = False):
        n_harm=self.exp[0].pulses[0].n_harm
        ##recuperation des indices
        plot=np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window-ante+1))
        fit = np.array([4.66745471*2, 5.80567673])

        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                
                val=np.zeros((self.exp[0].pulses[0].n_window,9))#ayo n_harm*2+2
                for k in range(3):#ayo
                    val[:,k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[k+1] for a in range(n)] ,axis=0))
                for k in range(5):#ayo
                    val[:,3+k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_w[k] for a in range(n)] ,axis=0))
                data=[]
                for po,elt in enumerate(val):
                    elt[-1]=j//rep + nbit[0]
                    fill=np.concatenate([val[po-ij] for ij in range(ante)])
                    if po>ante-2:
                        data.append(fill)                       
                #print(np.shape(np.array(data)))
                d= model.predict(data)
                plot[i,j,:]=d
                    
        fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
        for i,ax in enumerate(axes.flat):
            im = ax.imshow(plot[i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[i])[1]*self.size_decal/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]],vmin=vmin,vmax=vmax)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
            ax.title.set_text(legend[i].upper() + ' logarithmic')
            ax.title.set_fontweight('bold')
            ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
            if i==0:
                ax.set_ylabel('Pressure (KPa)')
        fig.colorbar(im,ax=axes.ravel().tolist())
        nom = "\\c_map_all.png"
        plt.savefig(chemin+nom,bbox_inches='tight')
        plt.clf()
                

        if save :        
            np.save(chemin+"\\map_complete.npy",plot)
        return plot
    
    
    # def creat_concentration_map(self,model,legend,chemin,nbit,n=1,rep=30):
    #     n_harm=self.exp[0].pulses[0].n_harm
    #     ##recuperation des indices
    #     plot=np.zeros((self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
    #     fit=np.array([29.95777559,  6.81690626])

    #     for i in range(self.n_exp):
    #         for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                
    #             val=np.zeros((self.exp[0].pulses[0].n_window,n_harm*2+2))
    #             for k in range(n_harm):
    #                 val[:,k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[k+1] for a in range(n)] ,axis=0))
    #             for k in range(n_harm+1):
    #                 val[:,n_harm+k] =  20*np.log10(np.mean([self.exp[i].pulses[j*n+a].indice_Uharm_w[k] for a in range(n)] ,axis=0))
    #             for elt in val:
    #                 elt[-1]=j//rep + nbit[0]
    #             plot[i,j,:]=model.predict(val)
                
    #     fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
    #     for i,ax in enumerate(axes.flat):
    #         im = ax.imshow(plot[i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[i])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]])  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
    #         ax.title.set_text(legend[i].upper() + ' logarithmic')
    #         ax.title.set_fontweight('bold')
    #         ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
    #         if i==0:
    #             ax.set_ylabel('Pressure (KPa)')
    #     nom = "\\concentration_map_all.png"
    #     fig.colorbar(im,ax=axes.ravel().tolist())
    #     #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
    #     plt.savefig(chemin+nom,bbox_inches='tight')
    #     plt.clf()
    #     return plot
    
    
    def plot_ratio_bulles_windowed(self,harm_num,harm_den,nom,chemin,nbit,n=1,legend = ['data '+str(int(i))for i in range(1,11)],):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##recuperation des indices
        plot=np.zeros((3,self.n_exp,self.exp[0].n_pulse//n,self.exp[0].pulses[0].n_window))
        
        for i in range(self.n_exp):
            for j in range(self.exp[0].n_pulse//n):#self.exp[i].pulses[j].indice_harm_w[1]+
                num = np.zeros(self.exp[0].pulses[0].n_window)
                for elt in harm_num:
                    num += np.mean([self.exp[i].pulses[j*n+a].indice_harm_w[elt] for a in range(n)] ,axis=0)
                den = np.zeros(self.exp[0].pulses[0].n_window)
                for elt in harm_den:
                    den +=  np.mean([ self.exp[i].pulses[j*n+a].indice_harm_w[elt] for a in range(n)] ,axis=0)
                ratio = num / den
                plot[0,i,j,:]=20*np.log10(ratio) #20*np.log10

        self.ratio=plot
        fig=plt.figure(figsize=(20,11))
        fit = np.array([4.66745471*2, 5.80567673])
        amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
        pression=[fit[0]*j+fit[1] for j in amp]
        self.n_window = self.exp[0].pulses[0].n_window

        mini_value = [10,20,30]
        maxi_value = [99,99,99]
        nom = ["\\ratio_"+nom+"_perLow{}_perHigh{}.png".format(mini_value[k],maxi_value[k]) for k in range(3)]
        for k in range(3):
            fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
            mini=np.percentile(plot[0],mini_value[k])
            maxi=np.percentile(plot[0],maxi_value[k])
            for i,ax in enumerate(axes.flat):
                im = ax.imshow(plot[0,i,:,:],aspect='auto',interpolation='none',cmap='turbo',extent=[0,np.shape(plot[k,i])[1]*self.size_window/self.freq_ech*1000.,(nbit[1])*fit[0]+fit[1],(nbit[0])*fit[0]+fit[1]], vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
                ax.title.set_text(legend[i].upper() + ' logarithmic')
                ax.title.set_fontweight('bold')
                ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
                if i==0:
                    ax.set_ylabel('Pressure (KPa)')
            
            fig.colorbar(im,ax=axes.ravel().tolist())
            #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
            plt.savefig(chemin+nom[k],bbox_inches='tight')
            plt.clf()
        plt.close("all")      
        
    
        
    def plot_indice_together_grp_vivo(self,harm_num,harm_den,nom,chemin,nbit,n=1,legend = ['data '+str(int(i))for i in range(1,11)],):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        chemin_log=chemin+'\\log\\'
        chemin_lin=chemin#+'\\linear\\'
# =============================================================================
#         if not os.path.isdir(chemin_log): # check if folder exists, otherwise create it
#             os.mkdir(chemin_log)
#         if not os.path.isdir(chemin_lin): # check if folder exists, otherwise create it
#             os.mkdir(chemin_lin)
# =============================================================================
        ##HARMONIQUES
        n_pulse_max=0
        for i in range(self.n_exp):
            if self.exp[i].n_pulse//n>n_pulse_max:
                n_pulse_max=self.exp[i].n_pulse
        plot=np.zeros((6,self.n_exp,n_pulse_max,3))
        
        #fenetre regardees :
        n0,n1, nf = 2,50,3
        n_plot = 6
        temp = np.zeros((2,n_plot,n))
        #temp_bis = np.zeros((2,n_plot,n,len(harm_num)))
        for i in range(self.n_exp):
            posi=0
            
            temp = np.zeros((2,n_plot,n))
            temp_bis = np.zeros((2,n_plot,n,nf))
            for j in range(self.exp[i].n_pulse//n):
                temp[0,0] =[np.mean(self.exp[i].pulses[j*n+a].indice_harm) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm[1])+
                temp[0,1] =[np.mean(self.exp[i].pulses[j*n+a].indice_Uharm) for a in range(n)]#np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3])+
                temp[0,2] =[np.mean(self.exp[i].pulses[j*n+a].indice_BB) for a in range(n)] #np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11])+
                for elt in harm_num:
                    temp[0,3] +=[self.exp[i].pulses[j*n+a].indice_harm[elt] for a in range(n)]
                for elt in harm_den:
                    temp[1,3] +=[self.exp[i].pulses[j*n+a].indice_harm[elt] for a in range(n)]
                for elt in harm_num:
                    #print(np.shape(self.exp[i].pulses[j*n+3].indice_harm_w[elt,n0:n0 + nf]))
                    temp_bis[0,4] +=[self.exp[i].pulses[j*n+a].indice_harm_w[elt,n0:n0 + nf] for a in range(n)]
                for elt in harm_den:
                    temp_bis[1,4] +=[self.exp[i].pulses[j*n+a].indice_harm_w[elt,n0:n0 + nf] for a in range(n)]
                for elt in harm_num:
                    temp_bis[0,5] +=[self.exp[i].pulses[j*n+a].indice_harm_w[elt,n1:n1 + nf] for a in range(n)]
                for elt in harm_den:
                    temp_bis[1,5] +=[self.exp[i].pulses[j*n+a].indice_harm_w[elt,n1:n1 + nf] for a in range(n)]
# =============================================================================
#                 temp[3] =[(np.mean(self.exp[i].pulses[j*n+a].indice_harm[2])+np.mean(self.exp[i].pulses[j*n+a].indice_harm[1]))/np.mean(self.exp[i].pulses[j*n+a].fondamental_w) for a in range(n)]
#                 temp[4] =[(np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[2,1:3])+np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,1:3]))/np.mean(self.exp[i].pulses[j*n+a].fondamental_w[1:3]) for a in range(n)]
#                 temp[5] =[(np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[2,9:11])+np.mean(self.exp[i].pulses[j*n+a].indice_harm_w[1,9:11]))/np.mean(self.exp[i].pulses[j*n+a].fondamental_w[1:3]) for a in range(n)]
# =============================================================================
                for k in range(3,4):
                    plot[k,i,j,0]=np.mean(temp[0,k]/temp[1,k])
                    plot[k,i,j,1]=plot[k,i,j,0]+np.std(temp[0,k]/temp[1,k])
                    plot[k,i,j,2]=plot[k,i,j,0]-np.std(temp[0,k]/temp[1,k])
                for k in range(4,6):
                    plot[k,i,j,0]=np.mean(temp_bis[0,k]/temp_bis[1,k])
                    plot[k,i,j,1]=plot[k,i,j,0]+np.std(np.mean(temp_bis[0,k]/temp_bis[1,k],axis=1))
                    plot[k,i,j,2]=plot[k,i,j,0]-np.std(np.mean(temp_bis[0,k]/temp_bis[1,k],axis=1))
        plot_legend=["Indices représentant les harmoniques","Indices représentant les ultra-harmoniques","Indices représentant le bruit de bande","Indices représentant le ratio de bulles "+nom+" sur tout le pulse","Indices représentant le ratio de bulles "+nom+" en début de pulse","Indices représentant le ratio de bulles "+nom+" en fin de pulse"]
        nom_img=["harm","U_harm","BB","ratio_"+nom+"_all","ratio_"+nom+"_start","ratio_"+nom+"_end"]
        colors=['black','red','peru','forestgreen','dodgerblue','gold']
        line= ['-','-','-','-','-','--','-']
        fit = np.array([4.66745471*2, 5.80567673])
        amp=np.arange(nbit[0],nbit[1])
        pression=[fit[0]*j+fit[1] for j in amp]
        redblue=['b','g','r','black','black','black']
        # Stable cavitation dose
        
        for i in range(3,6):
            
            fig, ax1 = plt.subplots(figsize=(20,11))
            color = 'tab:red'
            j=1
            ax1.set_xlabel('Pulses shot',fontsize=20)
            ax1.set_ylabel('Magnitude [a.u.]', color = color)  
            ax1.plot(plot[i,j,:self.exp[j].n_pulse//n,0], c=color,label=legend[j])
            ax1.tick_params(axis ='y', labelcolor = color)  
            j=0
            ax2 = ax1.twinx()  
            color = 'tab:green'
            ax2.set_ylabel('Magnitude [a.u.]', color = color)  
            ax2.plot(plot[i,j,:self.exp[j].n_pulse//n,0], c=color,label=legend[j])
            ax2.tick_params(axis ='y', labelcolor = color)
            plt.legend(fontsize=20)
            plt.title(plot_legend[i],color=redblue[i],fontsize=30, fontweight = 'bold')
            #plt.xlabel('Pression (KPa)',fontsize=20)
            #plt.ylabel('Magnitude [a.u.]', color=redblue[i],fontsize=20)
            #plt.plot([t[0], t[-1]], [1.0, 1.0],  ls='--',linewidth=2, c=redblue[i]) # plt.plot((x1, x2), (y1, y2), 'k-')
            Ymax = 1.01*np.amax(plot[i,:,:,1])
            Ymin = 0.99*np.amin(plot[i,:,:,2])
            #plt.ylim([Ymin, Ymax])    
            ax1.grid(True)
            plt.tight_layout()
            plt.savefig(chemin_lin+nom_img[i]+'_linear.png',bbox_inches='tight')
            plt.close("all")  

    def plot_ratio_bulles_windowed_vivo(self,harm_num,harm_den,nom,chemin,nbit,n=1,legend = ['data '+str(int(i))for i in range(1,11)],):
        if not os.path.isdir(chemin): # check if folder exists, otherwise create it
            os.mkdir(chemin)
        ##recuperation des indices
        n_pulse_max,n_window_max=0 , 0
        for i in range(self.n_exp):
            if self.exp[i].n_pulse>n_pulse_max:
                n_pulse_max=self.exp[i].n_pulse
            if self.exp[i].pulses[0].n_window>n_window_max:
                n_window_max=self.exp[i].pulses[0].n_window
        plot=np.zeros((6,self.n_exp,n_pulse_max,3))
        plot=np.zeros((3,self.n_exp,n_pulse_max,n_window_max))
        for i in range(self.n_exp):
            for j in range(self.exp[i].n_pulse):#self.exp[i].pulses[j].indice_harm_w[1]+
                num = np.zeros(self.exp[i].pulses[0].n_window)
                for elt in harm_num:
                    num += self.exp[i].pulses[j].indice_harm_w[elt]
                den = np.ones(self.exp[i].pulses[0].n_window)
                for elt in harm_den:
                    den += self.exp[i].pulses[j].indice_harm_w[elt]
                ratio = num / den
                plot[0,i,j,:self.exp[i].pulses[0].n_window]=20*np.log10(ratio) #20*np.log10

        self.ratio=plot
        fig=plt.figure(figsize=(20,11))
        fit=np.array([29.95777559,  6.81690626])
        amp=np.arange(nbit[0],(nbit[1]-1)*30+1)
        pression=[fit[0]*j+fit[1] for j in amp]
        self.n_window = self.exp[0].pulses[0].n_window

        mini_value = [10,25,30]
        maxi_value = [99,95,90]
        nom = ["\\ratio_"+nom+"_perLow{}_perHigh{}.png".format(mini_value[k],maxi_value[k]) for k in range(3)]
        for k in range(3):
            fig, axes = plt.subplots(nrows=1, ncols=self.n_exp, figsize=(20,11))
            mini=np.percentile(plot[0],mini_value[k])
            maxi=np.percentile(plot[0],maxi_value[k])
            for i,ax in enumerate(axes.flat):
                print(np.shape(plot[k,i]))
                im = ax.imshow(plot[0,i,:self.exp[i].n_pulse,:self.exp[i].pulses[0].n_window],aspect='auto',interpolation='none',cmap='turbo', vmin=mini,vmax=maxi)  #,extent=[0,pulse_time_ms,seq_time_s,0]   ,mapval  , vmin=mapval[i*2],vmax=mapval[i*2+1]
                ax.title.set_text(legend[i].upper() + ' logarithmic')
                ax.title.set_fontweight('bold')
                ax.set_xlabel('Intra-Pulse Time'.upper()+' (ms)')
                if i==0:
                    ax.set_ylabel('Pressure (KPa)')
            
            fig.colorbar(im,ax=axes.ravel().tolist())
            #plt.savefig(chemin+"\\bulles_H{}_fonda.png".format(k),bbox_inches='tight')
            plt.savefig(chemin+nom[k],bbox_inches='tight')
            plt.clf()
        plt.close("all")
