import numpy as np
import matplotlib.pyplot as pl

def xlabel_choice(ax,label_x='windows'):
    if label_x == 'time':
        ax.set_xlabel('Intra-Pulse Time (ms)'.upper())
    else:
        ax.set_xlabel('# TEMPORAL WINDOW'.upper())

def Plot_intrapulse_map_lin(F_mat,title,seq_time_s,pulse_time_ms,fsize):
    fig = pl.figure(figsize=(fsize[0],fsize[1]))
    pl.imshow((np.array(F_mat)),aspect='auto',interpolation='none',cmap='jet',extent=[0,pulse_time_ms,seq_time_s,0])
    pl.title(title.upper() + ' (lin scale)'.upper())
    pl.xlabel('Intrapulse Time (ms)'.upper())
    pl.ylabel('Sequence Time (s)'.upper())
    pl.colorbar()
    return fig

def Plot_intrapulse_map_dB(F_mat,title,seq_time_s,pulse_time_ms,c_min,c_max,fsize):
    fig = pl.figure(figsize=(fsize[0],fsize[1]))
    pl.imshow(20*np.log10(np.array(F_mat)),aspect='auto',interpolation='none',cmap='jet',vmin=c_min,vmax=c_max,extent=[0,pulse_time_ms,seq_time_s,0])
    pl.title(title.upper() + ' (dB scale)'.upper())
    pl.xlabel('Intrapulse Time (ms)'.upper())
    pl.ylabel('Sequence Time (s)'.upper())
    pl.colorbar()
    return fig


def plot_events(arrayPr,Pr_max,matComp1,matComp2,matComp3,arrayComp,titles,t_sequence,b_fen,born,fsize):

    t=np.linspace(0,t_sequence,len(arrayPr))
    
    fig, (ax1, ax2, ax3, ax4, ax5) = pl.subplots(1,5, gridspec_kw={'width_ratios': [1,2,2,2,1]}, figsize=(fsize[0],fsize[1]))
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    ax1.set_ylim(t_sequence,0)
    # ax1.set_xlim(1.1*Pr_max,0)
    ax1.fill_betweenx(t,arrayPr)
    ax1.plot([Pr_max,Pr_max],[t[0],t[-1]],color='r')
    ax1.set_ylabel('SEQUENCE TIME (s)')
    ax1.set_xlabel(titles[0])
    ax1.grid('on')
    im = ax2.imshow(matComp1,aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax2.set_title(titles[1][0])
#    ax2.scatter()
    #xlabel_choice(ax2)    
    im = ax3.imshow(matComp2,aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax3.set_title(titles[1][1])
#    ax2.scatter()
    #xlabel_choice(ax2)    
    im = ax4.imshow(matComp3,aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax4.set_title(titles[1][2])
#    ax2.scatter()
    #xlabel_choice(ax2)
    ax5.plot((arrayComp),t)
    # ax3.plot(20*np.log10(array1),t)
    ax5.invert_yaxis()
    ax5.set_ylim(t_sequence,0)
    # ax3.set_xlim(0,15)
    ax5.set_xlabel(titles[2])
    ax5.grid('on')
    # fig.colorbar(im,ax=ax2)
    # fig.tight_layout()
    # ##
    return fig, (ax1, ax2, ax3, ax4, ax5)


def Subplot_intrapulse_4maps_dB(Mats,titles,seq_time_s,b_fen,born,fsize):  
    i=0
    fig, axes = pl.subplots(nrows=1, ncols=4, figsize=(fsize[0],fsize[1]))
    for ax in axes.flat:
        im = ax.imshow(20*np.log10(np.array(Mats[i])),aspect='auto',interpolation='none',cmap='jet',vmin=born[0],vmax=born[1],extent=[b_fen[0],b_fen[1],seq_time_s,0])
        ax.title.set_text(titles[i].upper() + ' (dB scale)'.upper())
        xlabel_choice(ax)
        if i==0:
            ax.set_ylabel('Sequence Time (s)'.upper())
        i+=1
        
    fig.colorbar(im,ax=axes.ravel().tolist())
    pl.show()
    return fig, axes

def Subplot_intrapulse_3maps_dB(Mats,titles,seq_time_s,pulse_time_ms,c_min,c_max,fsize):   
    i=0
    fig, axes = pl.subplots(nrows=1, ncols=3, figsize=(fsize[0],fsize[1]))
    for ax in axes.flat:
        im = ax.imshow(20*np.log10(np.array(Mats[i])),aspect='auto',interpolation='none',cmap='jet',vmin=c_min,vmax=c_max,extent=[0,pulse_time_ms,seq_time_s,0])
        ax.title.set_text(titles[i].upper() + ' (dB scale)'.upper())
        xlabel_choice(ax)
        if i==0:
            ax.set_ylabel('Sequence Time (s)'.upper())
        i+=1
        
    fig.colorbar(im,ax=axes.ravel().tolist())
    pl.show()
    return fig, axes

def Subplot_intrapulse_2maps_dB(Mats,titles,seq_time_s,b_fen,c_min,c_max,fsize):   
    i=0
    fig, axes = pl.subplots(nrows=1, ncols=2, figsize=(fsize[0],fsize[1]))
    for ax in axes.flat:
        im = ax.imshow(20*np.log10(np.array(Mats[i])),aspect='auto',interpolation='none',cmap='jet',vmin=c_min,vmax=c_max,extent=[b_fen[0],b_fen[1],seq_time_s,0])
        ax.title.set_text(titles[i].upper() + ' (dB scale)'.upper())
        xlabel_choice(ax)
        if i==0:
            ax.set_ylabel('Sequence Time (s)'.upper())
        i+=1
        
    fig.colorbar(im,ax=axes.ravel().tolist())
    pl.show()
    return fig, axes





def Subplot_Map_Hist(Mat,Mat_mean,titles,Dur_seq,Dur_pulse,min_val,max_val,fsize):
    fig, (ax2, ax3) = pl.subplots(1,2, gridspec_kw={'width_ratios': [3,1]}, figsize=(fsize[0],fsize[1]))
    
    im = ax2.imshow((Mat),aspect='auto',interpolation='none',cmap='jet',vmin=min_val,vmax=max_val,extent=[0,Dur_pulse,Dur_seq,0])
    ax2.set_title(titles[0])
    xlabel_choice(ax2)
    
    t=np.linspace(0,Dur_seq,len(Mat_mean))
    
    ax3.fill_betweenx(t,Mat_mean,min_val)
    # ax3.plot(20*np.log10(array1),t)
    ax3.invert_yaxis()
    ax3.set_ylim(Dur_seq,0)
    ax3.set_xlim(min_val,max_val)
    ax3.set_xlabel(titles[1])
    ax3.grid('on')
    
    fig.colorbar(im,ax=ax2)
    fig.tight_layout()
    return fig


def plot_IndexMap_and_Pressure(arrayPr,Pr_max,matComp,arrayComp,titles,t_sequence,b_fen,born,fsize):

    t=np.linspace(0,t_sequence,len(arrayPr))
    
    fig, (ax1, ax2, ax3) = pl.subplots(1,3, gridspec_kw={'width_ratios': [1,3,1]}, figsize=(fsize[0],fsize[1]))
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    ax1.set_ylim(t_sequence,0)
    # ax1.set_xlim(1.1*Pr_max,0)
    ax1.fill_betweenx(t,arrayPr)
    ax1.plot([Pr_max,Pr_max],[t[0],t[-1]],color='r')
    ax1.set_ylabel('SEQUENCE TIME (s)')
    ax1.set_xlabel(titles[0])
    ax1.grid('on')
    
    im = ax2.imshow(20*np.log10(matComp),aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax2.set_title(titles[1]+' - dB')
    xlabel_choice(ax2)
    
    ax3.plot((arrayComp),t)
    # ax3.plot(20*np.log10(array1),t)
    ax3.invert_yaxis()
    ax3.set_ylim(t_sequence,0)
    # ax3.set_xlim(0,15)
    ax3.set_xlabel(titles[2])
    ax3.grid('on')
    
    fig.colorbar(im,ax=ax2)
    fig.tight_layout()
    ##
    return fig, (ax1, ax2, ax3)

def plot_Pressure_IndexMap_MeanMax(arrayPr,Pr_max,matComp,titles,t_sequence,b_fen,born,fsize):

    meanComp=np.mean(matComp,axis=1)
    maxComp=np.max(matComp,axis=1)
    t=np.linspace(0,t_sequence,len(meanComp))
    
    fig, (ax1, ax2, ax3) = pl.subplots(1,3, gridspec_kw={'width_ratios': [1,3,1]}, figsize=(fsize[0],fsize[1]))
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    ax1.set_ylim(t_sequence,0)
    # ax1.set_xlim(1.1*Pr_max,0)
    ax1.fill_betweenx(t,arrayPr)
    ax1.plot([Pr_max,Pr_max],[t[0],t[-1]],color='r')
    ax1.set_ylabel('SEQUENCE TIME (s)')
    ax1.set_xlabel(titles[0])
    ax1.grid('on')
    
    im = ax2.imshow(20*np.log10(matComp),aspect='auto',interpolation='none',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax2.set_title(titles[1]+' - dB')
    xlabel_choice(ax2)
    
    ax3.plot((meanComp),t,color='k')
    ax3.plot((maxComp),t,color='r')
    # ax3.plot(20*np.log10(array1),t)
    ax3.invert_yaxis()
    ax3.set_ylim(t_sequence,0)
    # ax3.set_xlim(0,15)
    ax3.set_xlabel(titles[2])
    ax3.grid('on')
    
    fig.colorbar(im,ax=ax2)
    fig.tight_layout()
    ##
    return fig, (ax1, ax2, ax3)

def plot_IndexMap(matComp,arrayComp,titles,t_sequence,b_fen,born,fsize):

    t=np.linspace(0,t_sequence,len(arrayComp))
    
    fig, (ax2, ax3) = pl.subplots(1,2, gridspec_kw={'width_ratios': [3,1]}, figsize=(fsize[0],fsize[1]))

    im = ax2.imshow(20*np.log10(matComp),aspect='auto',interpolation='bilinear',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax2.set_title(titles[0]+' - dB')
    xlabel_choice(ax2)
    ax2.set_ylabel('SEQUENCE TIME (s)')
    
    ax3.plot((arrayComp),t)
    # ax3.plot(20*np.log10(array1),t)
    ax3.invert_yaxis()
    ax3.set_ylim(t_sequence,0)
    # ax3.set_xlim(0,15)
    ax3.set_xlabel(titles[1])
    ax3.grid('on')
    
    fig.colorbar(im,ax=ax2)
    fig.tight_layout()
    ##
    return fig, ax2, ax3

def plot_IndexMap_mean_max(matComp,titles,t_sequence,b_fen,born,fsize):

    meanComp=np.mean(matComp,axis=1)
    maxComp=np.max(matComp,axis=1)
    t=np.linspace(0,t_sequence,len(meanComp))
    
    fig, (ax2, ax3) = pl.subplots(1,2, gridspec_kw={'width_ratios': [3,1]}, figsize=(fsize[0],fsize[1]))

    im = ax2.imshow(20*np.log10(matComp),aspect='auto',interpolation='bilinear',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax2.set_title(titles[0]+' - dB')
    xlabel_choice(ax2)
    ax2.set_ylabel('SEQUENCE TIME (s)')
    
    ax3.plot((meanComp),t,color='k')
    ax3.plot((maxComp),t,color='r')
    # ax3.plot(20*np.log10(array1),t)
    ax3.invert_yaxis()
    ax3.set_ylim(t_sequence,0)
    # ax3.set_xlim(0,15)
    ax3.set_xlabel(titles[1])
    ax3.grid('on')
    
    fig.colorbar(im,ax=ax2)
    fig.tight_layout()
    ##
    return fig, ax2, ax3

def plot_IndexMap_lin(matComp,arrayComp,titles,t_sequence,b_fen,born,fsize):

    t=np.linspace(0,t_sequence,len(arrayComp))
    
    fig, (ax2, ax3) = pl.subplots(1,2, gridspec_kw={'width_ratios': [3,1]}, figsize=(fsize[0],fsize[1]))

    im = ax2.imshow((matComp),aspect='auto',interpolation='bilinear',cmap='jet',extent=[b_fen[0],b_fen[1],t_sequence,0],vmin=born[0],vmax=born[1])
    ax2.set_title(titles[0]+' - dB')
    ax2.set_xlabel('# INTRA-PULSE TEMPORAL WINDOW')
    ax2.set_ylabel('SEQUENCE TIME (s)')
    
    ax3.plot((arrayComp),t)
    # ax3.plot(20*np.log10(array1),t)
    ax3.invert_yaxis()
    ax3.set_ylim(t_sequence,0)
    # ax3.set_xlim(0,15)
    ax3.set_xlabel(titles[1])
    ax3.grid('on')
    
    fig.colorbar(im,ax=ax2)
    fig.tight_layout()
    ##
    return fig, ax2, ax3


# 
def Scatter_and_Hists_wProc(BN_norm,UH2_norm,list_born,Nb_bins,fsize,list_thr=[]):
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.025
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    
    Val0=np.array(BN_norm).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts=Val
    (Xic,Yic,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    pl.close();

    Val0=np.array(UH2_norm).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts=Val
    (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    pl.close();
         
    # start with a square Figure
    fig_scatter = pl.figure(figsize=(fsize[0],fsize[1]))
    
    ax = fig_scatter.add_axes(rect_scatter)
    ax_histx = fig_scatter.add_axes(rect_histx, sharex=ax)
    ax_histy = fig_scatter.add_axes(rect_histy, sharey=ax)
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
                 
    ax_histx.bar(Yic[1:],np.log10(Xic+1),0.2,color='k')
    ax_histy.barh(Ysc2[1:],np.log10(Xsc2+1),0.2,color='k')
    ax.plot(20*np.log10(BN_pts),20*np.log10(iUH2_pts),color='k',linestyle='None',marker='+',markersize=8)
    ax.grid('on')  
    
    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    ax.set_xlim([bxmin,bxmax])
    ax_histx.set_xlim([bxmin,bxmax])
    ax.set_ylim([bymin,bymax])
    ax_histy.set_ylim([bymin,bymax])
    if len(list_thr)!=0:
        thr_y=list_thr[1]
        thr_x=list_thr[0]
        ax.plot([bxmin,bxmax],[thr_y,thr_y],color='k')
        ax.plot([thr_x,thr_x],[bymin,bymax],color='k')   
        
    return fig_scatter,(ax,ax_histx,ax_histy)  

def Scatter_simple_mult(ax,BN_norm,UH2_norm,list_born,list_thr=[]):
    
    Val0=np.array(BN_norm).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts=Val

    Val0=np.array(UH2_norm).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts=Val
   
    # start with a square Figure
    ax.plot(20*np.log10(BN_pts),20*np.log10(iUH2_pts),linestyle='None',marker='.',markersize=5,color='k')      
    ax.grid('on')

    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    ax.set_xlim([bxmin,bxmax])
    ax.set_ylim([bymin,bymax])
    if len(list_thr)!=0:
        thr_y=list_thr[1]
        thr_x=list_thr[0]
        ax.plot([bxmin,bxmax],[thr_y,thr_y],color='k')
        ax.plot([thr_x,thr_x],[bymin,bymax],color='k')   

def Scatter_simple_2(BN_norm,UH2_norm,list_born,fsize,legend=['',''],list_thr=[]):  
  
    Val0=np.array(BN_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts=Val

    Val0=np.array(UH2_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts=Val

    Val0=np.array(BN_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_2=Val

    Val0=np.array(UH2_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_2=Val
         
    # start with a square Figure
    fig_scatter,ax = pl.subplots(1,1,figsize=(fsize[0],fsize[1]))
                
    ax.plot(20*np.log10(BN_pts),20*np.log10(iUH2_pts),label=legend[0],linestyle='None',marker='.',markersize=5,color='k')
    ax.plot(20*np.log10(BN_pts_2),20*np.log10(iUH2_pts_2),label=legend[1],linestyle='None',marker='+',markersize=14,mew=3,color='indianred')      
    ax.grid('on')  
    
    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    ax.set_xlim([bxmin,bxmax])
    ax.set_ylim([bymin,bymax])
    if len(list_thr)!=0:
        thr_y=list_thr[1]
        thr_x=list_thr[0]
        ax.plot([bxmin,bxmax],[thr_y,thr_y],color='k')
        ax.plot([thr_x,thr_x],[bymin,bymax],color='k')   
        
    return fig_scatter,ax   


def Scatter_simple_2_mult(ax,BN_norm,UH2_norm,list_born,legend=['',''],list_thr=[]):  
  
    Val0=np.array(BN_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts=Val

    Val0=np.array(UH2_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts=Val

    Val0=np.array(BN_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_2=Val

    Val0=np.array(UH2_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_2=Val
         
    # start with a square Figure
                
    ax.plot(20*np.log10(BN_pts),20*np.log10(iUH2_pts),label=legend[0],linestyle='None',marker='.',markersize=5,color='k')
    ax.plot(20*np.log10(BN_pts_2),20*np.log10(iUH2_pts_2),label=legend[1],linestyle='None',marker='+',markersize=14,mew=3,color='indianred')      
    ax.grid('on')  
    
    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    ax.set_xlim([bxmin,bxmax])
    ax.set_ylim([bymin,bymax])
    if len(list_thr)!=0:
        thr_y=list_thr[1]
        thr_x=list_thr[0]
        ax.plot([bxmin,bxmax],[thr_y,thr_y],color='k')
        ax.plot([thr_x,thr_x],[bymin,bymax],color='k')   


def Scatter_simple_Paul(ax,BN_norm,UH2_norm,list_born,legend=['',''],list_thr=[]):  
  
    Val0=np.array(BN_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts= 20 * np.log10(abs(Val))

    Val0=np.array(UH2_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts= 20 * np.log10(abs(Val))

    Val0=np.array(BN_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_2= 20 * np.log10(abs(Val))

    Val0=np.array(UH2_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_2= 20 * np.log10(abs(Val))
         
    # start with a square Figure
    ax.plot(BN_pts,iUH2_pts,label=legend[0],linestyle='None',marker='o',markersize=5,color='indianred', alpha = 0.2)
    ax.plot(BN_pts_2,iUH2_pts_2,label=legend[1],linestyle='None',marker='o',markersize=5,color='indianred', alpha = 0.2)    
    ax.grid('on')  
    
    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    # ax.set_xlim([bxmin,bxmax])
    # ax.set_ylim([bymin,bymax])
    # if len(list_thr)!=0:
    #     thr_y=list_thr[1]
    #     thr_x=list_thr[0]
    #     ax.plot([bxmin,bxmax],[thr_y,thr_y],color='k')
    #     ax.plot([thr_x,thr_x],[bymin,bymax],color='k')   
        
    
def Scatter_and_Hists_wProc_2(BN_norm,UH2_norm,list_born,Nb_bins,fsize,legend=['',''],list_thr=[]):
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.025
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
      
    Val0=np.array(BN_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts=Val
    (Xic,Yic,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    # (Xic,Yic,p)=pl.hist(20*np.log10(Val),200);
    pl.close();

    Val0=np.array(UH2_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts=Val
    (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    # (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),400);
    pl.close();
                  
    Val0=np.array(BN_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_2=Val
    (Xic_2,Yic_2,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    # (Xic,Yic,p)=pl.hist(20*np.log10(Val),200);
    pl.close();

    Val0=np.array(UH2_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_2=Val
    (Xsc2_2,Ysc2_2,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    # (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),400);
    pl.close();
   
    # start with a square Figure
    fig_scatter = pl.figure(figsize=(fsize[0],fsize[1]))
    
    ax = fig_scatter.add_axes(rect_scatter)
    ax_histx = fig_scatter.add_axes(rect_histx, sharex=ax)
    ax_histy = fig_scatter.add_axes(rect_histy, sharey=ax)
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
   
    ax_histx.bar(Yic[1:],np.log10(Xic+1),0.2,color='k')
    ax_histy.barh(Ysc2[1:],np.log10(Xsc2+1),0.2,color='k')
    ax.plot(20*np.log10(BN_pts),20*np.log10(iUH2_pts),label=legend[0],linestyle='None',marker='.',markersize=5,color='k')
    
    ax_histx.bar(Yic_2[1:],np.log10(Xic_2+1),0.2,color='indianred')
    ax_histy.barh(Ysc2_2[1:],np.log10(Xsc2_2+1),0.2,color='indianred')
    ax.plot(20*np.log10(BN_pts_2),20*np.log10(iUH2_pts_2),label=legend[1],linestyle='None',marker='+',markersize=14,mew=3,color='indianred')    
    
    ax.grid('on')
    ax_histx.grid('on')
    ax_histy.grid('on')    
    
    # fig_scatter.tight_layout
    
    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    ax.set_xlim([bxmin,bxmax])
    ax_histx.set_xlim([bxmin,bxmax])
    ax.set_ylim([bymin,bymax])
    ax_histy.set_ylim([bymin,bymax])
    if len(list_thr)!=0:
        thr_y=list_thr[1]
        thr_x=list_thr[0]
        ax.plot([bxmin,bxmax],[thr_y,thr_y],color='k')
        ax.plot([thr_x,thr_x],[bymin,bymax],color='k')   
     
    BN_return=[BN_pts[np.isfinite(BN_pts)],BN_pts_2[np.isfinite(BN_pts_2)]]
    iUH_return=[iUH2_pts[np.isfinite(iUH2_pts)],iUH2_pts_2[np.isfinite(iUH2_pts_2)]]
        
    return fig_scatter,(ax,ax_histx,ax_histy),BN_return, iUH_return

def Scatter_and_Hists_wProc_5(BN_norm,UH2_norm,labels,list_born,Nb_bins,fsize,list_thr=[]):
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.025
    mrkr_size= 4
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
      
    Val0=np.array(BN_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts=Val
    (Xic,Yic,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    # (Xic,Yic,p)=pl.hist(20*np.log10(Val),200);
    pl.close();

    Val0=np.array(UH2_norm[0]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts=Val
    (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    # (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),400);
    pl.close();
                  
    Val0=np.array(BN_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_2=Val
    (Xic_2,Yic_2,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    # (Xic,Yic,p)=pl.hist(20*np.log10(Val),200);
    pl.close();

    Val0=np.array(UH2_norm[1]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_2=Val
    (Xsc2_2,Ysc2_2,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    # (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),400);
    pl.close();
    
    Val0=np.array(BN_norm[2]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_3=Val
    (Xic_3,Yic_3,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    # (Xic,Yic,p)=pl.hist(20*np.log10(Val),200);
    pl.close();

    Val0=np.array(UH2_norm[2]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_3=Val
    (Xsc2_3,Ysc2_3,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    # (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),400);
    pl.close();
 
    
    Val0=np.array(BN_norm[3]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_4=Val
    (Xic_4,Yic_4,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    # (Xic,Yic,p)=pl.hist(20*np.log10(Val),200);
    pl.close();

    Val0=np.array(UH2_norm[3]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_4=Val
    (Xsc2_4,Ysc2_4,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    # (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),400);
    pl.close();
    
    Val0=np.array(BN_norm[4]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan)
    Val=np.where(Val1!=1,Val1,np.nan)
    BN_pts_5=Val
    (Xic_5,Yic_5,p)=pl.hist(20*np.log10(Val),Nb_bins[0]);
    # (Xic,Yic,p)=pl.hist(20*np.log10(Val),200);
    pl.close();

    Val0=np.array(UH2_norm[4]).flatten()
    Val1=np.where(Val0!=0,Val0,np.nan) 
    Val=np.where(Val1!=1,Val1,np.nan)
    iUH2_pts_5=Val
    (Xsc2_5,Ysc2_5,p)=pl.hist(20*np.log10(Val),Nb_bins[1]);
    # (Xsc2,Ysc2,p)=pl.hist(20*np.log10(Val),400);
    pl.close();
        
    # start with a square Figure
    fig_scatter = pl.figure(figsize=(fsize[0],fsize[1]))
    
    ax = fig_scatter.add_axes(rect_scatter)
    ax_histx = fig_scatter.add_axes(rect_histx, sharex=ax)
    ax_histy = fig_scatter.add_axes(rect_histy, sharey=ax)
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    ax_histx.bar(Yic[1:],np.log10(Xic+1),0.2,color='k')
    ax_histy.barh(Ysc2[1:],np.log10(Xsc2+1),0.2,color='k')
    ax.plot(20*np.log10(BN_pts),20*np.log10(iUH2_pts),label=labels[0],linestyle='None',marker='.',markersize=mrkr_size,color='k')
    
    ax_histx.bar(Yic_2[1:],np.log10(Xic_2+1),0.2,color='dimgrey')
    ax_histy.barh(Ysc2_2[1:],np.log10(Xsc2_2+1),0.2,color='dimgrey')
    ax.plot(20*np.log10(BN_pts_2),20*np.log10(iUH2_pts_2),label=labels[1],linestyle='None',marker='.',markersize=mrkr_size,color='dimgrey')    

    ax_histx.bar(Yic_3[1:],np.log10(Xic_3+1),0.2,color='grey')
    ax_histy.barh(Ysc2_3[1:],np.log10(Xsc2_3+1),0.2,color='grey')
    ax.plot(20*np.log10(BN_pts_3),20*np.log10(iUH2_pts_3),label=labels[2],linestyle='None',marker='.',markersize=mrkr_size,color='grey')    

    ax_histx.bar(Yic_4[1:],np.log10(Xic_4+1),0.2,color='darkgrey')
    ax_histy.barh(Ysc2_4[1:],np.log10(Xsc2_4+1),0.2,color='darkgrey')
    ax.plot(20*np.log10(BN_pts_4),20*np.log10(iUH2_pts_4),label=labels[3],linestyle='None',marker='.',markersize=mrkr_size,color='darkgrey')    

    ax_histx.bar(Yic_5[1:],np.log10(Xic_5+1),0.2,color='lightgrey')
    ax_histy.barh(Ysc2_5[1:],np.log10(Xsc2_5+1),0.2,color='lightgrey')
    ax.plot(20*np.log10(BN_pts_5),20*np.log10(iUH2_pts_5),label=labels[4],linestyle='None',marker='.',markersize=mrkr_size,color='lightgrey')    

    ax.grid('on')
    
    # fig_scatter.tight_layout
    
    bxmin=list_born[0]#-6
    bxmax=list_born[1]#24
    bymin=list_born[2]#-14
    bymax=list_born[3]#34
    ax.set_xlim([bxmin,bxmax])
    ax_histx.set_xlim([bxmin,bxmax])
    ax.set_ylim([bymin,bymax])
    ax_histy.set_ylim([bymin,bymax])
    if len(list_thr)!=0:
        thr_y=list_thr[1]
        thr_x=list_thr[0]
        ax.plot([bxmin,bxmax],[thr_y,thr_y],color='k')
        ax.plot([thr_x,thr_x],[bymin,bymax],color='k')   
     
    
    pl.tight_layout()
        
    return fig_scatter,(ax,ax_histx,ax_histy)

def Subplot_Cavitation_Occurence(T_INTRA, P_VEC_f, T_VEC, T_VEC_f, arrays, figsize=(18, 10)):

    t_vec = arrays
    [t_vec, arrayPr] = arrays
    
    [T_INTRA_IC, T_INTRA_UH, T_INTRA_BOTH] = T_INTRA
    [P_VEC_IC_f, P_VEC_UH_f, P_VEC_BOTH_f] = P_VEC_f
    [T_VEC_IC, T_VEC_UH, T_VEC_BOTH] = T_VEC
    [T_VEC_IC_f, T_VEC_UH_f, T_VEC_BOTH_f] = T_VEC_f
    
    (N_UH,)=T_VEC_UH_f.shape
    (N_IC,)=T_VEC_IC_f.shape
    (N_both,)=T_VEC_BOTH_f.shape
    
    (N_UH_t,)=T_VEC_UH.shape
    (N_IC_t,)=T_VEC_IC.shape
    (N_both_t,)=T_VEC_BOTH.shape

    title_t="Total events - i-ICD: "+str(N_IC_t) +" - i-UCD: "+str(N_UH_t) +" - both: "+str(N_both_t)
    title_f="First appearing events - i-ICD: "+str(N_IC) +" - i-UCD: "+str(N_UH) +" - both: "+str(N_both)
    
    fig,(ax1,ax2)=pl.subplots(2,1,figsize=figsize)
    ax1.plot(t_vec,arrayPr,color='k',linestyle='dotted')
    ax1.plot(T_VEC_UH_f,P_VEC_UH_f,marker='.',linestyle="None",markersize=8,color='red',label='i-UCD Events')
    ax1.plot(T_VEC_IC_f,P_VEC_IC_f,marker='.',linestyle="None",markersize=8,color='k',label='i-ICD Events')
    ax1.plot(T_VEC_BOTH_f,P_VEC_BOTH_f,marker='x',linestyle="None",markersize=8,color='k',label='i-UCD & i-ICD Events')
    # ax.plot(t_vec,arrayPr,color='k')
    # pl.grid('on')
    ax1.set_xlim([0,1.01*t_vec[-1]])
    ax1.set_ylabel('Peak negetive pressure (MPa)')
    # pl.xlabel('Sequence time (s)')
    ax1.legend(loc='lower center',ncol=3)
    ax1.grid('on')
    ax1.set_title(title_f)
    # pl.savefig('Pressure_Evnts.png')


    ax2.plot(T_VEC_UH,T_INTRA_UH*1e-3,marker='.',linestyle="None",markersize=4,color='red',label='i-UCD Events')
    ax2.plot(T_VEC_IC,T_INTRA_IC*1e-3,marker='.',linestyle="None",markersize=4,color='k',label='i-ICD Events')
    ax2.plot(T_VEC_BOTH,T_INTRA_BOTH*1e-3,marker='x',linestyle="None",markersize=8,color='k',label='i-UCD & i-ICD Events')
    # ax.plot(t_vec,arrayPr,color='k')
    # pl.grid('on')
    # pl.axis([0,t_sequence,0,Pr_max*1.1])
    ax2.set_ylabel('Intrapulse time (ms)')
    ax2.set_xlabel('Sequence time (s)')
    ax2.set_xlim([0,1.01*t_vec[-1]])
    ax2.grid('on')
    ax2.set_title(title_t)
    # ax.legend(bbox_to_anchor=(0.93, 1.2),ncol=3, fancybox=True)
    pl.tight_layout
    
    return fig,(ax1,ax2)

def Subplot_Cavitation_Occurence_mult(axes,T_INTRA, P_VEC_f, T_VEC, T_VEC_f, arrays):

    t_vec = arrays
    [t_vec, arrayPr] = arrays
    
    [T_INTRA_IC, T_INTRA_UH, T_INTRA_BOTH] = T_INTRA
    [P_VEC_IC_f, P_VEC_UH_f, P_VEC_BOTH_f] = P_VEC_f
    [T_VEC_IC, T_VEC_UH, T_VEC_BOTH] = T_VEC
    [T_VEC_IC_f, T_VEC_UH_f, T_VEC_BOTH_f] = T_VEC_f
    
    (N_UH,)=T_VEC_UH_f.shape
    (N_IC,)=T_VEC_IC_f.shape
    (N_both,)=T_VEC_BOTH_f.shape
    
    (N_UH_t,)=T_VEC_UH.shape
    (N_IC_t,)=T_VEC_IC.shape
    (N_both_t,)=T_VEC_BOTH.shape

    title_t="Total events - i-ICD: "+str(N_IC_t) +" - i-UCD: "+str(N_UH_t) +" - both: "+str(N_both_t)
    title_f="First appearing events - i-ICD: "+str(N_IC) +" - i-UCD: "+str(N_UH) +" - both: "+str(N_both)
    
    ax1=axes[0]
    ax1.plot(t_vec,arrayPr,color='k',linestyle='dotted')
    ax1.plot(T_VEC_UH_f,P_VEC_UH_f,marker='.',linestyle="None",markersize=8,color='red',label='i-UCD Events')
    ax1.plot(T_VEC_IC_f,P_VEC_IC_f,marker='.',linestyle="None",markersize=8,color='k',label='i-ICD Events')
    ax1.plot(T_VEC_BOTH_f,P_VEC_BOTH_f,marker='x',linestyle="None",markersize=8,color='k',label='i-UCD & i-ICD Events')
    # ax.plot(t_vec,arrayPr,color='k')
    # pl.grid('on')
    ax1.set_xlim([0,1.01*t_vec[-1]])
    ax1.set_ylabel('Peak negetive pressure (MPa)')
    # pl.xlabel('Sequence time (s)')
    ax1.legend(loc='lower center',ncol=3)
    ax1.grid('on')
    ax1.set_title(title_f)
    # pl.savefig('Pressure_Evnts.png')

    ax2=axes[1]
    ax2.plot(T_VEC_UH,T_INTRA_UH*1e-3,marker='.',linestyle="None",markersize=4,color='red',label='i-UCD Events')
    ax2.plot(T_VEC_IC,T_INTRA_IC*1e-3,marker='.',linestyle="None",markersize=4,color='k',label='i-ICD Events')
    ax2.plot(T_VEC_BOTH,T_INTRA_BOTH*1e-3,marker='x',linestyle="None",markersize=8,color='k',label='i-UCD & i-ICD Events')
    # ax.plot(t_vec,arrayPr,color='k')
    # pl.grid('on')
    # pl.axis([0,t_sequence,0,Pr_max*1.1])
    ax2.set_ylabel('Intrapulse time (ms)')
    ax2.set_xlabel('Sequence time (s)')
    ax2.set_xlim([0,1.01*t_vec[-1]])
    ax2.grid('on')
    ax2.set_title(title_t)
    # ax.legend(bbox_to_anchor=(0.93, 1.2),ncol=3, fancybox=True)
    