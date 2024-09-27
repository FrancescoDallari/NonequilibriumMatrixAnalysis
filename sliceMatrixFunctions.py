# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:27:19 2024

@author: francesco.dallari
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def sliceMatrix(lagmin,lagmax,TT_in,mode='fast'):
    """
    Function to extract an effective g2 from a single correlation matrix
    the slice is extracted between le age times lagmin and lagmax
    
    Input:
        lagmin = smallest age
        lagmax = largest age
        TT_in = two-times correlation matrix
        mode = way in which the cut is performed:
                'fast' -> t_w = t1 ; t = |t2-t1|
                'ACS' -> t_w=(t1+t2)/2 ; t = |t2-t1|
                'classic' -> calculate the g2 over a submatrix between lagmin and lagmax
    
    output:
        g2 : 1d array , might contain NaN when in 'ACS' mode
    """
    
    
    TTcopy =  np.copy(np.squeeze(TT_in))#
    
    if mode =='ACS':  
        g2=np.zeros(TTcopy.shape[-1])*np.nan
        tw = np.round((lagmax+lagmin)/2).astype(int)
        dt = np.round((lagmax-lagmin)/2).astype(int)
        #dummy = np.zeros_like(tweaked)
        for j in range(TTcopy.shape[0]):
            t1 = np.floor(tw-dt-j/2).astype(int)
            t2 = np.ceil(tw+dt -j/2).astype(int)
                    
            if (t1<0) & (t2>0):
                g2[j] = np.nanmean(np.diag(TTcopy,j)[0:t2])
                        

            elif t2 >TTcopy.shape[0]:
                g2[j] = np.nanmean(np.diag(TTcopy,j)[t1:])
            
            elif (t2 <=TTcopy.shape[0]) & (t1>=0):
                g2[j] = np.nanmean(np.diag(TTcopy,j)[t1:t2])
                        
            
            else:
                g2[j] = np.nan
        
    elif mode=='fast':
        tw = lagmin#np.round((lagmax+lagmin)/2).astype(int)
        dt = (lagmax-lagmin)#.astype(int)
        g2=np.asarray([np.nanmean(np.diag(TTcopy,t_i)[tw:(tw+dt)]) for t_i in range(TTcopy.shape[0]-(tw+1)) ])            
    else:
        #classic
        tweaked =  np.squeeze(TT_in[lagmin:lagmax,lagmin:lagmax])#
        g2 = np.asarray([np.nanmean(np.diag(tweaked,j)) for j in range(tweaked.shape[0]) ])
     
    return g2


def analyzeAgedepTT(TT_in,tStep,lastTw,firstTw=0,usefulBaseline = None,method = 'fast',showPlot=False,userbounds=None):
    """
    Function to obtain the age evolution from a single matrix 
    ATTENTION! ALL VALUES CALCULATED IN FRAMES
    
    input:
        TT_in = Autocorellation matrix; it should be the g2(t1,t2)-1 !
        tStep = step between cuts
        lastTw= maximum age to be used. If the mode is not 'classic' it's better to keep id a bit smaller than TT_in.shape[0]
     optional:
         firstTw = first frame to be considered
         usefulBaseline = can be None, 'auto' or a number; if None the fit is done normally, if 'Auto' a baseline value is calculated taking the average of the corner (last 20 lags) of the TT_in; if it's a number the baseline in the g2 willbe fixed at that value 
         method = way in which the cut is performed:
                'fast' -> t_w = t1 ; t = |t2-t1|
                'ACS' -> t_w=(t1+t2)/2 ; t = |t2-t1|
                'classic' -> calculate the g2 over a submatrix between lagmin and lagmax
         userbounds = boundaries for the fitting function, if it's not None it has to be a touple with the bounds for: contrast, relax. rate, exponent, baseline
    
    output:
        outDset =  a dataset containing the results of a KWW fit and the age axis
    """
    allResults = []
    commonBaseline = []
    col = plt.cm.jet(np.linspace(0,1,len(range(firstTw,lastTw,tStep))) )
    if showPlot:
        plt.figure()
    
    
    
    if usefulBaseline =='auto':
        dummy = []
        for j in range(TT_in.shape[0]-20,TT_in.shape[0]):
            dummy.append(np.nanmean(np.diag(TT_in,j)))
        bsl = np.mean(dummy)
    
    for t_i,t_v in enumerate(range(firstTw,lastTw,tStep)):
        
        gg = sliceMatrix(t_v, t_v+tStep, TT_in, method)
        tAx = np.arange(0,len(gg),1)
        
        startVals = [np.nanmean(gg[2:10]),1/tStep,1,0.002]
        
        if userbounds:
            bounds = userbounds 
        else:
            bounds = ([0, 0, 0.5, -.001],
                      [0.8, np.inf, 2, np.nanmean(gg[2:10])*0.8 ])
        
        
        for i in range(len(startVals)):
            if (startVals[i]<bounds[0][i]) | (startVals[i]>bounds[1][i]):
                startVals[i] = (bounds[0][i]+bounds[1][i])/2         
        startVals=tuple(startVals)
        
        x_out,y_out,dx_out,dy_out = binandclean(tAx[1:], gg[1:], nbins = 100) 
        if showPlot:
            plt.semilogx(x_out,y_out,'.-',color=col[t_i],alpha =.25)
            plt.xlabel('frame')
            plt.ylabel('$g_2$-1')

        try:
            if usefulBaseline == None:
                rep2,covp2 = optimize.curve_fit(expStretch,x_out[~np.isnan(y_out)],y_out[~np.isnan(y_out)],
                                                        p0=startVals,
                                                             bounds=bounds)
                rep3 = rep2
                err3 = np.sqrt(np.diag(covp2))
                
            elif usefulBaseline =='auto':
                
                #bsl = np.mean(usedbsl[t_i:t_i+tStep+overlap])
                #if np.isnan(bsl):
                #    bsl = prevbsl
                #else:
                #    prevbsl = bsl
                                        
                fixblsexpstr = lambda x,c,g,b : expStretch(x,c,g,b,bsl)
                
                rep2,covp2 = optimize.curve_fit(fixblsexpstr,x_out ,y_out ,
                                                        p0=startVals[:-1],
                                                             bounds= (bounds[0][:-1],bounds[1][:-1]) ) #([0, 0, 0.5 ],
                                                                          #[0.8, np.inf, 2  ]))
                
                rep3 = np.asarray([rep2[0],rep2[1],rep2[2],bsl])
                err2 = np.sqrt(np.diag(covp2))
                err3 = np.asarray([err2[0],err2[1],err2[2],0])

                
            else:
                # if usefulBaseline is not a number here will fail
                fixblsexpstr = lambda x,c,g,b : expStretch(x,c,g,b,usefulBaseline)
                rep2,covp2 = optimize.curve_fit(fixblsexpstr,x_out ,y_out ,
                                                        p0=startVals[:-1],
                                                             bounds=(bounds[0][:-1],bounds[1][:-1]) )#([0, 0, 0.5 ],
                                                                     #[0.8, np.inf, 2  ]))
                
                rep3 = np.asarray([rep2[0],rep2[1],rep2[2],bsl])
                err2 = np.sqrt(np.diag(covp2))
                err3 = np.asarray([err2[0],err2[1],err2[2],0])

#                 rep2,covp2 = optimize.curve_fit(expStretch,x_out ,y_out ,
#                                                         p0=(np.nanmean(gg[2:10]),1/tStep,1,usefulBaseline),
#                                                              bounds=([0, 0, 0.5, usefulBaseline - np.abs(usefulBaseline*0.01)],
#                                                                      [0.8, np.inf, 2, usefulBaseline + np.abs(usefulBaseline*0.01) ]))
            #commonBaseline.append(rep2[3])
        except Exception as err:
            print('\(O_o)/')
            #print(err)
            rep3,err3 = [np.NaN,np.NaN,np.NaN,np.NaN],[np.NaN,np.NaN,np.NaN,np.NaN] 
        
        if showPlot:
            plt.plot(x_out,expStretch(x_out,*rep3),'--k')
        
        allResults.append((rep3,err3))
            
    
    gamm = np.asarray( [allResults[i][0][1] for i in range(len(allResults))] )
    err_gamm = np.asarray([allResults[i][1][1] for i in range(len(allResults))] )

    betas = np.asarray( [allResults[i][0][2] for i in range(len(allResults))] )
    err_betas = np.asarray([allResults[i][1][2] for i in range(len(allResults))])

    contrasts = np.asarray( [allResults[i][0][0] for i in range(len(allResults))] )
    err_contrasts = np.asarray([allResults[i][1][0] for i in range(len(allResults))])

    baselines = np.asarray( [allResults[i][0][3] for i in range(len(allResults))] )
    err_baselines = np.asarray([allResults[i][1][3] for i in range(len(allResults))])
    if showPlot:
        plt.show()
    
    
    
    outDset={}
    outDset['relaxRate'] = gamm
    outDset['err_relaxRate'] = err_gamm
    outDset['exponent'] = betas
    outDset['err_exponent'] = err_betas
    outDset['contrast'] = contrasts
    outDset['err_contrast'] = err_contrasts
    outDset['baselines'] = baselines
    outDset['err_baselines'] = err_baselines
    outDset['lastTw'] = lastTw
    outDset['tw_in_frames'] = (tStep/2+np.arange(firstTw,lastTw,tStep))
                               
    return outDset 
    
    
    
    
#%%

def purifier(x,y,dx,dy):
    x=x[~np.isnan(y)]
    dx=dx[~np.isnan(y)]
    dy=dy[~np.isnan(y)]
    y=y[~np.isnan(y)]
    return x,y,dx,dy

def binandclean( tAx, g2, nbins = 50):
    x,y,dx,dy = bin_data_3(tAx[1:],(g2[1:]),
                               np.logspace(np.log10(tAx[1]),np.log10(tAx[-1]),nbins) )
    x,y,dx,dy = purifier(x,y,dx,dy)
    return x,y,dx,dy
    

def bin_data_3(x,y,kbins):
    # kbins are the bin edges
    digitized = np.digitize(x, kbins)
    # digitized=0 are all x<kbin[0] and digitized=len(kbins) are all x>kbin[-1]
    y_out = np.asarray([np.nanmean(y[digitized == i]) for i in range(1, len(kbins))])
    x_out = np.asarray([x[digitized == i].mean() for i in range(1, len(kbins))])
    dy_out = np.asarray([np.nanstd(y[digitized == i]) for i in range(1, len(kbins))])
    dx_out = np.asarray([x[digitized == i].std() for i in range(1, len(kbins))])
    y_out = y_out[np.where(~np.isnan(x_out))]
    dy_out = dy_out[np.where(~np.isnan(x_out))]
    dx_out = dx_out[np.where(~np.isnan(x_out))]
    x_out = x_out[np.where(~np.isnan(x_out))]
    
    return x_out, y_out, dx_out,dy_out

def expStretch(x,cont,Gamm,bet,baseline):
    return cont*np.abs(np.exp(-(Gamm*x)**bet))**2+baseline
    
