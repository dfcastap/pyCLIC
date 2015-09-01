# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:31:36 2013

@author: Diego
"""

import numpy as np
import glob
import scipy
import scipy.signal
import scipy.integrate as integrate

incl = [0.0,10.,20.,30.,40.,50.,60.,70.,80.,90.]
sampling = 0.02
boxsize = 50

def boxfilter(flux,boxsize,dlam):
    global window
    window = scipy.signal.get_window('boxcar', boxsize/dlam)
    window = window/(boxsize/dlam)
    filtered = scipy.signal.convolve(flux,window, mode='same')
    return filtered

def luminosity(flux):
    lum = np.zeros(len(flux[0,:])-1)
    for i in range(len(flux[0,:])-1):
        lum[i]=integrate.trapz(flux[:,i+1],flux[:,0])
    return lum
    
def check_for_pieces(incl):
    fluxcontainer = []
    for i in range(len(incl)):
        catstr1="np.concatenate(("
        catstr2="np.concatenate(("
        fluxtempcontainer = []
        temp = glob.glob("outputflux_i"+str(incl[i])+"*")
        if len(temp)>0:
            for j in range(len(temp)):
                another = np.genfromtxt(temp[j])
                #plot(another[:,0],another[:,1])
                fluxtempcontainer.append(another)
                catstr1+="fluxtempcontainer["+str(j)+"][:,0]"
                catstr2+="fluxtempcontainer["+str(j)+"][:,1]"
                if j<len(temp)-1:
                    catstr1+=","
                    catstr2+=","
            
            catstr1+="))"
            catstr2+="))"
            wav = eval(catstr1)
            flux= eval(catstr2)
            fluxcontainer.append([wav,flux])
    return fluxcontainer

def resample(newx,flux):
    resampled=np.zeros((len(newx[:]),len(flux[0,:])))
    resampled[:,0]=newx
    for i in range(len(flux[0,:])-1):
        resampled[:,i+1]=np.interp(newx,flux[:,0],flux[:,i+1])
    return resampled
    
def filteroutput(flux,boxsize,dlam):
    filtered = np.zeros_like(flux)
    filtered[:,0] = flux[:,0]
    for i in range(len(flux[0,:])-1):
        filtered[:,i+1] = boxfilter(flux[:,i+1],boxsize,dlam)
    return filtered
    
def build_outputflux(nincl):
    files = glob.glob("outputflux_i*.final")
    tempout = []
    for i in range(nincl):
        tempout.append(np.genfromtxt("outputflux_i"+str(float(i*10))+".final"))
    return tempout

d=14.68
nincl = len(incl)
#tempout = check_for_pieces(incl)
tempout = build_outputflux(nincl)
newout=np.zeros((len(tempout[0][:,0]),nincl+1))
if nincl>1:
    newout[:,0]=tempout[0][:,0]
    newout[:,1]=tempout[0][:,1]
    for i in range(nincl-1):
        newout[:,i+2]=(tempout[i+1][:,1])
        
lum = luminosity(newout)*(4*np.pi*(d*3.08567758e18)**2)/(3.9e33)

###############################################FILTER
#resampled_x = np.arange(1000.,20000.,sampling)
#resampledSED = resample(resampled_x,newout)
#filteredSED=filteroutput(resampledSED,boxsize,sampling)
######################################################

np.savetxt("scaling/luminosity.dat",np.transpose(lum),fmt="%.4f")
np.savetxt("scaling/outputflux",newout,fmt="%.8e")
#np.savetxt("scaling/outputfluxFiltered",filteredSED,fmt="%.8e")

