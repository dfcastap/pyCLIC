# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 13:28:33 2015

@author: diego
"""
import numpy as np
import scipy.signal as sig
import pylab as plt
import pandas as pd
import scipy.integrate as scpi


"""
If we take the observed spectrum of Vega:
   Passband      Energy Flux       Photon "flux"
               (erg/sq.cm/sec)   (photons/sq.cm/sec)
  ---------------------------------------------------
      U          2.96e-06            550,000
      B          5.27e-06          1,170,000
      V          3.16e-06            866,000
      R          3.39e-06          1,100,000
      I          1.68e-06            675,000
  ---------------------------------------------------
"""

def filter(flux,resp_f):
    filtered = sig.convolve(flux,resp_f, mode='full')
    return filtered
    

flux = np.genfromtxt("outputflux_i0.0.final")
johnson = np.genfromtxt("johnson.dat",skip_header=2)
johsondf = pd.DataFrame(johnson)
johsondf.columns = ['wav_U','flux_U','wav_B','flux_B','wav_V','flux_V','wav_I','flux_I','wav_R','flux_R']

#johnson_U = np.transpose(np.array([johnson[:,0]*1e4,johnson[:,1]]))
johnson_U = np.array(johsondf[johsondf['wav_U']>0][['wav_U','flux_U']])
johnson_U[:,0] *= 1e4

#johnson_B = np.transpose(np.array([johnson[:,0]*1e4,johnson[:,1]]))
johnson_B = np.array(johsondf[johsondf['wav_B']>0][['wav_B','flux_B']])
johnson_B[:,0] *= 1e4

#johnson_V = np.transpose(np.array([johnson[:,0]*1e4,johnson[:,1]]))
johnson_V = np.array(johsondf[johsondf['wav_V']>0][['wav_V','flux_V']])
johnson_V[:,0] *= 1e4

#johnson_R = np.transpose(np.array([johnson[:,0]*1e4,johnson[:,1]]))
johnson_R = np.array(johsondf[johsondf['wav_R']>0][['wav_R','flux_R']])
johnson_R[:,0] *= 1e4

#johnson_I = np.transpose(np.array([johnson[:,0]*1e4,johnson[:,1]]))
johnson_I = np.array(johsondf[johsondf['wav_I']>0][['wav_I','flux_I']])
johnson_I[:,0] *= 1e4


johnson_U_i = np.interp(flux[:,0],johnson_U[:,0],johnson_U[:,1])
johnson_B_i = np.interp(flux[:,0],johnson_B[:,0],johnson_B[:,1])
johnson_V_i = np.interp(flux[:,0],johnson_V[:,0],johnson_V[:,1])
johnson_R_i = np.interp(flux[:,0],johnson_R[:,0],johnson_R[:,1])
johnson_I_i = np.interp(flux[:,0],johnson_I[:,0],johnson_I[:,1])

u = scpi.trapz(johnson_U_i*flux[:,1],flux[:,0])
b = scpi.trapz(johnson_B_i*flux[:,1],flux[:,0])
v = scpi.trapz(johnson_V_i*flux[:,1],flux[:,0])
r = scpi.trapz(johnson_R_i*flux[:,1],flux[:,0])
i = scpi.trapz(johnson_I_i*flux[:,1],flux[:,0])

"""
j=0
umag=0
for i in range(len(flux[:,0])):
     if (flux[i,0]>=johnson_U[j,0] and flux[i,0] < johnson_U[j+1,0]):
        fract = (flux[i,0]-johnson_U[j,0])/(johnson_U[j+1,0]-johnson_U[j,0])
        f = fract*johnson_U[j+1,1]+(1.-fract)*johnson_U[j,1]

        umag = umag + flux[i,1]*f*(flux[i+1,0]-flux[i-1,0])/2.

     if (flux[i,0]>johnson_U[j+1,0]): j = j+1
     if (j>11): j = 11

newout = np.empty((len(flux[:,0]),11))
newout[:,0] = 1.*flux[:,0]
newout[:,1:11] = np.transpose([1.*flux[:,1] for x in range(10)])
np.savetxt("outputflux",newout,fmt="%.8e")

plt.plot(flux[:,0],flux[:,1],color="0.5")
plt.plot(flux[:,0],johnson_U_i*flux[:,1])
plt.plot(flux[:,0],johnson_B_i*flux[:,1])
plt.plot(flux[:,0],johnson_V_i*flux[:,1])
plt.plot(flux[:,0],johnson_R_i*flux[:,1])
plt.plot(flux[:,0],johnson_I_i*flux[:,1])
"""