# -*- coding: utf-8 -*-
"""
Created on Thu May 15 13:09:18 2014

@author: castaned
"""
import numpy as np
import readpy
#from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
#import pylab as pl


dr_arr = []
phi = []
model = []
def modelread(model_name,dr_arr,par):
    model = readpy.genfromfort(model_name)
    if len(model[:,0])==10:
        temp1 = model[::-1]
        temp2 =  np.tile(0.,(20,5))
        for i in range(20):
            if i<10:
                temp2[i] = model[i]
            else:
                temp2[i,0] = i*(9.)+4.5
                temp2[i,1:5] = temp1[i-10,1:5]
        model=1.*temp2
    
    if par=="ODD":
        model = readpy.genfromfort(model_name)
    
    dr_arr[0] = 0.
    dr_arr[1:20] = (np.diff(model[:,1]))
    dr_arr[20] = 0.
    #interp = modelinterpolate(model,dr_arr,nzth,nzphi)
    return model
    
def  modelinterpolate(model,dr_arr,nzth,nzphi):
    xi = model[:,0]
    #xi_dr = np.linspace(13.5,166.5,len(dr_arr))
    xi_dr = np.linspace(9,171,len(dr_arr))
    y = []
    x1 = np.linspace(0,180,nzth+1)
    midx = x1[0:nzth]+np.diff(x1)/2.
    #phi = np.linspace(0,360,nzphi)
    # spline order: 1 linear, 2 quadratic, 3 cubic ... 
    order = 1
    y.append(midx)
    # do inter/extrapolation
    for i in range(1,len(model[0,:])):
        #s = InterpolatedUnivariateSpline(xi, model[:,i], k=order)
        s = np.interp(midx,xi,model[:,i])
        #y.append(s(midx))
        y.append(s)
        
    #s = InterpolatedUnivariateSpline(xi_dr,dr_arr, k=order)
    y.append(np.interp(midx,xi_dr,dr_arr))
    y = np.array(y)
    y = np.transpose(y)
    return y