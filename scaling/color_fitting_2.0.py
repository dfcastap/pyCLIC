# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 15:36:31 2012

@author: Diego
"""
import numpy as np
import readpy
from pylab import *
from scipy import interpolate


#which = "\M2-1975V229p3"
#wdir = "C:\Users\Diego\Documents\SaintMarys\StellarRotation\M2-1975V229p3_SUITE" +str(which)
jmagout = readpy.genfromfort("jmagout")
#outFile = readpy.genfromfort("outputfluxFiltered")


fig_width_pt = 1280.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'png',
          'axes.labelsize': 16,
          'text.fontsize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'figure.figsize': fig_size}
rcParams.update(params)

def bc(fv):
    cV = -12.999909395430052
    d=14.68
    lum = np.genfromtxt("luminosity.dat")
    Mbol = 4.74 -2.5*np.log10(lum)
    Mv = -2.5*np.log10(fv*(d/10.)**2)+cV
    return Mbol-Mv
    
def gint(target, lim1, lim2, comp1, comp2):
    i = lim1
    counter=0
    flag = 0
    comp = comp1
    while flag==0:
        if(i>target):
            comp1=comp
            lim1=i
            i = (lim2+lim1)/2.
            comp = (comp1+comp2)/2
                    
        else:
            comp2 = comp
            lim2=i
            i = (lim2+lim1)/2.
            comp = (comp1+comp2)/2
            
        counter+=1
        if(np.floor(i*100)==np.floor(target*100)):
            flag=1
            break
        
    return comp
    
bv_logg30 = np.genfromtxt("magnitudes/BV_logg30")
bv_logg33 = np.genfromtxt("magnitudes/BV_logg33")
bv_logg36 = np.genfromtxt("magnitudes/BV_logg36")
bv_logg40 = np.genfromtxt("magnitudes/BV_logg40")
bv_logg43 = np.genfromtxt("magnitudes/BV_logg43")

cBV = -0.73315739785601064

modelg = np.genfromtxt('../model',usecols=(3))
modelr = np.genfromtxt('../model',usecols=(1))

bv_model = np.array(np.zeros(10))
deduced_temps = np.array(np.zeros(10))

loggs = [4.333,4.0,3.666,3.333,3.0]

#bolcorr = np.zeros(10)
bolcorr = bc(jmagout[:,2])

np.savetxt("bc.dat",bolcorr)

for i in range(10):
    bv_model[i] = -2.5*np.log10(jmagout[i,1]/jmagout[i,2])-cBV
   
x = np.arange(0,100,10)
np.savetxt("B-V.dat",np.transpose([x,bv_model]))
for i in range(10):
    
    for j in loggs:
        if(j>=modelg[i]):
            lim1=j
            
        if(j<modelg[i]):
            lim2=j
            break
        
    print "for modelg="+str(modelg[i])+": "+str(lim1)+", "+str(lim2)
    
    if(np.floor(lim2*10)==30):
        comp1 = bv_logg33
        comp2 = bv_logg30
    elif(np.floor(lim2*10)==33):
        comp1 = bv_logg36
        comp2 = bv_logg33
    elif(np.floor(lim2*10)==36):
        comp1 = bv_logg40
        comp2 = bv_logg36
    elif(np.floor(lim2*10)==40):
        comp1 = bv_logg43
        comp2 = bv_logg40     
    
    l=0
    if(len(comp1[:,0])>=len(comp2[:,0])):
        bv_interpolated = np.transpose(np.array([np.zeros(len(comp2[:,0])),np.zeros(len(comp2[:,0]))]))
        for j in range(len(comp1[:,0])):
            for k in range(len(comp2[:,0])):
                if(comp1[j,0]==comp2[k,0]):
                    tmpcomp1 = comp1[j,1]
                    tmpcomp2 = comp2[k,1]
                    bv_interpolated[l,0] = comp2[k,0]
                    bv_interpolated[l,1] = gint(modelg[i],lim1,lim2,tmpcomp1,tmpcomp2)
                    l+=1
    else:
        bv_interpolated = np.transpose(np.array([np.zeros(len(comp1[:,0])),np.zeros(len(comp1[:,0]))]))
        for j in range(len(comp2[:,0])):
            for k in range(len(comp1[:,0])):
                if(comp1[k,0]==comp2[j,0]):
                    tmpcomp1 = comp1[k,1]
                    tmpcomp2 = comp2[j,1]
                    bv_interpolated[l,0] = comp1[k,0]
                    bv_interpolated[l,1] = gint(modelg[i],lim1,lim2,tmpcomp1,tmpcomp2)
                    l+=1
    
    l=0
    f = interpolate.interp1d(-1*bv_interpolated[:,1],bv_interpolated[:,0])        
    #xnew = np.arange(8000,11000,333)
    #xnew = np.arange(-0.55,-0.8,-0.05)
    ynew = f(-1*bv_model[i])
    #print f(-1*bv_model[i])
    #ynew = f(-1*xnew)
    deduced_temps[i] = ynew
    #plot(ynew,xnew,"k-")
    

"""
plot(bv_logg33[:,0],bv_logg33[:,1],"b-",)
plot(bv_logg36[:,0],bv_logg36[:,1],"g-",)
plot(bv_logg40[:,0],bv_logg40[:,1],"r-")
plot(bv_logg43[:,0],bv_logg43[:,1],"c-")
"""
modeltemps = genfromtxt('../model')
modeltemps = modeltemps[0:10]
#itemps = genfromtxt(wdir+'/array_temp_2.0.txt')
x = np.arange(0,100,10)
full_temp_arr = np.transpose(np.array([x,deduced_temps[:]]))
savetxt("full_temp_arr_2.0.dat",full_temp_arr,'%.4f')
plot(modeltemps[:,0],modeltemps[:,2],"b-",label="Latitude temps")
plot(x,deduced_temps[:],"g-",label="Inclination temps_BV")
#plot(itemps[:,0],itemps[:,1],"r-",label="Inclination temps_Vfit")
title("-- Rp/Re="+str(np.round(modelr[0]/modelr[9],3)))
xlabel("Inclination (degree)")
ylabel("Temperature (K)")
grid()
legend()
plt.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.07)
savefig('inclination_temps_overplot.png')
clf()
