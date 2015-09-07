# -*- coding: utf-8 -*-
"""
Created on Thu May 15 14:02:08 2014

@author: castaned
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 15 13:40:03 2014

@author: castaned
"""
import numpy as np
import glob,os,subprocess,sys
from readmodel import modelread,modelinterpolate

atmdir = "/home/castaned/atmospheres/"
### Stellar grid definition #####
nzphi = 400	# Phi zones
nzth = 200	# Theta zones
dth = 180./nzth
dphi = 360./nzphi

# Atmospheric Grid Specs:
if (len(sys.argv) > 4):
   dt_grid = float(sys.argv[4])	# Step in Teff [Teff]
else:
   dt_grid = 250. 
dg_grid = 1./3.	# Step in log(g)[cgs]

if (len(sys.argv) > 1):
	dlambda = float(sys.argv[3])
else:
	dlambda = 1.0 # [angstrom]

if (len(sys.argv) > 1):
	ilambda = float(sys.argv[1])
	flambda = float(sys.argv[2])
else:
	ilambda = 8000. # [angstrom]
	flambda = 19999. # [angstrom]

#is there a doppler flag?
dopp = 0
if (len(sys.argv) > 5):
    dopp = float(sys.argv[5])


if(dopp!=0):
    flambda=ilambda+99
    



def genatmgridvals(model,dt_grid,dg_grid):
    global factor_t
    tup = dt_grid*np.floor(model[:,2]/dt_grid)+dt_grid
    tdown = dt_grid*np.floor(model[:,2]/dt_grid)
    factor_t = (np.log10((model[:,2]))-np.log10((tdown)))/(np.log10(tup)-np.log10(tdown))
    gup = (dg_grid*1e5*np.floor(model[:,3]/(dg_grid))+dg_grid*1e5)/1e5
    gdown = (gup*1e5-dg_grid*1e5)/1e5
    factor_g = (model[:,3]-gdown)/(gup-gdown)
    return np.array(np.transpose([tup,tdown,gup,gdown,factor_t,factor_g]))
    
def genfilenames(tandg):
    fnames = []
    for i in range(len(tandg)):
        if int(tandg[i,2]*10) == 39: tandg[i,2] = 4.0
        if int(tandg[i,3]*10) == 39: tandg[i,3] = 4.0
        if int(tandg[i,2]*10) == 29: tandg[i,2] = 3.0
        if int(tandg[i,3]*10) == 29: tandg[i,3] = 3.0
        fnames.append("nlte"+str(int(tandg[i,0]))+"K"+str(int(tandg[i,2]*10))+".70")
        fnames.append("nlte"+str(int(tandg[i,0]))+"K"+str(int(tandg[i,3]*10))+".70")
        fnames.append("nlte"+str(int(tandg[i,1]))+"K"+str(int(tandg[i,2]*10))+".70")
        fnames.append("nlte"+str(int(tandg[i,1]))+"K"+str(int(tandg[i,3]*10))+".70")
    return np.sort(list(set(fnames)))
    
def getatmkeys(tandg,atmfiles):
    keys = np.zeros((len(tandg[:,0]),4))
    atmvals = []
    for s in atmfiles:
            atmvals.append([float(s.split('K')[0].split('e')[1]),float(s.split('K')[1].split('.')[0])])
    for i in range(len(tandg[:,0])):
        for s in range(len(atmvals[:])):
            if(tandg[i,0]==atmvals[s][0] and np.floor(tandg[i,2]*10.)==atmvals[s][1]): keys[i,0]=s
            if(tandg[i,1]==atmvals[s][0] and np.floor(tandg[i,2]*10.)==atmvals[s][1]): keys[i,1]=s
            if(tandg[i,0]==atmvals[s][0] and np.floor(tandg[i,3]*10.)==atmvals[s][1]): keys[i,2]=s
            if(tandg[i,1]==atmvals[s][0] and np.floor(tandg[i,3]*10.)==atmvals[s][1]): keys[i,3]=s
                

    return np.array(keys)

def interpatms(atmfiles,ilambda,flambda,dlambda):
    which_os = "linux"
    f = open('wavelength.nml','w')
    f.write("$lambdas \n")
    f.write("\n")
    f.write(" blambda = "+str(ilambda)+"\n")
    f.write(" tlambda = "+str(flambda)+"\n")
    f.write(" x = "+str(dlambda)+"d0"+" \n")
    f.write("\n")
    f.write("$end \n")
    f.close()
    if(which_os=="windows"):
        subprocess.call(["wavelengthinterp.exe"])
    elif(which_os=="linux"):
        subprocess.call(["./wavelengthinterplog.exe"])
    #if len(glob.glob("*.i"))==len(atmfiles):
        #toremove = glob.glob("*.70")
        #for i in range(len(toremove)):
            #os.remove(toremove[i])
    return
    
## Generate the interpolated model of the integration grid:
dr_arr = np.empty(21,dtype=float)
phi = np.linspace(0,360,nzphi)
model = modelread(dr_arr)
interpmodel = modelinterpolate(model,dr_arr,nzth,nzphi)
#Calculate the grid points in temp and grav for each point in the grid:
tandg = genatmgridvals(interpmodel,dt_grid,dg_grid)
#Generate the names of the atmospheres required
atmfiles = genfilenames(tandg)
    
#copyatmfiles(atmdir,atmfiles)
interpatms(atmfiles,ilambda,flambda,dlambda)