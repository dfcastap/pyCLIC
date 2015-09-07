# -*- coding: utf-8 -*-
"""
Created on Thu May 15 13:40:03 2014

@author: castaned
"""
import numpy as np
import glob,os,subprocess,sys
from readmodel import modelread,modelinterpolate

atmdir = "/scratch2/castaned/atmospheres/"
### Stellar grid definition #####
nzphi = 400	# Phi zones
nzth = 200	# Theta zones
dth = 180./nzth
dphi = 360./nzphi

# Atmospheric Grid Specs:
if (len(sys.argv) > 1):
   dt_grid = float(sys.argv[1])	# Step in Teff [Teff]
else:
   dt_grid = 250. 
dg_grid = 1./3.	# Step in log(g)[cgs]

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

def copyatmfiles(atmdir,atmfiles):
    which_os = "linux"
    for i in range(len(atmfiles)):
		tmpfl = glob.glob(atmdir+"*/"+atmfiles[i]+".gz")
		print atmdir+"*/"+atmfiles[i]+".gz"
		os.system("cp "+tmpfl[0]+" "+os.getcwd())
		if(which_os=="windows"):
			subprocess.call(["bin/gzip","-d",atmfiles[i]+".gz"],shell=True)
		elif(which_os=="linux"):
			os.system("gunzip "+atmfiles[i]+".gz")
			#subprocess.call(["gunzip",atmfiles[i]+".gz"],shell=True)
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
np.savetxt("tmp/atmfilelist",atmfiles,fmt="%s")
    
copyatmfiles(atmdir,atmfiles)
