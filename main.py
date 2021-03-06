# -*- coding: utf-8 -*-
"""
Created on Mon Mar  11 14:07:46 2013

@author: Diego Castaneda
"""
#### MPI run?----
is_mpi = True
#### ------------

import numpy as np
import readpy
#from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
#import pylab as pl
import glob,shutil,os,subprocess,shlex
import fluxtable as f
import sys


if is_mpi==True:
	from mpi4py import MPI
	

### Stellar grid definition #####
nzphi = 400	# Phi zones
nzth = 200	# Theta zones
dth = 180./nzth
dphi = 360./nzphi
d = 14.68	# Distance to star [pc]
incl = [0.0,10.,20.,30.,40.,50.,60.,70.,80.,90.]	# Inclinations to compute [deg]
# Atmospheric Grid Specs:
if (len(sys.argv) > 4):
   dt_grid = float(sys.argv[4])	# Step in Teff [Teff]
else:
   dt_grid = 250. 
dg_grid = 1./3.	# Step in log(g)[cgs]
#Wavelength specs:
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


wavmax = 500. # [angstrom]

#Additional suffix for output? if yes, wavelength range will be pasted in filename
sffx = True
#Atmospheres directory
atmdir = "/home/castaned/atmospheres/"

which_os="linux"

if is_mpi==True:
################ MPI DEF ####################
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()

##############################################

gen_atm_files = True
only_interp_atm_files = False

# Definition of functions
dr_arr = []
phi = []
model = []
def modelread():
    global model, dr_arr, rotorc_dr, xi_dr
    model = readpy.genfromfort("model")
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
          
    dr_arr = np.empty(21,dtype=float)
    dr_arr[0] = 0.
    dr_arr[1:20] = (np.diff(model[:,1]))
    dr_arr[20] = 0.
    interp = modelinterpolate(model,dr_arr)
    return interp
    
def  modelinterpolate(model,dr_arr):
    global nzth,phi,nzphi, xi_dr
    xi = model[:,0]
    #xi_dr = np.linspace(13.5,166.5,len(dr_arr))
    xi_dr = np.linspace(9,171,len(dr_arr))
    y = []
    x1 = np.linspace(0,180,nzth+1)
    midx = x1[0:nzth]+np.diff(x1)/2.
    phi = np.linspace(0,360,nzphi)
    # spline order: 1 linear, 2 quadratic, 3 cubic ... 
    order = 1
    y.append(midx)
    # do inter/extrapolation
    for i in range(1,len(model[0,:])):
        #s = InterpolatedUnivariateSpline(xi, model[:,i], k=order)
        s = np.interp(midx,xi,model[:,i])
        #y.append(s(midx))
        y.append(s)
        
    s = InterpolatedUnivariateSpline(xi_dr,dr_arr, k=order)
    y.append(s(midx))
    y = np.array(y)
    y = np.transpose(y)
    return y

def xi(i,r1,theta,dr,phi):
    global dth,dr_arr,alambda
    dx1 = np.sin(np.deg2rad(i))
#    dy1 = 0.
    dz1 = np.cos(np.deg2rad(i))
    x1 = r1*np.sin(np.deg2rad(theta))*np.cos(np.deg2rad(phi))
    y1 = r1*np.sin(np.deg2rad(theta))*np.sin(np.deg2rad(phi))
    z1 = r1*np.cos(np.deg2rad(theta))
    pert = 0.1*r1
    alambda = dr/(r1*np.deg2rad(9))
    r3 = pert*np.tan(alambda)
    r2 = np.sqrt((r1+pert)**2+r3**2)
    psi = np.rad2deg(np.arcsin(r3/r2))
    th2 = theta-psi
    x2 = r2*np.sin(np.deg2rad(th2))*np.cos(np.deg2rad(phi))
    y2 = r2*np.sin(np.deg2rad(th2))*np.sin(np.deg2rad(phi))    
    z2 = r2*np.cos(np.deg2rad(th2))
    dx2 = x2-x1
    dy2 = y2-y1
    dz2 = z2-z1
    v2 = np.sqrt(dx2**2+dy2**2+dz2**2)
    cosang = (dx1*dx2+dz1*dz2)/v2
    return cosang
    
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
    global which_os
    for i in range(len(atmfiles)):
		tmpfl = glob.glob(atmdir+atmfiles[i]+".gz")
		print tmpfl
		os.system("cp "+tmpfl[0]+" "+os.getcwd())
		if(which_os=="windows"):
			subprocess.call(["bin/gzip","-d",atmfiles[i]+".gz"],shell=True)
		elif(which_os=="linux"):
			os.system("gunzip "+atmfiles[i]+".gz")
			#subprocess.call(["gunzip",atmfiles[i]+".gz"],shell=True)
    return
    
def interpatms(atmfiles,ilambda,flambda,dlambda):
    global which_os
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
###############################################################################
##########    ##    ##    ##    ## ##    ##    ################################
##########    ###  ###   #  #   ## ###   ##    ##### Calls to every ###########
##########    ## ## ##  ######  ## ## ## ##    ##### Function #################
##########    ##    ## ##    ## ## ##   ###    ################################
###############################################################################

## Generate the interpolated model of the integration grid:
interpmodel = modelread()
cossquiggle = []
darea = []

#ROTORC model dtheta:
model_dth = abs(model[0,0]-model[1,0])

#Geometric correction factor:
correct = np.sqrt(1.+interpmodel[:,-1]**2/(interpmodel[:,1]*np.deg2rad(model_dth))**2)

#Cossquiggle and darea calculation for every point in the grid:
for angle in phi:
	if is_mpi == True:
		xitemp = xi(incl[rank],interpmodel[:,1],interpmodel[:,0],interpmodel[:,-1],angle)
	else:
		xitemp = xi(incl,interpmodel[:,1],interpmodel[:,0],interpmodel[:,-1],angle)
	cossquiggle.append(xitemp)
	radius = interpmodel[:,1]*6.9598e10
	darea.append(correct*np.sin(np.deg2rad(interpmodel[:,0]))*dth*dphi*(np.pi*radius)**2/(180.**2))

#Convert lists to numpy.array:    
cossquiggle = np.array(cossquiggle)
darea = np.array(darea)

#Calculate the grid points in temp and grav for each point in the grid:
tandg = genatmgridvals(interpmodel,dt_grid,dg_grid)
#Generate the names of the atmospheres required
atmfiles = genfilenames(tandg)
#Generate the index correlation between atmosphere file and the temp-grav grid
atmkeys = getatmkeys(tandg,atmfiles)

if is_mpi == True:
	if rank==0:
		try:
			os.mkdir("tmp")
		except OSError:
			print "tmp exists"
		# Save auxiliary files for Fortran for this inclination
		np.savetxt("tmp/atmfilelist",atmfiles,fmt="%s")
		np.savetxt("tmp/darea",darea)
		np.savetxt("tmp/tandg",tandg)
		np.savetxt("tmp/atmkeys",atmkeys,fmt="%i")
		sendmsg = True
	else:
		sendmsg = None
	recvmsg = comm.bcast(sendmsg, root=0)
	np.savetxt("tmp/cossquiggle_i"+str(incl[rank]),cossquiggle)
else:
	try:
		os.mkdir("tmp")
	except OSError:
		print "tmp exists"
	# Save auxiliary files for Fortran for this inclination
	np.savetxt("tmp/atmfilelist",atmfiles,fmt="%s")
	np.savetxt("tmp/cossquiggle_i"+str(incl),cossquiggle)
	np.savetxt("tmp/darea",darea)
	np.savetxt("tmp/tandg",tandg)
	np.savetxt("tmp/atmkeys",atmkeys,fmt="%i")
##############################################################
if is_mpi == True:
	if rank==0:
		if gen_atm_files == False:
			if len(glob.glob("*.i"))!=len(atmfiles):
				f = open('check_error', 'w')
				f.close()
				exit("Inconsistent number of .i files")
			check = glob.glob("*.i")[0]
		else:
			if only_interp_atm_files == False:
				## Copy atmosphere .70 files
				copyatmfiles(atmdir,atmfiles)
			if len(glob.glob("*.70"))!=len(atmfiles):
				f = open('check_error', 'w')
				f.close()
				exit("Inconsistent number of .70 files")
			## Interpolate them
			interpatms(atmfiles,ilambda,flambda,dlambda)
		sendmsg = True
	else:
		sendmsg = None
	recvmsg = comm.bcast(sendmsg, root=0)

# Number of intensity angles in     
ntheta = 32
angles = np.empty(ntheta)
atmfname = "atmfilelist"
nfnames = len(atmfiles)
ext = str(incl[rank])
dist = d*3.08567758e18
nwav = int((flambda-ilambda)/dlambda)
wav = np.empty(nwav)
angles = np.empty(ntheta)
ftable = np.empty([nwav,ntheta,nfnames])
ftableout = f.fluxtable(ftable,wav,angles,atmfname,nfnames,ntheta,nwav)
ftable = ftableout[0]
wav = ftableout[1]
wavmax = wavmax/dlambda 
done = 0
if(int((flambda-ilambda)/dlambda)>wavmax):
    runs = int(((flambda-ilambda)/dlambda)/wavmax)
    modulo = ((flambda-ilambda)/dlambda)%wavmax
    if(modulo>0):
        po = 1
    else:
        po = 0
    for r in range(runs+po):
        if r<runs:
            nilam = r*wavmax*dlambda + ilambda
            nflam = nilam+wavmax*dlambda
        else:
            nilam = runs*wavmax*dlambda + ilambda
            nflam = nilam+modulo*dlambda
        nwavt = int((nflam-nilam)/dlambda)
        ul = r*wavmax+nwavt
        ll = r*wavmax
        wavt = wav[ll:ul]
        tempflux = np.empty([nwavt,2])
        #iout = np.arange(nwav*nzphi*nzth).reshape(nwav,nzphi,nzth)
        sed = f.interp(ftable[ll:ul,:,:],wavt,angles,ext,dist,nfnames,ntheta,nwavt,nzth,nzphi,incl)
        if(r<10):
            np.savetxt("outputflux_i"+str(incl[rank])+".r0"+str(r),sed)
        elif(r>=10):
            np.savetxt("outputflux_i"+str(incl[rank])+".r"+str(r),sed)
else:
    tempflux = np.empty([nwav,2])
    #iout = np.arange(nwav*nzphi*nzth).reshape(nwav,nzphi,nzth)
    sed = f.interp(ftable,wav,angles,ext,dist,nfnames,ntheta,nwav,nzth,nzphi,incl)
    np.savetxt("outputflux_i"+str(incl[rank]),sed)


done = 1
if rank!=0:
    comm.send(done,dest=0,tag=rank*10)
if rank==0:
    for r in range(size-1):
        recv = comm.recv(source=r+1,tag=(r+1)*10)
    if(int((flambda-ilambda)/dlambda)>wavmax):
        if(which_os=="linux"):
            for i in incl:
                cmd1 = "cat outputflux_i"+str(i)+".r*"+" > "+"outputflux_i"+str(i)+"_"+str(int(ilambda))+"-"+str(int(flambda))+"A"
                cmd2 = "rm outputflux_i"+str(i)+".r*"
                exit = subprocess.call(cmd1,shell=True)
                if exit==0:
                   subprocess.call(cmd2,shell=True)
                if glob.glob("outputflux_i"+str(i)+"_*")>1:
                    cmd1 = "cat outputflux_i"+str(i)+"_*"+" > "+"outputflux_i"+str(i)+".final"
                    subprocess.call(cmd1,shell=True)
"""
##################JASON TABLE##############################
f_r = interpolate.interp1d(interpmodel[:,0],interpmodel[:,1])
jangles = np.genfromtxt('Jtableangles.dat')
values = []
k=0
for i in range(len(interpmodel[:,0])):
    for j in range(len(phi)):
        #print np.round(interpmodel[i,0]),np.round(jangles[k,0]),np.round(phi[j]),np.round(jangles[k,1])
        if(np.round(interpmodel[i,0])==np.round(jangles[k,0]) and np.round(phi[j])==np.round(jangles[k,1])):
            values.append(np.array([interpmodel[i,0],phi[j],np.log10(iout[0,j,i]),darea[j,i],cossquiggle[j,i],correct[i]]))
            #print interpmodel[i,0],interpmodel[i,1]
            k+=1
            if (k >= len(jangles[:,0])):
                k=0

values = np.array(values)
np.savetxt('jTable.dat',values)
#########################################################33333
"""
#theta = np.linspace(0,180,nzth)
#pl.plot(theta,alambda)
#oldClic = readpy.genfromfort("outputflux",ncols=(0,1))
#pl.plot(fout[:,0],fout[:,1],"x")
#pl.plot(fout[:,0],fout[:,1]/oldClic[0:200,1])
#pl.plot(oldClic[:,0],oldClic[:,1])
#pl.show()
