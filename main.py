# -*- coding: utf-8 -*-
"""
Created on Mon Mar  11 14:07:46 2013

___  _   _ ____ _    _ ____    
|__]  \_/  |    |    | |       
|      |   |___ |___ | |___    


@author: Diego Castaneda

Based on the Fortran code of Catherine Lovekin (2006ApJ...643..460L), 
it calculates the observed SED of a rotating stellar model at any inclination. 
It relies on PHOENIX atmospheric models (more specifically, the *.70 
intensity files).

"""
#### MPI run?----
is_mpi = False
#### ------------

import numpy as np

import glob,os,subprocess
import fluxtable as f
import sys
from readmodel import modelread,modelinterpolate


if is_mpi==True:
	from mpi4py import MPI
	
def xi(i,r1,theta,dr,phi):
    """
    Function that calculates xi (squiggle).
    (see eq. 8, Lovekin et. al. 2006 for details!)
    Parameters:
        i: inclination
        r1: radius
        theta: colatitude 
        dr: 
        phi: 

    Returns:
        Cos(xi) 
    """
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
    """
    Based on the atmospheric grid specification, it calculates
    the interpolation factors required to map the surface of the
    star with the given PHOENIX intensity files

    Parameters:
        model: stellar surface properties
        dt_grid: step in T of atmosphere grid
        dg_grid: step in log(g) of atmosphere grid
    Returns:
        array of interpolation weight factors
    """
    global factor_t
    tup = dt_grid*np.floor(model[:,2]/dt_grid)+dt_grid
    tdown = dt_grid*np.floor(model[:,2]/dt_grid)
    factor_t = (np.log10((model[:,2]))-np.log10((tdown)))/(np.log10(tup)-np.log10(tdown))
    gup = (dg_grid*1e5*np.floor(model[:,3]/(dg_grid))+dg_grid*1e5)/1e5
    gdown = (gup*1e5-dg_grid*1e5)/1e5
    factor_g = (model[:,3]-gdown)/(gup-gdown)
    return np.array(np.transpose([tup,tdown,gup,gdown,factor_t,factor_g]))
        

def genfilenames(tandg):
    """
    This function, based on the model, generate the filename list
    of PHOENIX atmosphere files required for pyCLIC to run

    Parameters:
        tandg: matrix of temperatures and log(g) needed for the run
    Returns:
        sorted list of PHOENIX filenames
    """
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
    """
    This function relates the PHOENIX file names to the interpolating
    factors.

    Parameters:
        tandg: matrix of temperatures and log(g) needed for the run
        atmfiles: sorted list of PHOENIX filenames
    Returns:
        Array to match the interpolation weight factors with the index
        of the atmosphere filename list
    """
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
    """
    Copy the atmospheric files to the cwd

    Parameters:
        atmdir: path to *.tar.gz PHOENIX *.70 files
        atmfiles: list of files to copy
    Return:
        Nothing...
    """
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
    
def copyatmfiles_i(atmdir,atmfiles):
    """
    Copy the (already interpolated *.70.i) atmospheric files to the cwd

    Parameters:
        atmdir: path to PHOENIX *.70.i files
        atmfiles: list of files to copy
    Return:
        Nothing...
    """
    global which_os
    
    for i in range(len(atmfiles)):
		tmpfl = glob.glob(atmdir+atmfiles[i]+".i")
		os.system("cp "+tmpfl[0]+" "+os.getcwd())
		#if(which_os=="windows"):
		#	subprocess.call(["bin/gzip","-d",atmfiles[i]+".gz"],shell=True)
		#elif(which_os=="linux"):
			#os.system("gunzip "+atmfiles[i]+".gz")
			#subprocess.call(["gunzip",atmfiles[i]+".gz"],shell=True)
    return
    
def interpatms(atmfiles,ilambda,flambda,dlambda):
    """
    Set up the interpolation input file. This uses Aaron Gillich
    PHOENIX files interpolation routine
    Parameters:
        atmfiles: list of files to interpolate
        ilambda: target initial wavelength
        flambda: target final wavelength
        dlambda: target wavelength step
    Return:
        Nothing...
    """
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
    
def doppler(vel,incl):
    """
    Calculates the Doppler factor at every point on the surface mesh
    Parameters:
        vel: rotational velocity
        incl: inclination
    Returns:
        Array of doppler factors
    """
    global phi,c,deltal
    #Speed of light
    c= 2.99792458e10
    dopp_arr = []
    for angle in phi:
        deltal = vel*10**5*np.sin(np.deg2rad(incl))*np.sin(np.deg2rad(angle))/c
        dopp_arr.append(deltal)
    dopp_arr = np.array(dopp_arr)
    return dopp_arr
    

def run_CLIC(model_name,incl,is_mpi,ilambda,flambda,dlambda,par,dt_grid=250,dg_grid = 1./3.,dopp = False):
    """
    Main routine. Defines most of the run. Should probably use an external parameter file...
    Currently one run per inclination (it's fairly fast though: ~5 minutes / 10000 wavelength points)
    If you're desperate you can do an MPI run: One inclination / core
    Parameters:
        model_name: model filename
        incl: array, target inclination(s) - len > 1 requires MPI
        is_mpi: Bool, is it mpi?
        ilambda: target initial wavelength
        flambda: target final wavelength
        dlambda: target wavelength step
        par: (for mode visibility) is it even or odd mode ---
        dt_grid: step in T of atmosphere grid
        dg_grid: step in log(g) of atmosphere grid
        dopp: Bool. include doppler effects? 
    Returns:
        Nothing, it writes out a 
    """
    global phi,c,deltal
    global dth,dr_arr,alambda
    global which_os
    global factor_t
    ### Stellar grid definition #####
    nzphi = 400	# Phi zones
    nzth = 200	# Theta zones
    dth = 180./nzth
    dphi = 360./nzphi
    d = 1.	# Distance to star [pc]
    #incl = [0.0,10.,20.,30.,45.,50.,60.,70.,80.,90.]	# Inclinations to compute [deg]
    # Atmospheric Grid Specs:
    
    #dt_grid = 250.
    #dg_grid = 1./3.	# Step in log(g)[cgs]
    """
    if (len(sys.argv) > 4):
       dt_grid = float(sys.argv[4])	# Step in Teff [Teff]
    else:
       dt_grid = 250.
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
    """
    #include doppler effects?
    if dopp==False:
        dopp = 0
    else:
        dopp = 1


    ofln="outputflux"
    if(dopp!=0):
        ofln="lineflux" 
        flambda=ilambda+99
        
        
    wavmax = 1000. # [angstrom]
    
    
    
    #Additional suffix for output? if yes, wavelength range will be pasted in filename
    sffx = True
    #Atmospheres directory
    atmdir = "/home/diego/Documents/pyCLIC/atms/"
    
    which_os="linux"
    
    ################ MPI DEF ####################
    if is_mpi==True:
    
    	comm = MPI.COMM_WORLD
    	rank = comm.Get_rank()
    	size = comm.Get_size()
    
    else:
        rank = 0
        size = 1
        incl = [incl[0]]
    ##############################################
    
    
   
    ###############################################################################
    ##########    ##    ##    ##    ## ##    ##    ################################
    ##########    ###  ###   #  #   ## ###   ##    ##### Calls to every ###########
    ##########    ## ## ##  ######  ## ## ## ##    ##### Function #################
    ##########    ##    ## ##    ## ## ##   ###    ################################
    ###############################################################################
    
    ## Generate the interpolated model of the integration grid:
    dr_arr = np.empty(21,dtype=float)
    phi = np.linspace(0,360.-360./nzphi,nzphi)
    model = modelread(model_name,dr_arr,par)
    interpmodel = modelinterpolate(model,dr_arr,nzth,nzphi)
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
    		xitemp = xi(incl[0],interpmodel[:,1],interpmodel[:,0],interpmodel[:,-1],angle)
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
    # Generate doppler array
    
    if dopp!=0:
        if is_mpi == True:
            doppler_arr = doppler(interpmodel[:,4],incl[rank])
        else:
            doppler_arr = doppler(interpmodel[:,4],incl[0])
    else:
        doppler_arr = np.zeros((len(cossquiggle[:,0]),len(cossquiggle[0,:])))
    
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
    	np.savetxt("tmp/doppler_arr_i"+str(incl[rank]),doppler_arr)
    else:
    	try:
    		os.mkdir("tmp")
    	except OSError:
    		print "tmp exists"
    	# Save auxiliary files for Fortran for this inclination
    	np.savetxt("tmp/atmfilelist",atmfiles,fmt="%s")
    	np.savetxt("tmp/cossquiggle_i"+str(incl[0]),cossquiggle)
    	np.savetxt("tmp/doppler_arr_i"+str(incl[0]),doppler_arr)
    	np.savetxt("tmp/darea",darea)
    	np.savetxt("tmp/tandg",tandg)
    	np.savetxt("tmp/atmkeys",atmkeys,fmt="%i")
    ##############################################################
    """
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
    """
    
    copyatmfiles_i(atmdir,atmfiles)
    # Number of intensity angles in     
    ntheta = 32
    angles = np.empty(ntheta)
    atmfname = "atmfilelist"
    nfnames = len(atmfiles)
    if is_mpi==True:
        ext = str(incl[rank])
    else:
        ext = str(incl[0])
        
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
    
    print (flambda-ilambda)/dlambda,wavmax 
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
            if dopp!=0:
                if is_mpi == True:
                    doppler_arr = doppler(interpmodel[:,4],incl[rank])
                else:
                    doppler_arr = doppler(interpmodel[:,4],incl[0])
                    sed = f.interp(ftable[ll:ul,:,:],wavt,angles,ext,dist,nfnames,ntheta,nwavt,nzth,nzphi,incl)
            else:
                doppler_arr = np.zeros((len(cossquiggle[:,0]),len(cossquiggle[0,:])))
                sed = f.interp_nodopp(ftable[ll:ul,:,:],wavt,angles,ext,dist,nfnames,ntheta,nwavt,nzth,nzphi,incl)
            
            
            if(r<10):
                np.savetxt(ofln+"_i"+str(incl[rank])+".r0"+str(r),sed)
            elif(r>=10):
                np.savetxt(ofln+"_i"+str(incl[rank])+".r"+str(r),sed)
                
    else:
        tempflux = np.empty([nwav,2])
        #iout = np.arange(nwav*nzphi*nzth).reshape(nwav,nzphi,nzth)
        if dopp!=0:
            if is_mpi == True:
               doppler_arr = doppler(interpmodel[:,4],incl[rank])
            else:
               doppler_arr = doppler(interpmodel[:,4],incl[0])
            sed = f.interp(ftable,wav,angles,ext,dist,nfnames,ntheta,nwav,nzth,nzphi,incl)
        else:
            doppler_arr = np.zeros((len(cossquiggle[:,0]),len(cossquiggle[0,:])))
            sed = f.interp_nodopp(ftable,wav,angles,ext,dist,nfnames,ntheta,nwav,nzth,nzphi,incl)
            #sed = f.interp(ftable,wav,angles,ext,dist,nfnames,ntheta,nwav,nzth,nzphi,incl)
        np.savetxt(ofln+"_i"+str(incl[rank]),sed)
    
    
    done = 1
    if rank!=0:
        comm.send(done,dest=0,tag=rank*10)
    if rank==0:
        for r in range(size-1):
            recv = comm.recv(source=r+1,tag=(r+1)*10)
        if(int((flambda-ilambda)/dlambda)>wavmax):
            if(which_os=="linux"):
                for i in incl:
                    cmd1 = "cat "+ofln+"_i"+str(i)+".r*"+" > "+ofln+"_i"+str(i)+"_"+str(int(ilambda))+"-"+str(int(flambda))+"A"
                    cmd2 = "rm "+ofln+"_i"+str(i)+".r*"
                    exit = subprocess.call(cmd1,shell=True)
                    if exit==0:
                       subprocess.call(cmd2,shell=True)
                    if glob.glob(ofln+"_i"+str(i)+"_*")>1:
                        cmd1 = "cat "+ofln+"_i"+str(i)+"_*"+" > "+ofln+"_i"+str(i)+".final"
                        subprocess.call(cmd1,shell=True)
    

    return
    
#run_CLIC("2p5Msun_V0_CLIC_s",[0.],is_mpi,3000.,12099.,2.0)
