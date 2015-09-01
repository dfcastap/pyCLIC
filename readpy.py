import numpy as np

def genfromfort(flname, ncols = (0,0)):
	
    if ncols[1]==0:
		array = np.genfromtxt(flname,dtype=('string'))
    else:
		array = np.genfromtxt(flname,dtype=('string'),usecols=ncols)

    for i in range(len(array[:,0])):
			for j in range(len(array[0,:])):
				array[i,j]=array[i,j].replace("D","E")

    array = array.astype(np.float64)

    return array
