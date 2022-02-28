"""Functions for reading spectra and peak picking.

This includes all the function definitions used for reading the NMR ft files and peak picking.
Includes the following functions:
	* readFt - Read ft files using nmrglue.
	* shifting - Shift the spectra as needed for proper referencing.
	* frange - Generate a range of floats.
	* pick - Peak picking function.
"""

import math
import numpy as np
import nmrglue as ng

def readFt (ftfile):
	"""Read ft files using nmrglue.

	Args:
		ftfile (str): ft filename.

	Returns:
		ndarray : A numpy array with the intensities.
		1D-array : A 1D array with the 13C ppm values.
		1D-array : A 1D array with the double quantum ppm values.
	"""

	ft_dic,ft_data = ng.pipe.read(ftfile)
	In=ft_data.transpose()

	DQ = ng.pipe.make_uc(ft_dic, ft_data, 0)
	Ys=DQ.ppm_scale()

	CS = ng.pipe.make_uc(ft_dic, ft_data, 1)
	Xs=CS.ppm_scale()
	return(In,Xs,Ys)

def shifting (In,pX,fullX,fullY,direction):
	"""Shift the spectra as needed for proper referencing.

	This function can be used to shift all the peaks along both the 13C and the DQ axis to reference them properly.
	Shift occurs with user provided units along the 13C axis and by double the units along DQ axis.
	The empty cells are zero filled to complete the output array.

	Args:
		In (ndarray): Array with the peak intensities.
		pX (int): Units to shift all the peaks (in points; 200 ppm ~ 4000 points).
		fullX (int): Length of the 13C ppm array (4096 for a typical INADEQUATE spectrum).
		fullY (int): Length of the DQ ppm array (8192 for a typical INADEQUATE spectrum).
		direction (str): 'pos' indicates shifting the peaks towards a higher ppm value while 'neg' indicates shifting peaks towards a lower ppm value.
	
	Returns:
		ndarray : Intensity array after shifting of peaks.
	"""
	ratio=fullY/fullX
	pY=math.ceil(pX*ratio)
	if direction.lower() == "pos":      # For increasing ppm axes (eg: peak at 35 ppm is now at 40 ppm)
		padX=np.zeros((pX,fullY))   # fullY= 8192 for INADEQUATE
		padY=np.zeros((fullX,pY))   # fullX=4096 for INADEQUATE
		In=np.hstack((In,padY))
		In=In[:,pY:]
		In=np.concatenate((In,padX))
		In=In[pX:,:]
	elif direction.lower() == "neg":        # For decreasing ppm axes (eg: peak at 40 ppm is now at 35 ppm)
		padX=np.zeros((pX,fullY))   # fullY= 8192 for INADEQUATE
		padY=np.zeros((fullX,pY))   # fullX=4096 for INADEQUATE
		In=np.hstack((padY,In)) 
		In=In[:,:-pY]
		In=np.concatenate((padX,In))
		In=In[:-pX,:]
	return(In)

def frange (start,end,parts):
	"""Generate a range of floats.

	Function to generate a list of floats from a given value to an end value.
	Generates list in reverse order if start is a larger value than end.

	Args:
		start (float): Starting value for the range list.
		end (float): End value for the range list.
		parts (int): Indicates how many parts to split the range of start-end.

	Returns:
		ndarray : A 1D array of floating point values ranging from start to end of size parts.
	"""

	duration=abs(end-start)
	part_duration = duration / (parts-1)
	return [start+(i * part_duration) for i in range(parts-1,-1,-1)]

def pick (In,xppm,yppm,PPmin,PPmax,steps):
	"""Peak picking function.

	This function iterates over a provided range of intensities and collects all cells with intensity values that fall within that range.

	Args:
		In (ndarray): Array with the peak intensities.
		xppm (ndarray): A 1D array with the 13C ppm values.
		yppm (ndarray): A 1D array with the double quantum ppm values.
		PPmin (float): Minimum intensity value to be considered a peak.
		PPmax (float): Maximum intensity value to be considered a peak.
		steps (int): Number of iterations to find peaks within the PPmin to PPmax range.

	Returns:
		dict : A dict mapping an array of points (x,y) to the iteration number it was found in.
		dict : A dict mapping an array of x axis values (13C ppm values) to the iteration number it was found in.
		dict : A dict mapping an array of y axis values (DQ ppm values) to the iteration number it was found in.
	"""

	points=np.copy(In)
	j=0
	X={}
	Y={}
	Pts={}
	for i in frange(PPmin,PPmax,steps):
		(Xind, Yind) = np.where(abs(points) > int(i))
		points[abs(points) > i] =0    # This removes hit peaks from the previous iteration
		selX=np.around(xppm[Xind],decimals=2)
		selY=np.around(yppm[Yind],decimals=2)
		selX=xppm[Xind]
		selY=yppm[Yind]
		X[j]=selX
		Y[j]=selY
		Pts[j]=list(zip(selX,selY))
		j=j+1
	return (Pts,X,Y)