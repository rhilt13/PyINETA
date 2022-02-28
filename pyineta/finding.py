"""Functions for finding INETA networks.

This contains functions used by PyINETA to construct networks in the processed spectra.
Includes the following functions:
	* findClosestPoints - Uses KDTree to find close points.
	* mergeLevels - Merge points found in different iterations of peak picking.
	* horzAlign - Find horizontally aligned points.
	* getPairs - Get pairs.
	* buildNetwork - Build network.
	* listPairs - List all pairs of points that are part of a network.
"""

import itertools
import numpy as np
import networkx as nx
from collections import defaultdict
from scipy import spatial as spatial
import pyineta.clustering as clustering

def findClosestPoints (Pts,Qry):
	"""Uses KDTree to find close points.
	
	This function uses KDTree for nearest neighbor lookup to merge close points identified in multiple iterations.

	Args:
		Pts (ndarray): A 2D array of points.
		Qry (ndarray): A 1D array with a single point.
	
	Returns:
		tuple : A tuple with the distance and index of the closest point in Pts to Qry.
	"""

	PtsList = spatial.KDTree(Pts)
	QryRes=PtsList.query(Qry)
	return(QryRes)


def mergeLevels (Clist,levdist):
	"""Merge points found in different iterations of peak picking.

	This function is used to merge points found in multiple iterations of peak picking into a single array.
	Points that are closer to each other than a defined threshold are merged into a single point.

	Args:
		Clist (dict): A dict mapping the iteration number of peak picking to a list of cluster centers for points found in that iteration.
		levdist (float): Threshold for merging cluster centers from different iterations into a single point.
	
	Returns:
		ndarray : A 2D array of merged points. 
	"""

	mergedClistAll=np.vstack(list(Clist.values()))
	mergedClist=Clist[list(Clist.keys())[0]]
	for k, C in list(Clist.items()):
		for P in C:
			res=findClosestPoints(mergedClist,P)
			if (res[0]>float(levdist)):
				mergedClist=np.vstack([mergedClist, P])
			else:
				avgP=np.mean([mergedClist[res[1]], P],0)
				mergedClist = np.delete(mergedClist, (res[1]), axis=0)
				mergedClist=np.vstack([mergedClist, avgP])	
	mergedClist.view('i8,i8').sort(order=['f0'], axis=0) # Sort numerically based on 1st column
	print("Step3.1==> Merging levels: ",mergedClistAll.shape[0]," Points merged to ",mergedClist.shape[0])
	return(mergedClist)

def horzAlign (P,threshold1,threshold2,threshold3):
	"""Find horizontally aligned points.

	This function starts building the INETA networks by finding the horizontally aligned points, i.e. points with same double quantum ppm values as defined by the thresholds.
	To find these points, this function also filters points using basic INADEQUATE spectrum rules to only keeps points that satisfy these rules.

	Args:
		P (ndarray): A 2D array of final merged points.
		threshold1 (float): Threshold to split the merged points along the DQ axis.
		threshold2 (float): Threshold to check if the sum of 13C ppm values for two points equals their DQ ppm or not.
		threshold3 (float): Threshold to check of the 2 points are equidistant from the diagonal.

	Returns:
		dict : A dict mapping horizontally aligned peaks to their indices.
	"""

	ct=0
	Pd=dict()
	Pd[0]=P
	AlnY=clustering.splitY(Pd,threshold1,1)
	filtI=defaultdict(list)
	for i, lst in list(AlnY.items()):
		if (len(lst)>=2):
			for (a, b) in itertools.combinations(enumerate(lst), 2):
				meanY=(a[1][1]+b[1][1])/2
				sumX=a[1][0]+b[1][0]
				if (abs(meanY-sumX)<=float(threshold2)):
					D1=abs(a[1][0]-a[1][1]/2) # {diag:y=2x;mid-point:(y/2,y) so dist=x-y/2}
					D2=abs(b[1][0]-b[1][1]/2)
					if (abs(D1-D2)<=float(threshold3)):
						filtI[ct].append(a[1])
						filtI[ct].append(b[1])
						ct+=1
	return (filtI)

def getPairs(lst):
	"""Get pairs.

	Get pairs of points.

	Args:
		lst (tuple): Input points.

	Yields:
		tuple : Pairs of points.
	"""

	i = iter(lst)
	second = first = curr = next(i)
	for curr in i:
		yield first, curr
		first = curr
	yield curr, second

def buildNetwork (HorzPts,VertPts):	
	"""Build network.

	This function connects the horizontally and vertically connected points to build complete INETA networks.

	Args:
		HorzPts (dict): A dict mapping horizontally aligned peaks to their indices.
		VertPts (dict): A dict mapping vertically aligned peaks to their indices.

	Returns:
		list : list of lists with points belonging to a network.
	"""

	AllPts=[]
	for i,v in list(VertPts.items()):
		AllPts.append(v)
	for i,v in list(HorzPts.items()):
		AllPts.append(v)
	i=0
	for v in AllPts:
		i+=1
		g = nx.Graph()
	for lst in AllPts:
		lstt=tuple(map(tuple, lst))
		for edge in getPairs(lstt):
			g.add_edge(*edge)
	Net=list(nx.connected_components(g))
	j=0
	for i in Net:
		Net[j]=list(i)
		j+=1
	return (Net)

def listPairs (horzPts,vertPts):
	"""List all pairs of horizontally connected points that are part of a network.

	This function is used by PyINETA to generate a final list of all the points that are included in a network.

	Args:
		horzPts (dict): A dict mapping horizontally aligned peaks to their indices.
		vertPts (dict): A dict mapping vertically aligned peaks to their indices.

	Returns:
		list : A list of all horizontally connected pairs of points included in a network.
	"""

	allPairs=[]
	for i,v in list(vertPts.items()):
		if (len(v)>=2):
			for a, b in itertools.combinations(enumerate(v),2):
				allPairs.append(list((a[1],b[1])))
		else:
			continue

	for i,v in list(horzPts.items()):
		if (len(v)>=2):
			for a, b in itertools.combinations(enumerate(v),2):
				allPairs.append(list((a[1],b[1])))
		else:
			continue
	return (allPairs)