"""Functions for clustering and centering a list of peaks.

This includes all the functions used by PyINETA to cluster the list of peaks in a spectra to find a peak center.
Includes the following functions:
	* gather - Gather points along a given axis.
	* splitY - Split groups of points based on threshold.
	* centerMass - Find the center of mass for clustered points.
"""

import numpy as np
from collections import defaultdict

def gather (Points, thres, ax):
	"""Gather points along a given axis.

	This function clusters points within the threshold of the first point in the cluster.
	
	Args:
		Points (list): A list of points (x,y).
		thres (float): Points found within this threshold are grouped into a single cluster. 
		ax (int): '0' indicates gathering along the 13C axis while '1' indicates gathering along the DQ axis.
	
	Returns:
		dict : A dict mapping lists of clustered points to their respective cluster numbers.
	"""

	sets = defaultdict(list)
	sPts = sorted(Points, key=lambda tup: tup[ax]) # ax=0 for X, 1 for Y
	i=0
	sets[i].append(sPts.pop(0))
	for elem in sPts:
		minval=min(sets[i], key = lambda t: t[0])
		if (abs(minval[ax]-elem[ax]) < float(thres)):
			sets[i].append(elem)
		else:
			i=i+1
			sets[i].append(elem)
	return(sets)

def splitY (Points, thres, ax):
	"""Split groups of points based on threshold.

	This function continuously adds points to the same cluster as long as that point is within the threshold distance of the closest point in the cluster.
	This function is used by PyINETA to split the points gathered using the gather function.

	Args:
		Points (dict): A dict mapping lists of clustered points to their cluster number.
		thres (float): A distance threshold for splitting points into a different cluster.
		ax (int): '0' indicates splitting along the 13C axis while '1' indicates splitting along the DQ axis.
	
	Returns:
		dict : A dict mapping lists of clustered points to their reorganized cluster numbers after splitting.
	"""

	sets = defaultdict(list)
	i=0
	for k, Group in list(Points.items()):
		sPts = sorted(Group, key=lambda tup: tup[ax])	
		sets[i].append(sPts[0])
		for first, second in zip(sPts, sPts[1:]):
			if (abs(first[1]-second[1]) < float(thres)):
				sets[i].append(second)
			else:
				i=i+1
				sets[i].append(second)
		i=i+1
	return(sets)

def centerMass (Clusters,metric):
	"""Find the center of mass for clustered points.

	Args:
		Clusters (dict): A dict mapping list of points to their cluster number.
		metric (str): use "mean" or "median" as the method to find the center of mass for the clusters.

	Returns:
		ndarray : An array of center of masses (peak centers) for all the clusters.
	"""

	Center = []
	C=[]
	for k, clust in list(Clusters.items()):
		if (metric == "mean"):
			C1=sum(v[0] for v in clust) / float(len(clust))
			C2=sum(v[1] for v in clust) / float(len(clust))
		elif (metric =="median"):
			C1=np.median([v[0] for v in clust])
			C2=np.median([v[1] for v in clust])
		C=(C1,C2)
		Center.append(C)
		C=[]
	return(np.asarray(Center))