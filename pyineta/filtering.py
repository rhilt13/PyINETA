from collections import defaultdict
import numpy as np

def sortX (Points, thres, ax):
	sets = defaultdict(list)
	sPts = sorted(Points, key=lambda tup: tup[ax]) # ax=0 for X, 1 for Y
	# print sPts
	i=0
	sets[i].append(sPts.pop(0))
	# print sets[i]
	for elem in sPts:
		minval=min(sets[i], key = lambda t: t[0])
		# minval=min(sets[i], key = lambda t: t[ax])
		# print "MinVal for sets",i,"=>",minval,"//",elem
		# print minval[ax],elem[ax]
		if (abs(minval[ax]-elem[ax]) < float(thres)):
			sets[i].append(elem)
		else:
			i=i+1
			sets[i].append(elem)
	# print sets
	return(sets)

def sortY (Points, thres, ax):
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
	Center = []
	C=[]
	for k, clust in list(Clusters.items()):
		if (metric == "mean"):
			C1=sum(v[0] for v in clust) / float(len(clust))
			C2=sum(v[1] for v in clust) / float(len(clust))
		elif (metric =="median"):
			C1=np.median([v[0] for v in clust])
			C2=np.median([v[1] for v in clust])
		# print C1, C2
		C=(C1,C2)
		Center.append(C)
		C=[]
	return(np.asarray(Center))