from scipy import spatial as spatial
import numpy as np
import pyineta.filtering as filtering
from collections import defaultdict
import itertools
import networkx as nx

def findClosestPoints (Pts,Qry):
	PtsList = spatial.KDTree(Pts)
	QryRes=PtsList.query(Qry)
	return(QryRes)

def mergeLevels (Clist,levdist):
	mergedClistAll=np.vstack(list(Clist.values()))
	mergedClist=Clist[0]
	for k, C in list(Clist.items()):
		for P in C:
			res=findClosestPoints(mergedClist,P)
			# print(res,P,mergedClist[res[1]])
			if (res[0]>float(levdist)):
				mergedClist=np.vstack([mergedClist, P])
				# Clistmerge2[q]=
			else:
				closestP=np.vstack([mergedClist[res[1]], P])
				# print "ClP=>",closestP, res[1], mergedClist[res[1]], mergedClist
				avgP=np.mean([mergedClist[res[1]], P],0)
				# print mergedClist[res[1]], P, avgP
				mergedClist = np.delete(mergedClist, (res[1]), axis=0)
				mergedClist=np.vstack([mergedClist, avgP])	
	mergedClist.view('i8,i8').sort(order=['f0'], axis=0) # Sort numerically based on 1st column
	print("Step3.1==> Merging levels: ",mergedClistAll.shape[0]," Points merged to ",mergedClist.shape[0])
	return(mergedClist)


def horzAlign (P,threshold1,threshold2,threshold3):
	ct=0
	Pd=dict()
	Pd[0]=P
	AlnY=filtering.sortY(Pd,threshold1,1)
	# AlnY=Sort(P,threshold1,1)
	filtI=defaultdict(list)
	# print "AlnY==> ",AlnY
	for i, lst in list(AlnY.items()):
		# print lst,"==",len(lst)
		if (len(lst)>=2):
			for (a, b) in itertools.combinations(enumerate(lst), 2):
				# print a,b,a[1],b[1]
				meanY=(a[1][1]+b[1][1])/2
				sumX=a[1][0]+b[1][0]
				# print "Y=",a[1][1],b[1][1],"X=",a[1][0],b[1][0],"meanY=",meanY,"sumX=",sumX,"//",abs(meanY-sumX),"//",threshold2
				if (abs(meanY-sumX)<=float(threshold2)):
					D1=abs(a[1][0]-a[1][1]/2) # {diag:y=2x;mid-point:(y/2,y) so dist=x-y/2}
					D2=abs(b[1][0]-b[1][1]/2)
					if (abs(D1-D2)<=float(threshold3)):
						filtI[ct].append(a[1])
						filtI[ct].append(b[1])
						# print filtI[ct]
						ct+=1
						# print "===>Passed Diag", D1,D2,threshold3,"Y=",a[1][1],b[1][1],"X=",a[1][0],b[1][0]
					# else:
						# print "#####Did not pass Diag", D1,D2,threshold3,"Y=",a[1][1],b[1][1],"X=",a[1][0],b[1][0]
	return (filtI)

def getPairs(lst):
	i = iter(lst)
	second = first = curr = next(i)
	for curr in i:
		yield first, curr
		first = curr
	yield curr, second

def buildNetwork (HorzPts,VertPts):	# Make network from selected peaks
	# Input filtered set of Xpeak points from above function
	# print HorzPts
	# print VertPts
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
			# print type(edge)
			# print type(lstt)
			g.add_edge(*edge)
	Net=list(nx.connected_components(g))
	j=0
	Netf=Net
	# print "NetBefore=",Net
	for i in Net:
		Net[j]=list(i)
		# for k in list(i):
		# 	print k, VertPts
			# if k in VertPts:
		j+=1
	# print "NetAfter=",Net
	return (Net)

def listPairs (horzPts,vertPts):
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