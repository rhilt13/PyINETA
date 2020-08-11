#! /usr/bin/python

## Author:	Rahil Taujale
## Version:	1.11
## Date:	10/27/2017

## Description:
# Set of functions used by INETA

## Version:
# 0.2 = Changed the split function to extend minimum and maximum 
# 1.0 = Named _v1.py : Finalized
# 1.1 = Changed xrange and yrange of plots to 0-200 and 0-400 respectively
# 1.11 = Added major and minor gridlines

## To-Do
# Change HorzAlign and 3findNetwork functions
# Remove Xs, Ys
# Change variable names to sth appropriate for all functions

import itertools
from collections import defaultdict
import numpy as np 
import pylab as pl
import matplotlib.pyplot as plt
import networkx # standard import is "import networkx as nx"
from scipy import spatial as sp # just use "spatial", sp is not normal
# from statistics import mean, median
float_formatter = lambda x: "%.3f" % x
np.set_printoptions(formatter={'float_kind':float_formatter})

# def plotSingle (X,Y,Xlim,Ylim,label,imgname):
# 	if (type(label)==float):
# 		title="{:.2e}".format(label[i])
# 	else:
# 		title=label
# 	fig=plt.figure(figsize=(30, 30)) # This increases resolution
# 	ax = fig.add_subplot(111)
# 	ax.plot(X,Y,'k.',)
# 	ax.set_xlim(Xlim)
# 	ax.set_ylim(Ylim)
# 	ax.set_ylabel('Double Quantum', fontsize=30)
# 	ax.set_xlabel('13C', fontsize=30)
# 	ax.set_title(title, fontsize=30)
# 	ax.tick_params(axis='both', which='major', labelsize=25)
# 	ax.invert_xaxis()
# 	ax.invert_yaxis()
# 	plt.tight_layout(pad=5, w_pad=0.5, h_pad=1.0)
# 	plt.minorticks_on()
# 	plt.grid(b=True, which='major', color='0.45', linestyle=':', linewidth=0.25)
# 	plt.grid(b=True, which='minor', color='0.75', linestyle=':', linewidth=0.2)
# 	fig.savefig(imgname,format='eps',dpi=300)
# 	plt.close(fig)
# 	return(fig, ax)

# def plotFigSep (X,Y,Xlim,Ylim,label,imgname):
# 	r=np.ceil(float(len(X))/2)
# 	c=2
# 	l=c*30
# 	w=r*30
# 	fig, axs = plt.subplots(nrows=int(r), ncols=c, figsize=(l,w))
# 	for i, ax in enumerate(fig.axes):
# 		if i in X:
# 			title="{:.2e}".format(label[i])
# 			ax.plot(X[i],Y[i],'k.')
# 			ax.set_xlim(Xlim)
# 			ax.set_ylim(Ylim)
# 			ax.set_ylabel('Double Quantum', fontsize=60)
# 			ax.set_xlabel('13C', fontsize=60)
# 			ax.set_title(title, fontsize=60)
# 			ax.tick_params(axis='both', which='major', labelsize=40)
# 			ax.invert_xaxis()
# 			ax.invert_yaxis()
# 			plt.minorticks_on()
# 			plt.grid(b=True, which='major', color='0.45', linestyle=':', linewidth=0.25)
# 			plt.grid(b=True, which='minor', color='0.75', linestyle=':', linewidth=0.2)
# 	plt.tight_layout(pad=10, w_pad=10, h_pad=10)
# 	if (len(X)%2!=0):
# 		axs[-1,-1].axis('off')
# 	fig.savefig(imgname,format='eps',dpi=300)
# 	return(fig,axs)

# def PlotAllNet(HorzPts,VertPts,ax,rad):
# 	for i,v in list(HorzPts.items()):
# 		for q in v:
# 			circ=plt.Circle(q,radius=rad,color='b',fill=False,linestyle='dotted',linewidth=0.8)
# 			ax.add_patch(circ)
# 		ax.plot([v[0][0], v[1][0]], [v[0][1], v[1][1]], color='#f72d1c', linestyle='--', linewidth=2)
# 	for i,v in list(VertPts.items()):
# 		for a, b in itertools.combinations(enumerate(v),2):
# 			ax.plot([a[1][0], b[1][0]], [a[1][1], b[1][1]], color='#f72d1c', linestyle='--', linewidth=2)
# 	return (ax)

def getPairs(lst):
	i = iter(lst)
	second = first = curr = next(i)
	for curr in i:
		yield first, curr
		first = curr
	yield curr, second

# def peakPick (points,xppm,yppm,PPmin,PPmax,steps):
# 	j=0
# 	X={}
# 	Y={}
# 	P={}
# 	Pts={}
# 	for i in pl.frange(int(PPmax),int(PPmin),npts=steps):
# 		# print i
# 		(Xind, Yind) = np.where(abs(points) > i)
# 		points[abs(points) > i] =0          ## This line removes hit peaks from the previous iteration (do we do this??)
# 		selX=np.around(xppm[Xind],decimals=2)
# 		selY=np.around(yppm[Yind],decimals=2)
# 		selX=xppm[Xind]
# 		selY=yppm[Yind]
# 		selP=np.vstack((selX,selY))
# 		P[j]=selP
# 		X[j]=selX
# 		Y[j]=selY
# 		Pts[j]=list(zip(selX,selY))
# 		j=j+1
# 	return (Pts,X,Y,P)

# def frange(start,end,parts):
#     duration=abs(end-start)
#     part_duration = duration / (parts-1)
#     return [start+(i * part_duration) for i in range(parts-1,-1,-1)]

# def peakPick (points,xppm,yppm,PPmin,PPmax,steps):
#         j=0
#         X={}
#         Y={}
#         P={}
#         Pts={}
#         for i in frange(PPmin,PPmax,steps):
#                 (Xind, Yind) = np.where(abs(points) > int(i))
#                 points[abs(points) > i] =0          ## This line removes hit peaks from the previous iteration (do we do this??)
#                 selX=np.around(xppm[Xind],decimals=2)
#                 selY=np.around(yppm[Yind],decimals=2)
#                 selX=xppm[Xind]
#                 selY=yppm[Yind]
#                 selP=np.vstack((selX,selY))
#                 P[j]=selP
#                 X[j]=selX
#                 Y[j]=selY
#                 Pts[j]=list(zip(selX,selY))
#                 j=j+1
#         return (Pts,X,Y,P) 

def GetPoints (List,ind):
	Points=[]
	for i in ind:
		Points.append(List[i])
	return(Points)

# @profile
def Sort (Points, thres, ax):
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
		if (abs(minval[ax]-elem[ax]) < thres):
			sets[i].append(elem)
		else:
			i=i+1
			sets[i].append(elem)
	# print sets
	return(sets)

def YSort (Points, thres, ax):
	sets = defaultdict(list)
	i=0
	for k, Group in list(Points.items()):
		sPts = sorted(Group, key=lambda tup: tup[ax])	
		sets[i].append(sPts[0])
		for first, second in zip(sPts, sPts[1:]):
			if (abs(first[1]-second[1]) < thres):
				sets[i].append(second)
			else:
				i=i+1
				sets[i].append(second)
		i=i+1
	return(sets)

def HorzAlign2 (P,threshold1,threshold2,threshold3,ct):
	Pd=dict()
	Pd[0]=P
	AlnY=YSort(Pd,threshold1,1)
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
				if (abs(meanY-sumX)<=threshold2):
					D1=abs(a[1][0]-a[1][1]/2) # {diag:y=2x;mid-point:(y/2,y) so dist=x-y/2}
					D2=abs(b[1][0]-b[1][1]/2)
					if (abs(D1-D2)<=threshold3):
						filtI[ct].append(a[1])
						filtI[ct].append(b[1])
						# print filtI[ct]
						ct+=1
						# print "===>Passed Diag", D1,D2,threshold3,"Y=",a[1][1],b[1][1],"X=",a[1][0],b[1][0]
					# else:
						# print "#####Did not pass Diag", D1,D2,threshold3,"Y=",a[1][1],b[1][1],"X=",a[1][0],b[1][0]
	return (filtI,ct)
	# return (filtI)

def CenterMass2 (Clusters,metric):
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

def FindClosestPoints (Pts,Qry):
	PtsList = sp.KDTree(Pts)
	QryRes=PtsList.query(Qry)
	return(QryRes)

def HorzAlign (X,Y,threshold1,threshold2,threshold3,ct): 	# Find peaks horizontally aligned
# GetPoints
	(AlnY,AlnYI)=Gather(Y,threshold1)
	j=0
	XLI=[]
	filtI=defaultdict(list)
	filtX=[]
	filtY=[]
	for i,lst in AlnYI.items():
		# print "lst==>>",lst
		if (len(lst)>=2):
			for a, b in itertools.combinations(enumerate(lst), 2):
				# print "a,b==>",a,b
				meanY=(Y[a[1]] + Y[b[1]]) /2
				sumX=X[a[1]] + X[b[1]]
				if (abs(meanY-sumX)<=threshold2):
					D1=abs((X[a[1]]-(Y[a[1]])/2))
					D2=abs((X[b[1]]-(Y[b[1]])/2))
					if (abs(D1-D2)<=threshold3):
						filtI[ct].append((X[a[1]],Y[a[1]]))
						filtI[ct].append((X[b[1]],Y[b[1]]))
						ct+=1
	return (filtI,ct)

## OLD Unused functions ###

def Gather (Points, thres):	# Cluster groups of peaks along x-axis
	sets = defaultdict(list)
	all_list = defaultdict(list)
	index= defaultdict(list)
	list1=[]
	ind_list1=[]
	ind_list2=[]
	ind_list3=[]
	i=0
	init=0

	if (len(Points)==1):
		sets[i].append(Points[0])
		index[i].append(i)
	for a, b in itertools.combinations(enumerate(Points), 2):
		all_list[a[0]]=a[1];
		all_list[b[0]]=b[1];
		if init == a[0]:
			num=a[1]
			ind=a[0]
			if (abs(a[1]-b[1]) < thres):
				if a[0] not in ind_list1:
					sets[i].append(a[1])
					index[i].append(a[0])
					ind_list1.append(a[0])
				if b[0] not in ind_list2:
					sets[i].append(b[1])
					index[i].append(b[0])
					ind_list2.append(b[0])
			else:
				ind_list3.append(a[0])
				ind_list3.append(b[0])
		else:
			if i not in sets:
				sets[i].append(num)
				index[i].append(ind)
				ind_list1.append(ind)
			if a[0] in ind_list3 and a[0] not in ind_list2:
				init = a[0]
				i+=1
				if (abs(a[1]-b[1]) < thres):
					if a[0] not in ind_list1:
						sets[i].append(a[1])
						index[i].append(a[0])
						ind_list1.append(a[0])
					if b[0] not in ind_list2:
						sets[i].append(b[1])
						index[i].append(b[0])
						ind_list2.append(b[0])
				else:
					ind_list3.append(b[0])
	for ind, value in list(all_list.items()):
		if (ind not in ind_list1 and ind not in ind_list2):
			if i in sets:
				i +=1
			sets[i].append(value)
			index[i].append(ind)
			ind_list1.append(ind)
	return (sets, index)

def Split (Points,thres):	# Cluster groups of peaks along y-axis
	Psort = np.argsort(Points)
	Points.sort(key=int)
	sets2 = defaultdict(list)
	index2 = defaultdict(list)
	i=0
	for j in range(len(Points) - 1):
		curr, nxt = Points[j], Points[j + 1]
		if (j==0):
			sets2[i].append(curr)
			index2[i].append(Psort[j])
		lg = max(sets2[i])
		sm = min(sets2[i])
		if abs(lg-nxt) < thres or abs(sm-nxt) < thres:
			sets2[i].append(nxt)
			index2[i].append(Psort[j+1])
		else:
			i+=1
			sets2[i].append(nxt)
			index2[i].append(Psort[j+1])
	return (sets2, index2)


def ClustY (Xindex,Y,thres):
#Split
	YgroupIndex=defaultdict(list)
	for i in Xindex:
		Ylist=[]
		Yorig_ind=defaultdict(list)
		k=0
		for j in Xindex[i]:
			Ylist.append(Y[j])
			Yorig_ind[k]=j
			k+=1
		(Ysets,Yindex)=Split(Ylist,thres)
		for m in Yindex:
			for n,val in enumerate(Yindex[m]):
				Yindex[m][n] = Yorig_ind[Yindex[m][n]]
		YgroupIndex[i].append(Yindex)
	return(YgroupIndex)

def MergeLevels (indL):
	Lev=[]
	for i,val in indL.items():
		for j in val:
			for k,val2 in j.items():
				Lev.append(val2)	
	return(Lev)

def CenterMass (X,Y,indexList):	# Calculate center of mass for each Cluster
# GetPoints
	i=0
	Center=[]
	C=[]
	for val in indexList:
			XP=GetPoints(X,val)
			meanX=sum(XP) / float(len(XP))
			YP=GetPoints(Y,val)
			meanY=sum(YP) / float(len(YP))
			C.append(meanX)
			C.append(meanY)
			Center.append(C)
			C=[]
			i+=1
	return(np.asarray(Center))
	# return(Center)

## HorzAlign here

def BuildNetwork (HorzPts,VertPts):	# Make network from selected peaks
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
		g = networkx.Graph()
	for lst in AllPts:
		lstt=tuple(map(tuple, lst))
		for edge in getPairs(lstt):
			# print type(edge)
			# print type(lstt)
			g.add_edge(*edge)
	Net=list(networkx.connected_components(g))
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
	return (Net,Netf,j)

def ParseDbNetwork (Name,Net):
	IndDict={}
	IndNetLab=[]
	IndNetX=[]
	IndNetY=[]
	IndNetPts=[]
	# print Net
	for m in Net:	# json_db[i]['Networks']: # Each pair of points
		IndNetPts.extend([m[0][1], m[1][1]])
		for j in m:	# Each point in pair
			allX=[]
			allY=[]
			for k in j[1]:	# Only points without label
				allX.append(k[0])
				allY.append(k[1])
			Xval=sum(allX)/float(len(allX))
			Yval=sum(allY)/float(len(allY))
			IndNetLab.append(j[0])
			IndNetX.append(Xval)
			IndNetY.append(Yval)
	IndDict.setdefault(Name, []).append(IndNetLab)
	IndDict[Name].append(IndNetX)
	IndDict[Name].append(IndNetY)
	IndDict[Name].append(IndNetPts)
	IndDict[Name].append(Net)
	return(IndDict)


def MatchTag (match_peaks,candidate_peaks,unknown_conn):
	final_match={}
	no_match={}
	hits=0
	for i in unknown_conn:
		key=i[0]+"-"+i[1]
		if i[0] in match_peaks:
			if i[1] in match_peaks:
				final_match.setdefault(key, []).append(match_peaks[i[0]])
				final_match[key].append(match_peaks[i[1]])
			else:
				no_match.setdefault(key, []).append(match_peaks[i[0]])
				no_match[key].append("?")
		else:
			no_match.setdefault(key, []).append("?")
			no_match[key].append("?")
	for key,val in list(final_match.items()):
		if (val[0] != val[1]):
			hits=hits+1

	# hits=len(final_match.keys())
	# print "Match",len(final_match.keys()),hits,final_match
	return (final_match,no_match,hits)

def get_cmap(n, name='Set1'):
    return plt.cm.get_cmap(name, n)

def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
        
    return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

def MatchDatabase (json_db,Pval,tag,unknown_conn,amb_tol,near_tol,match_tol,top_tol,hit_tol,cov_tol,NetNum,Pts):
#MatchTag
	X = [float(i[0]) for i in Pval]
	Y = [float(i[1]) for i in Pval]
	hitList=[]
	hitList.append("Network"+ str(NetNum))
	hitList.append(X)
	hitList.append(Y)
	hitList.append(tag)
	hitList.append(Pts)
	# print X
	# print Y
	# print Pval
	for i in json_db:
		# print i
		matchCnt=0
		matchCS=[]
		if (len(json_db[i]['Networks'] ) >=1):
			# print json_db[i]['Networks'], len(json_db[i]['Networks'] )
			if (json_db[i]['Ambiguity'] <= amb_tol):
				# print json_db[i]['Ambiguity'],"===", amb_tol
				for xval in X:
					for CS in json_db[i]['ChemicalShifts']:
						if CS not in matchCS:
							CSlist=json_db[i]['ChemicalShifts'][CS]
							meanCS=sum(CSlist) / float(len(CSlist))
							# print i, xval, CS,"////", meanCS, matchCS
							if (abs(xval-meanCS)<=near_tol):
								# print i, xval, meanCS, "==", abs(xval-meanCS), "//", near_tol
								matchCS.append(CS)
								matchCnt+=1
				if (matchCnt>=match_tol):
					# print i, matchCS
					found={}
					candPeaks={}
					matchPeaks={}
					for j in json_db[i]['Networks']:
						# print "j==>",j
						candPeaks.setdefault(j[0][0], set()).add(j[1][0])
						for k in j:
							# print "k==>",k
							break_outer_loop=False
							for l in k[1]:
								# print "l==>",l
								for ind,xval in enumerate(X):
									mp=[]
									# print ind,xval,Y[ind],tag[ind]
									# print i,"=",k[0],l,"//",tag[ind],xval,Y[ind],"//",(xval-l[0])**2,(Y[ind]-l[1])**2,((xval-l[0])**2 + (Y[ind]-l[1])**2)**(1/2.0),top_tol
									if (((xval-l[0])**2 + (Y[ind]-l[1])**2)**(1/2.0) <= top_tol):										
										matchPeaks[tag[ind]]=k[0]
										# print "matchPeaks",len(tag),tag[ind],matchPeaks
										# print i,"HEE",l,xval,Y[ind],(xval-l[0])**2,(Y[ind]-l[1])**2,((xval-l[0])**2 + (Y[ind]-l[1])**2)**1/2,top_tol
										found[tag[ind]]=1
										break_outer_loop=True
										break
								if break_outer_loop: break
					# print i
					# print "matchPeaks==>",matchPeaks
					# print "candPeaks==>",candPeaks
					# print "unknown_conn==>",unknown_conn
					(fin_match,no_match,hitCount)=MatchTag(matchPeaks,candPeaks,unknown_conn)
					hitScore=float("{0:.3f}".format(float(hitCount)/float(len(json_db[i]['Networks']))))
					CovScore=float("{0:.3f}".format(float(len(found))/float(len(tag))))
					# print len(found),found,len(tag),tag 
					# print "fin_match==>",fin_match,
					# print "no_match==>",no_match,
					# print "hitCount=",hitCount," Total=",len(json_db[i]['Networks']),"=> hitScore=",hitScore
					if (hitScore >= hit_tol and CovScore >= cov_tol):
						hitOut=[]
						hitOut.append(i)
#						print json_db[i]['Networks']
						hitOut.append([json_db[i]['Networks']])
						# IndDict=ParseDbNetwork(i, json_db[i]['Networks'])
						# hitOut.append(IndDict)
						hitOut.append(fin_match)
						hitOut.append(no_match)
						hitOut.append(hitScore)
						hitOut.append(CovScore)
						hitList.append(hitOut)
						# print hitList
	return(hitList)

def MatchMetab (json_db,metab_list,Pval,tag,unknown_conn,amb_tol,near_tol,match_tol,top_tol,hit_tol,cov_tol,NetNum,Pts):
	X = [float(i[0]) for i in Pval]
	Y = [float(i[1]) for i in Pval]
	hitList=[]
	hitList.append("Network"+ str(NetNum))
	hitList.append(X)
	hitList.append(Y)
	hitList.append(tag)
	hitList.append(Pts)
	print("===")
	print("Network"+str(NetNum),":")
	print("Pvals:",np.asarray(Pval))
	print("Net:",np.asarray(X))
	# print X
	# print Y
	for i in json_db:
		name_parts=i.split('::')
		if name_parts[2] in metab_list:
			print("***Matching to",name_parts[2])	
			matchCnt=0
			matchCS=[]
			print("Length of network for",name_parts[2],"is",len(json_db[i]['Networks']))
			print("Number of peaks for",name_parts[2],"is",len(json_db[i]['ChemicalShifts']))
			print("Ambiguity is "+str(json_db[i]['Ambiguity']))
			print("DB:",json_db[i]['ChemicalShifts'])
			print("Net:",np.asarray(X))
			print("Peak by peak Matches: (NetPeak DBPeak Diff // CSMT",near_tol,")")
			if (len(json_db[i]['Networks'] ) >=1):
				# print json_db[i]['Networks'], len(json_db[i]['Networks'] )
				if (json_db[i]['Ambiguity'] <= amb_tol):
					# print json_db[i]['Ambiguity'],"===", amb_tol
					for xval in X:
						for CS in json_db[i]['ChemicalShifts']:
							if CS not in matchCS:
								CSlist=json_db[i]['ChemicalShifts'][CS]
								meanCS=sum(CSlist) / float(len(CSlist))
								# print " ",xval, CS,"////", meanCS, matchCS
								if (abs(xval-meanCS)<=near_tol):
									print(" ",xval,CS,meanCS, "==", abs(xval-meanCS), "//", near_tol, matchCS)
									matchCS.append(CS)
									matchCnt+=1
			
			if (matchCnt>=match_tol):
				print(i, matchCS)
				found={}
				candPeaks={}
				matchPeaks={}
				for j in json_db[i]['Networks']:
					# print "j==>",j
					candPeaks.setdefault(j[0][0], set()).add(j[1][0])
					for k in j:
						# print "k==>",k
						break_outer_loop=False
						for l in k[1]:
							# print "l==>",l
							for ind,xval in enumerate(X):
								mp=[]
								# print ind,xval,Y[ind],tag[ind]
								# print i,"=",k[0],l,"//",tag[ind],xval,Y[ind],"//",(xval-l[0])**2,(Y[ind]-l[1])**2,((xval-l[0])**2 + (Y[ind]-l[1])**2)**(1/2.0),top_tol
								if (((xval-l[0])**2 + (Y[ind]-l[1])**2)**(1/2.0) <= top_tol):										
									matchPeaks[tag[ind]]=k[0]
									# print "matchPeaks",len(tag),tag[ind],matchPeaks
									# print i,"HEE",l,xval,Y[ind],(xval-l[0])**2,(Y[ind]-l[1])**2,((xval-l[0])**2 + (Y[ind]-l[1])**2)**1/2,top_tol
									found[tag[ind]]=1
									break_outer_loop=True
									break
							if break_outer_loop: break
				# print i
				# print "matchPeaks==>",matchPeaks
				# print "candPeaks==>",candPeaks
				# print "unknown_conn==>",unknown_conn
				(fin_match,no_match,hitCount)=MatchTag(matchPeaks,candPeaks,unknown_conn)
				hitScore=float("{0:.3f}".format(float(hitCount)/float(len(json_db[i]['Networks']))))
				CovScore=float("{0:.3f}".format(float(len(found))/float(len(tag))))
				print(len(found),found,len(tag),tag) 
				# print "fin_match==>",fin_match,
				# print "no_match==>",no_match,
				print("hitCount=",hitCount," Total=",len(json_db[i]['Networks']),"=> hitScore=",hitScore,"=> CovScore=",CovScore)
				if (hitScore >= hit_tol and CovScore >= cov_tol):
					hitOut=[]
					hitOut.append(i)
#							print json_db[i]['Networks']
					hitOut.append([json_db[i]['Networks']])
					# IndDict=ParseDbNetwork(i, json_db[i]['Networks'])
					# hitOut.append(IndDict)
					hitOut.append(fin_match)
					hitOut.append(no_match)
					hitOut.append(hitScore)
					hitOut.append(CovScore)
					hitList.append(hitOut)
					print(hitList)
	return(hitList)


def MatchLines (query, json_db, tol1, tol2):
	hitList=list()
	for i in json_db:
		match_ct=0;
		match_query=list()
		match_db=list()
		hitOut=list()
		# print "DB==>",i,map(lambda x: isinstance(x, float) and round(x, 3) or x, json_db[i]['Lines'])
		for j in json_db[i]['Lines']:
			for k in query:
				if (abs(j-k)<=float(tol1) and k not in match_query and j not in match_db):
					# print "HITHIT",j, k, abs(j-k),type(abs(j-k)), type(tol1), tol1, query
					match_ct+=1
					match_query.append(k)
					match_db.append(j)
		# print len(json_db[i]['Lines']),"//", match_ct	
		# if (abs(len(json_db[i]['Lines'])-match_ct)<=float(tol2) and abs(len(query)-match_ct)<=float(tol2) and match_ct > 0):
		if (len(json_db[i]['Lines'])==match_ct and len(query)==match_ct and match_ct > 0):
			# print i, len(json_db[i]['Lines']), len(query),match_ct,type(abs(len(json_db[i]['Lines'])-match_ct)), type(tol2), type(len(query))
			hitOut.append(i)
			hitOut.append(len(json_db[i]['Lines']))
			hitOut.append(len(query))
			hitOut.append([isinstance(x, float) and round(x, 3) or x for x in match_query])
			hitOut.append([isinstance(x, float) and round(x, 3) or x for x in match_db])
			hitOut.append(match_ct)
			hitList.append(hitOut)
	return(hitList)

def plotDb (dbNet,ax,col,pos,name):
	# print "dbNet==>",dbNet
	X=[]
	Y=[]
	ID=[]
	name_parts=name.split('::')
	for Net in dbNet:
		Id1=Net[0][0]#.pop(0)
		Id2=Net[1][0]#.pop(0)
		# print "Ids==>",Id1,Id2
		# print "Net==>",len(Net[0][1]),Net[0][1]
		for j in range(0,len(Net[0][1])):
			x1=Net[0][1][j][0]
			x2=Net[1][1][j][0]
			y1=Net[0][1][j][1]
			y2=Net[1][1][j][1]
			if x1 in X:
				x3=x1
				y3=Y[X.index(x1)]
				# print "x3,y3=",x1,y1,x3,y3
				ax.plot([x1,x3],[y1,y3],color=col, linewidth=2)
				#get x1,y1
			X.append(x1)
			Y.append(y1)
			ID.append(Id1)
			if x2 in X:
				x3=x2
				y3=Y[X.index(x2)]
				# print "x3,y3=",x2,y2,x3,y3
				ax.plot([x2,x3],[y2,y3],color=col, linewidth=2)
				#get x2,y2
			X.append(x2)
			Y.append(y2)
			ID.append(Id2)
			# print Net[0][0][j][0],Net[1][0][j][0]
			# print Net[0][0][j][1],Net[1][0][j][1]
			# print x1,x2,y1,y2
			ax.plot([x1,x2],[y1,y2],color=col, linewidth=2)
		# print Net[0]
		# print Net[1]
		# print "---"
		# print "X=",X
		# print "Y=",Y
#	print X,Y,ID
	ax.scatter(X,Y,c=col, s=100)
	for label, x, y in zip(ID,X,Y):
		ax.annotate(label, xy=(x, y), xytext=(40, -40),fontsize='25',textcoords='offset points', ha='right', va='bottom',color=col)
	ax.text(2,pos*5+5,name_parts[2],fontsize='20',color=col,ha='right')
	return(ax)

#### Functions not used #####

# def peakPickNew (points,xppm,yppm,PPmin,PPmax,steps):
# 	j=0
# 	X={}
# 	Y={}
# 	P={}
# 	Pts={}
# 	for i in pl.frange(int(PPmax),int(PPmin),npts=steps):
# 		# print i
# 		(Xind, Yind) = np.where(abs(points) > i)
# 		points[abs(points) > i] =0          ## This line removes hit peaks from the previous iteration (do we do this??)
# 		selX=np.around(xppm[Xind],decimals=2)
# 		selY=np.around(yppm[Yind],decimals=2)
# 		selX=xppm[Xind]
# 		selY=yppm[Yind]
# 		selP=np.vstack((selX,selY))
# 		P[j]=selP
# 		X[j]=selX
# 		Y[j]=selY
# 		Pts[j]=list(zip(selX,selY))
# 		j=j+1
# 	return (Pts,X,Y,P)

# def printRound(y):
#     for j, dk in y.items():
#     	print(j,": ", end=' ')
#     	for d in dk:
#     		print("%.2f" % d, end=' ')
#     	print()







	