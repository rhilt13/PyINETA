import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from shutil import copyfile

from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')

def plotSingle (pyinetaObj,Xlim,Ylim,label,imgname=None,grid=True):
	try:
		X=pyinetaObj.Xlist[len(pyinetaObj.Xlist)-1]
		Y=pyinetaObj.Ylist[len(pyinetaObj.Ylist)-1]
	except AttributeError:
		X=pyinetaObj[0]
		Y=pyinetaObj[1]
	if (type(label)==float):
		title="{:.2e}".format(label[i])
	else:
		title=label
	fig=plt.figure(figsize=(30, 30)) # This increases resolution
	ax = fig.add_subplot(111)
	ax.plot(X,Y,'k.',)
	ax.set_xlim(Xlim)
	ax.set_ylim(Ylim)
	ax.set_ylabel('Double Quantum', fontsize=30)
	ax.set_xlabel('13C', fontsize=30)
	ax.set_title(title, fontsize=30)
	ax.tick_params(axis='both', which='major', labelsize=25)
	ax.invert_xaxis()
	ax.invert_yaxis()
	plt.tight_layout(pad=5, w_pad=0.5, h_pad=1.0)
	plt.minorticks_on()
	if grid:
		plt.grid(b=True, which='major', color='0.45', linestyle=':', linewidth=0.2)
		plt.grid(b=True, which='minor', color='0.75', linestyle=':', linewidth=0.1)
	if imgname is not None:
		fig.savefig(imgname,format='eps',dpi=300)
	plt.close(fig)
	return(fig, ax)

def plotFigSep (pyinetaObj,Xlim,Ylim,label,imgname=None):
	X=pyinetaObj.Xlist
	Y=pyinetaObj.Ylist
	r=np.ceil(float(len(X))/2)
	c=2
	l=c*30
	w=r*30
	fig, axs = plt.subplots(nrows=int(r), ncols=c, figsize=(l,w))
	for i, ax in enumerate(fig.axes):
		if i in X:
			title="{:.2e}".format(label[i])
			ax.plot(X[i],Y[i],'k.')
			ax.set_xlim(Xlim)
			ax.set_ylim(Ylim)
			ax.set_ylabel('Double Quantum', fontsize=60)
			ax.set_xlabel('13C', fontsize=60)
			ax.set_title(title, fontsize=60)
			ax.tick_params(axis='both', which='major', labelsize=40)
			ax.invert_xaxis()
			ax.invert_yaxis()
			plt.minorticks_on()
			plt.grid(b=True, which='major', color='0.45', linestyle=':', linewidth=0.25)
			plt.grid(b=True, which='minor', color='0.75', linestyle=':', linewidth=0.2)
	plt.tight_layout(pad=10, w_pad=10, h_pad=10)
	if (len(X)%2!=0):
		axs[-1,-1].axis('off')
	if imgname is not None:
		fig.savefig(imgname,format='eps',dpi=300)
	plt.close(fig)
	return(fig,axs)

def plotFilteredPoints (figSep,axsSep,figCom,axsCom,pyinetaObj,PPcs,PPdq,out_sep2,out_comp2):
	panels=len(pyinetaObj.Xlist)
	Pts=pyinetaObj.filteredPts
	AllC=[]
	for k, C in Pts.items():
		r=int(k/2)
		if (k%2==0):
			c=0
		else:
			c=1
		for i in C:
			circ=plt.Circle(i,radius=min(PPcs,PPdq),color='b',fill=False,linestyle='dotted',linewidth=0.8)
			circ2=plt.Circle(i,radius=min(PPcs,PPdq),color='b',fill=False,linestyle='dotted',linewidth=0.8)
			if (panels>1):
				axsSep[r,c].add_patch(circ)
			axsCom.add_patch(circ2)
		AllC.extend(C)
	if (panels>1):
		figSep.savefig(out_sep2,dpi=300) 	# Separate
	figCom.savefig(out_comp2,dpi=300) 		# Complete

def plotAllNet (horzPts,vertPts,ax,rad):
	for i,v in list(horzPts.items()):
		for q in v:
			circ=plt.Circle(q,radius=rad,color='b',fill=False,linestyle='dotted',linewidth=0.8)
			ax.add_patch(circ)
		ax.plot([v[0][0], v[1][0]], [v[0][1], v[1][1]], color='#f72d1c', linestyle='--', linewidth=2)
	for i,v in list(vertPts.items()):
		for a, b in combinations(enumerate(v),2):
			ax.plot([a[1][0], b[1][0]], [a[1][1], b[1][1]], color='#f72d1c', linestyle='--', linewidth=2)
	return (ax)

def plotNetwork (pyinetaObj,Xrng,Yrng,PPcs,PPdq,out_nets=None):
	(figFin,axsFin)=plotSingle(pyinetaObj,Xrng,Yrng,'Complete plot',grid=False)
	## Getting the Diagonal
	midXL=[]
	midYL=[]
	for i,v in pyinetaObj.horzPts.items():
		midX=(v[0][0]+v[1][0])/2
		midY=(v[0][1]+v[1][1])/2
		midXL.append(midX)
		midYL.append(midY)
	indMax = np.argmax(midXL)
	indMin = np.argmin(midXL)
	x1=midXL[indMin]
	y1=midYL[indMin]
	x2=midXL[indMax]
	y2=midYL[indMax]
	xlim=axsFin.get_xlim()
	ylim=axsFin.get_ylim()
	line_eqn = lambda x : ((y2-y1)/(x2-x1)) * (x - x1) + y1        
	## generate range of x values based on your graph
	xvals=np.arange(xlim[1],xlim[0],10)
	## plot the line with generate x ranges and created y ranges
	axsFin.plot(xvals, [ line_eqn(x) for x in xvals], color='g', linestyle='-', linewidth=2)
	axsFin=plotAllNet(pyinetaObj.horzPts,pyinetaObj.vertPts,axsFin,max(PPcs,PPdq))
	if out_nets is not None:
		figFin.savefig(out_nets,format='eps',dpi=300)

def colors(n):
	# Custom colormap for visually distinct colors:
	# Colors are: orangered, darkblue,forestgreen,maroon3,fuchsia,cornflower,lime,aqua,moccasin,yellow
	col_list=['#ff4500','#00008b','#228b22','#b03060','#ff00ff','#6495ed','#00ff00','#00ffff','#ffe4b5','#ffff00']
	return col_list[n]

def plotDb (dbNet,ax,col,pos,name):
	X=[]
	Y=[]
	ID=[]
	name_parts=name.split('::')
	for Net in dbNet:
		Id1=Net[0][0]
		Id2=Net[1][0]
		for j in range(0,len(Net[0][1])):
			x1=Net[0][1][j][0]
			x2=Net[1][1][j][0]
			y1=Net[0][1][j][1]
			y2=Net[1][1][j][1]
			if x1 in X:
				x3=x1
				y3=Y[X.index(x1)]
				ax.plot([x1,x3],[y1,y3],color=col, linewidth=2)
				#get x1,y1
			X.append(x1)
			Y.append(y1)
			ID.append(Id1)
			if x2 in X:
				x3=x2
				y3=Y[X.index(x2)]
				ax.plot([x2,x3],[y2,y3],color=col, linewidth=2)
				#get x2,y2
			X.append(x2)
			Y.append(y2)
			ID.append(Id2)
			ax.plot([x1,x2],[y1,y2],color=col, linewidth=2)
	ax.scatter(X,Y,c=col, s=100)
	for label, x, y in zip(ID,X,Y):
		ax.annotate(label, xy=(x, y), xytext=(40, -40),fontsize='25',textcoords='offset points', ha='right', va='bottom',color=col)
	ax.text(2,pos*5+5,name_parts[2],fontsize='20',color=col,ha='right')
	return(ax)

def plotMatches (pyinetaObj,outFolder,db_file,Xrng,Yrng):
	
	# Initialize folders and output files:
	outNetFigs=outFolder+"/Network_figs/"
	if not os.path.exists(outNetFigs):
		os.makedirs(outNetFigs)
	outMatchFigs=outFolder+"/Matches/"
	if not os.path.exists(outMatchFigs):
		os.makedirs(outMatchFigs)
	li = db_file.rsplit('/')
	li.pop()
	src="/".join(li)
	# Plotting individual network and matches
	for i in pyinetaObj.NetMatch:
		col=0
		pos=0
		netfig=outNetFigs+i[0]+".eps"
		pts=list()
		pts.append(i[1])
		pts.append(i[2])
		(figCom,axsCom)=plotSingle(pts,Xrng,Yrng,i[0],netfig,grid=False)
		axsCom.scatter(i[1],i[2], c='0.45', s=100,marker='o',facecolors='none')
		for label, x, y in zip(i[3],i[1],i[2]):	
			axsCom.annotate(label, xy=(x, y), xytext=(-20, 20),fontsize='25',textcoords='offset points', ha='right', va='bottom')
		for j in i[4]:
			axsCom.plot([j[0][0], j[1][0]], [j[0][1], j[1][1]], color='0.45', linestyle='--', linewidth=4)
		if len(i)>5:
			for j in range(5,len(i)):
				axsCom=plotDb(i[j][1][0],axsCom,colors(col),pos,i[j][0])
				for label, x, y in zip(i[3],i[1],i[2]):
					axsCom.annotate(label, xy=(x, y), xytext=(-20, 20),fontsize='25',textcoords='offset points', ha='right', va='bottom')
				name=i[j][0]
				name=name.replace('::','-')
				name=name.replace('/','-')
				name +='.eps'
				srcfile=src+"/db_images_constXYlim/"+name
				destfile=outMatchFigs+name
				copyfile(srcfile,destfile)
				pos+=1
				col+=1
				if col==10:
					col=0
		figCom.savefig(netfig, dpi=300)
