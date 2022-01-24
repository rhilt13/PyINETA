"""Functions for plotting INETA results

This includes all the plotting functions used by INETA.
Includes the following funcitons:
	* plot1D - Plot 1D spectra
	* plotSingle - Generate a single plot for a given spectrum.
	* plotFigSep - Generate image with spectra for each saved iteration in a separate panel.
	* plotClusteredPoints - Plot the clustered points and theirr centers in the spectra.
	* plotAllNet - Generate a plot with all the networks.
	* plotNetwork - Plot all the networks overlaid in the spectra figure.
	* colors - Get a specific color from a list of custom colormap.
	* plotDb - Generate plots for all networks in the INETA database.
	* plotMatches - Plot the matches for a network.
	* plotIndividualMatch - Plot individual network with individual database entries of choice.

"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from shutil import copyfile
from itertools import combinations
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')

def plot1D (Intsy,ppm,range,ax=None,title=None,net=None):
	sc=30
	ax = ax or plt.gca()
	ax.plot(ppm,Intsy,'k-')
	ax.tick_params(axis='x', which='major', length=5, labelsize=(sc*0.8))
	ax.xaxis.set_minor_locator(MultipleLocator(5))
	ax.tick_params(which='minor', length=2)
	if title is not None:
		ax.set_title(title)
	ax.set_xlabel("13C ppm")
	i=0
	for key,vals in range.items():
		j=0
		patch=list()
		for r in vals:
			if net is None:
				if (j==0):
					ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25, label=key)
				else:
					ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25)
				j+=1
			else:
				# If specific network is provided
				Netlist=net.lower().split(",")
				if key.lower() in Netlist:
					if (j==0):
						ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25, label=key)
					else:
						ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25)
					j+=1
		i+=1
		if i==10:
			i=0
	return(ax)

def plotNetWith1D (pyinetaObj):
	try:
		X=pyinetaObj.Xlist[len(pyinetaObj.Xlist)-1]
		Y=pyinetaObj.Ylist[len(pyinetaObj.Ylist)-1]
	except AttributeError:
		X=pyinetaObj[0]
		Y=pyinetaObj[1]
	x, y = np.random.randn(2, 100)
	fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True)
	ax1.xcorr(x, y, usevlines=True, maxlags=50, normed=True, lw=2)
	ax1.grid(True)

	ax2.acorr(x, usevlines=True, normed=True, maxlags=50, lw=2)
	ax2.grid(True)

	plt.show()

def plotJres (In,Xs,Ys,ax=None):
	sc=30
	ax = ax or plt.gca()

	ax.plot(X,Y,'k.',markersize=2)
	# ax.xaxis.grid(True, which='minor')
	ax.tick_params(axis='x', which='major', length=5, labelsize=(sc*0.8))
	ax.xaxis.set_minor_locator(MultipleLocator(5))
	ax.tick_params(which='minor', length=2)
	# ax.set_xlim(Xlim)
	# ax.set_yticklabels([])
	if title is not None:
		ax.set_title(title)
	ax.set_xlabel("13C ppm")
	# ax.set_xlim(200, 0)
	# ax.set_ylim(-80000, 2500000)
	i=0
	for key,vals in range.items():
		j=0
		patch=list()
		for r in vals:
			if net is None:
				if (j==0):
					ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25, label=key)
				else:
					ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25)
				j+=1
			else:
				# If specific network is provided
				Netlist=net.lower().split(",")
				if key.lower() in Netlist:
					if (j==0):
						ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25, label=key)
					else:
						ax.axvspan(r[0], r[1], facecolor=colors(i), alpha=0.25)
					j+=1
		i+=1
		if i==10:
			i=0
	return(ax)


def plotSingle (pyinetaObj,Xlim,Ylim,label,imgname=None,grid=False):
	"""Generate a single plot for a given spectrum.

	Args:
		pyinetaObj (pyineta object): The pyineta object with the spectra.
		Xlim (tuple): Tuple of (min,max) values of the X axis (13C spectrum); (0,200) for a typical INADEQUATE spectra.
		Ylim (tuple): Tuple of (min,max) values of the Y axis (DQ spectrum); (0,400) for a typical INADEQUATE spectra.
		label (str): Plot title
		imgname (str, optional): Filename to save the image. Not saved if not specified. Defaults to None.
		grid (bool, optional): Show grid in the plot. Defaults to False.

	Returns:
		figure and axis handles for the generated image.
	"""

	sc=30
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
	fig=plt.figure(figsize=(sc, sc)) # This changes figure size and resolution
	ax = fig.add_subplot(111)
	ax.plot(X,Y,'k.',markersize=2)
	ax.set_xlim(Xlim)
	ax.set_ylim(Ylim)
	ax.set_ylabel('Double Quantum', fontsize=sc)
	ax.set_xlabel('13C', fontsize=sc)
	ax.set_title(title, fontsize=sc)
	ax.tick_params(axis='both', which='major', labelsize=(sc*0.8))
	ax.invert_xaxis()
	ax.invert_yaxis()
	plt.tight_layout(pad=(sc*0.15), w_pad=0.5, h_pad=1.0)
	plt.minorticks_on()
	if grid:
		plt.grid(b=True, which='major', color='0.45', linestyle=':', linewidth=0.2)
		plt.grid(b=True, which='minor', color='0.75', linestyle=':', linewidth=0.1)
	if imgname is not None:
		fig.savefig(imgname,format='eps',dpi=300)
	plt.close(fig)
	return(fig, ax)

def plotFigSep (pyinetaObj,Xlim,Ylim,label,imgname=None):
	"""Generate image with spectra for each saved iteration in a separate panel.

	Args:
		pyinetaObj (pyineta object): The pyineta object with the spectra.
		Xlim (tuple): Tuple of (min,max) values of the X axis (13C spectrum); (0,200) for a typical INADEQUATE spectra.
		Ylim (tuple): Tuple of (min,max) values of the Y axis (DQ spectrum); (0,400) for a typical INADEQUATE spectra.
		label (str): Plot title.
		imgname (str, optional): Filename to save the image. Not saved if not specified. Defaults to None.

	Returns:
		figure and axis handles for the generated image.
	"""

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

def plotClusteredPoints (figSep,axsSep,figCom,axsCom,pyinetaObj,PPcs,PPdq,out_sep2,out_comp2):
	"""Plot the clustered points and theirr centers in the spectra.
	
	This function generates an image of the spectra with dotted circles around every center of mass point to indicate the range of clustering. 

	Args:
		figSep (figure): Figure handle for the plot with separate panels.
		axsSep (axes): Axes handle for the plot with separate panels.
		figCom (figure): Figure handle for the plot with a single panel.
		axsCom (axes): Axes handle for the plot with a single panel.
		pyinetaObj (pyineta object): The pyineta object with the spectra.
		PPcs (float): ppm threshold to cluster points along the 13C axis.
		PPdq (float): ppm threshold to cluster points along the DQ axis.
		out_sep2 (str): Filename to save the plot with separate panels.
		out_comp2 (str): Filename to save the plot with single panel.
	
	"""

	panels=len(pyinetaObj.Xlist)
	Pts=pyinetaObj.clusteredPts
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
		figSep.savefig(out_sep2,dpi=300) 	# Separate plots
	figCom.savefig(out_comp2,dpi=300) 		# Complete single plot with lowest level

def plotAllNet (horzPts,vertPts,ax,rad):
	"""Generate a plot with all the networks.

	Args:
		horzPts (dict): A dict mapping horizontally aligned peaks to their indices.
		vertPts (dict): A dict mapping vertically aligned peaks to their indices.
		ax ([type]): Axes handle for the plot with a single panel.
		rad ([type]): clustering radius shown using dotted circles.

	Returns:
		axes handle for the resulting figure
	"""

	for i,v in list(horzPts.items()):
		for q in v:
			circ=plt.Circle(q,radius=rad,color='b',fill=False,linestyle='dotted',linewidth=0.8)
			ax.add_patch(circ)
		ax.plot([v[0][0], v[1][0]], [v[0][1], v[1][1]], color='#f72d1c', linestyle='--', linewidth=2)
	for i,v in list(vertPts.items()):
		for a, b in combinations(enumerate(v),2):
			ax.plot([a[1][0], b[1][0]], [a[1][1], b[1][1]], color='#f72d1c', linestyle='--', linewidth=2)
	return (ax)

def plotNetwork (pyinetaObj,Xlim,Ylim,PPcs,PPdq,out_nets=None):
	"""Plot all the networks overlaid in the spectra figure.

	Args:
		pyinetaObj (pyineta object): The pyineta object with the spectra.
		Xlim (tuple): Tuple of (min,max) values of the X axis (13C spectrum); (0,200) for a typical INADEQUATE spectra.
		Ylim (tuple): Tuple of (min,max) values of the Y axis (DQ spectrum); (0,400) for a typical INADEQUATE spectra.
		PPcs (float): ppm threshold to cluster points along the 13C axis.
		PPdq (float): ppm threshold to cluster points along the DQ axis.
		out_nets (str, optional): Filename to save the plot. Defaults to None.
	
	Returns:
		axes handle for the resulting figure
	"""

	(figFin,axsFin)=plotSingle(pyinetaObj,Xlim,Ylim,'Complete plot',grid=False)
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
	line_eqn = lambda x : ((y2-y1)/(x2-x1)) * (x - x1) + y1        
	## generate range of x values based on your graph
	xvals=np.arange(xlim[1],xlim[0],10)
	## plot the line with generate x ranges and created y ranges
	axsFin.plot(xvals, [ line_eqn(x) for x in xvals], color='g', linestyle='-', linewidth=2)
	axsFin=plotAllNet(pyinetaObj.horzPts,pyinetaObj.vertPts,axsFin,max(PPcs,PPdq))
	if out_nets is not None:
		figFin.savefig(out_nets,format='eps',dpi=300)
	return(figFin,axsFin)

def colors(n):
	"""Get a specific color from a list of custom colormap.

	Args:
		n (int): Color number.

	Returns:
		str : Hex code for the specified color.
	"""

	# Custom colormap for visually distinct colors:
	# Colors are: orangered,darkblue,forestgreen,maroon3,fuchsia,cornflower,lime,aqua,moccasin,yellow
	col_list=['#ff4500','#00008b','#228b22','#b03060','#ff00ff','#6495ed','#00ff00','#00ffff','#ffe4b5','#ffff00']
	return col_list[n]

def plotDb (dbNet,ax,col,pos,name):
	"""Generate plots for all networks in the INETA database.

	Args:
		dbNet (dict): The INETA database in json format.
		ax (axes): axes handle.
		col (str): color.
		pos (int): text position.
		name (str): plot title (metabolite name).

	Returns:
		axes handle for the plot.
	"""

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
		ax.annotate(label, xy=(x, y), xytext=(20, -20),fontsize='15',textcoords='offset points', ha='right', va='bottom',color=col)
	ax.text(2,pos*5+5,name_parts[2],fontsize='20',color=col,ha='right')
	return(ax)

def plotMatches (pyinetaObj,outFolder,db_file,Xlim,Ylim):
	"""Plot the matches for a network.

	This function generates a plot for an unknown network with the matched networks overlaid on top.
	Useful for visually assessing the matches identified by pyineta.

	Args:
		pyinetaObj (pyineta object): The pyineta object with the spectra.
		outFolder (str): path to the output folder for the images.
		db_file (dict): The INETA database in json format.
		Xlim (tuple): Tuple of (min,max) values of the X axis (13C spectrum); (0,200) for a typical INADEQUATE spectra.
		Ylim (tuple): Tuple of (min,max) values of the Y axis (DQ spectrum); (0,400) for a typical INADEQUATE spectra.
	"""
	
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
		(figCom,axsCom)=plotSingle(pts,Xlim,Ylim,i[0],netfig,grid=False)
		axsCom.scatter(i[1],i[2], c='0.45', s=100,marker='o',facecolors='none')
		for label, x, y in zip(i[3],i[1],i[2]):	
			axsCom.annotate(label, xy=(x, y), xytext=(-5, 5),fontsize='15',textcoords='offset points', ha='right', va='bottom')
		for j in i[4]:
			axsCom.plot([j[0][0], j[1][0]], [j[0][1], j[1][1]], color='0.45', linestyle='--', linewidth=4)
		if len(i)>5:
			for j in range(5,len(i)):
				axsCom=plotDb(i[j][1][0],axsCom,colors(col),pos,i[j][0])
				for label, x, y in zip(i[3],i[1],i[2]):
					axsCom.annotate(label, xy=(x, y), xytext=(-5, 5),fontsize='15',textcoords='offset points', ha='right', va='bottom')
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

		try:
			X=pyinetaObj.Xlist[len(pyinetaObj.Xlist)-1]
			Y=pyinetaObj.Ylist[len(pyinetaObj.Ylist)-1]
		except AttributeError:
			X=pyinetaObj[0]
			Y=pyinetaObj[1]
		axsCom.plot(X,Y,'k.',markersize=2)
		figCom.savefig(netfig, dpi=300)

def plotIndividualMatch(pyinetaObj,db_file,netNum,id,Xlim,Ylim):
	"""Plot individual network with individual database entries of choice.

	This function allows users to plot any network of their choosing with a database entry.
	This is helpful to generate overlays of any database entry on top of the spectra for visual comparison.

	Args:
		pyinetaObj (pyineta object): The pyineta object with the spectra.
		db_file (dict): The INETA database in json format.
		netNum (str): The network number to plot (Eg: Network1).
		id (str): A string identifier that matches part or all of the Internal IDs used for the INETA database entries.
		Xlim (tuple): Tuple of (min,max) values of the X axis (13C spectrum); (0,200) for a typical INADEQUATE spectra.
		Ylim (tuple): Tuple of (min,max) values of the Y axis (DQ spectrum); (0,400) for a typical INADEQUATE spectra.
	
	Returns:
		figure and axes handle for the plot.
	"""

	# Generate a list of all network names provided
	if ',' in netNum:
		netList=netNum.split(',')
	else:
		netList=[netNum]
	netFirst=netList.pop()
	netfig=netNum+"_"+id+".eps"
	# Plot the first network name match
	for nets in pyinetaObj.NetMatch:
		if netFirst in nets:
			pts=list()
			pts.append(nets[1])
			pts.append(nets[2])
			(figCom,axsCom)=plotSingle(pts,Xlim,Ylim,nets[0],netfig,grid=False)
			axsCom.scatter(nets[1],nets[2], c='0.25', s=100,marker='o',facecolors='none')
			for label, x, y in zip(nets[3],nets[1],nets[2]):	
				axsCom.annotate(label, xy=(x, y), xytext=(-5, 5),fontsize='15',textcoords='offset points', ha='right', va='bottom')
			for j in nets[4]:
				axsCom.plot([j[0][0], j[1][0]], [j[0][1], j[1][1]], color='0.25', linestyle='--', linewidth=4)
	# Plot all other network name matches separated by comma in the -n option
	for n,netSingle in enumerate(netList):
		for nets in pyinetaObj.NetMatch:
			if netSingle in nets:
				c_num=0.45+(n*0.2)
				if (c_num>0.9):
					c_num=0.45
				print(c_num)
				axsCom.scatter(nets[1],nets[2], c=str(c_num), s=100,marker='o',facecolors='none')
				for label, x, y in zip(nets[3],nets[1],nets[2]):	
					axsCom.annotate(label, xy=(x, y), xytext=(-5, 5),fontsize='15',textcoords='offset points', ha='right', va='bottom')
				for j in nets[4]:
					axsCom.plot([j[0][0], j[1][0]], [j[0][1], j[1][1]], color=str(c_num), linestyle='--', linewidth=4)
	# Add title name with all matched network names
	axsCom.set_title(netNum, fontsize=25)
	# Plot all the database network name or id matches
	col=0
	pos=0
	with open(db_file, 'r') as json_file:
		json_db = json.load(json_file)
	for i in json_db:
		if id in i:
			axsCom=plotDb(json_db[i]['Networks'],axsCom,colors(col),pos,i)
			pos+=1
			col+=1
			if col==10:
				col=0
	# Plot all the picked peaks on top of the network and database matches
	try:
		X=pyinetaObj.Xlist[len(pyinetaObj.Xlist)-1]
		Y=pyinetaObj.Ylist[len(pyinetaObj.Ylist)-1]
	except AttributeError:
		X=pyinetaObj[0]
		Y=pyinetaObj[1]
	axsCom.plot(X,Y,'k.',markersize=2)
	# Save final figure
	figCom.savefig(netfig, dpi=300)
	return(figCom,axsCom)