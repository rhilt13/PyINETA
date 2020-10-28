"""The core classes and functions for pyINETA

This file includes the core Pyineta class and function definitions.
Includes the following:
	* Pyineta class - The core Pyineta class.
		> pickPeak - Peak picking the input spectra.
		> clusterPoints - Clustering a list of closely located points.
		> findNetwork - Find INETA networks in the Pyineta object spectrum.
		> writeNetwork - Write the found networks along with their points to a file.
		> matchDb - Match the found INETA networks to the database entries.
		> writeMatches - Write the matches report to a tab-delimited file.
		> summarize - Write a summary file with the main reports for the entire INETA run.
	* readConfig - Read the parameters and options from the config file.
	* stepError - Report error message.
"""

import sys
sys.dont_write_bytecode = True
import numpy as np
np.set_printoptions(precision=3)

from datetime import datetime
import pyineta.picking as picking
import pyineta.clustering as clustering
import pyineta.finding as finding
import pyineta.matching as matching

class Pyineta:
	"""The core Pyineta class.

	Attributes:
		In (ndarray): A numpy array with the intensities.
		Cppm (1D-array): A 1D array with the 13C ppm values.
		DQppm (1D-array): A 1D array with the double quantum ppm values. 
		Pts (dict): A dict mapping an array of points (x,y) to the iteration number it was found in.
        Xlist (dict) : A dict mapping an array of x axis values (13C ppm values) to the iteration number it was found in.
        Ylist (dict) : A dict mapping an array of y axis values (DQ ppm values) to the iteration number it was found in.
		clusteredPts (ndarray) : An array of center of masses (peak centers) for all the clusters.
		mergedPts (ndarray) : A 2D array of merged points.
		horzPts (dict): A dict mapping horizontally aligned peaks to their indices.
		vertPts (dict): A dict mapping lists of clustered points to their respective cluster numbers.
		Networks (list) : list of lists with points belonging to a network.
		Pairs (list) : A list of all horizontally connected pairs of points included in a network.
		NetTag (list) : List of Tags for naming the unknown peaks in a network.
		NetMatch (list): A list of networks with their corresponding hits and hit details.

	"""

	def __init__(self, spectrum):	# Read the ft file using nmrglue
		"""The __init__ method.

		Initialize the Pyineta object.

		Args:
			spectrum (str): input ft filename with the pre-processed NMR spectra.
		"""

		(self.In,self.Cppm,self.DQppm)=picking.readFt(spectrum)

	@classmethod
	def readMat (cls, spectrum, xmat, ymat):	# Read in numpy arrays
		"""Read the spectra from matrices instead of an ft file.

		This method allows defining the Pyineta class using 3 separate matrices.

		Args:
			spectrum (ndarray): A numpy array with the intensities.
			xmat (1D-array): A 1D array with the 13C ppm values.
			ymat (1D-array): A 1D array with the double quantum ppm values.

		Returns:
			pyineta object : The pyineta object with the spectra.
		"""

		In=spectrum
		Cppm=xmat
		DQppm=ymat
		return cls(In,Cppm,DQppm)
	
	# Step 1: Peak Picking
	
	def pickPeak (self,PPmin,PPmax,steps,shift=None):
		"""Peak picking the input spectra.

		Args:
			PPmin (float): Minimum intensity value to be considered a peak.
        	PPmax (float): Maximum intensity value to be considered a peak.
			steps (int): Number of iterations to find peaks within the PPmin to PPmax range.
			shift (list, optional): A list with 4 values: units to shift, length of 13C and DQ axes and the direction to shift. Defaults to None.
		"""

		if shift is not None:
			if type(shift) is list and len(shift)==4:
				if (shift[3] == 'pos'):
					directionShift="higher"
				else:
					directionShift="lower"
				ppmUnits=(200/shift[1])*shift[0]
				print("\nStep1.1==>Shifting spectra to %s ppm level by %.2f units (%.2f ppm on 13C axis)" % (directionShift,shift[0],ppmUnits))
				self.In=picking.shifting(self.In,*shift)
			else:
				exit("ERROR: Argument shift needs to be a list with 4 items:padding units, total size of X axis, total size of Y axis and direction (either pos or neg). Eg: [20,4096,8192,'pos']")
		(self.Pts,self.Xlist,self.Ylist)=picking.pick(self.In,self.Cppm,self.DQppm,PPmin,PPmax,steps)

	# Step 2: Cluster points

	def clusterPoints (self,PPcs,PPdq):
		"""Clustering a list of closely located points.

		Args:
			PPcs (float): Points found within this threshold along the 13C axis are grouped into a single cluster.
			PPdq (float): A distance threshold for splitting points into a different cluster along the DQ axis.
		"""

		self.clusteredPts={}
		for k, P in list(self.Pts.items()):
			sortedP=clustering.gather(P,float(PPcs),0)
			ysortedP=clustering.splitY(sortedP,float(PPdq),1)
			self.clusteredPts[k]=clustering.centerMass(ysortedP,"median")

	# Step 3: Find Networks

	def findNetwork (self,levdist,dqt,sumXY,sdt,cst,sel='all'):
		"""Find INETA networks in the Pyineta object spectrum.

		Args:
			levdist (float): Threshold for merging cluster centers from different iterations into a single point.
			dqt (float): Threshold to split the merged points along the DQ axis.
			sumXY (float): Threshold to check if the sum of 13C ppm values for two points equals their DQ ppm or not.
			sdt (float): Threshold to check of the 2 points are equidistant from the diagonal.
			cst (float): Points found within this threshold are grouped into a single cluster.
			sel (str, optional): 'all' or 'last'. Defaults to 'all'.
		"""

		## Find horizontally aligned peaks
		if sel=="all":
			self.mergedPts=finding.mergeLevels(self.clusteredPts,levdist)
		elif sel=="last":
			self.mergedPts=self.clusteredPts[len(self.clusteredPts)-1]
		else:
			exit("ERROR: Unknown value for sel. Use either all or last.")
		(self.horzPts)=finding.horzAlign(self.mergedPts,dqt,sumXY,sdt)
		
		## Find vertically aligned peaks from the set of horizontally aligned peaks
		horzMerged = [self.horzPts[k][0] for k in self.horzPts]
		horzMerged.extend([self.horzPts[k][1] for k in self.horzPts])
		self.vertPts=clustering.gather(horzMerged,cst,0) # 0 for x-axis, 1 of y-axis

		## Build network from aligned peaks
		self.Networks=finding.buildNetwork(self.horzPts,self.vertPts)
		print("Step3.2==> Built ",len(self.Networks)," networks.")

		## Generate a list of all connected pairs of peaks.
		self.Pairs=finding.listPairs(self.horzPts,self.vertPts)

	def writeNetwork (self,net_file):
		""" Write the found networks along with their points to a file.

		Args:
			net_file (str): Network output filename.
		"""

		out_file = open(net_file, 'w')
		j=1
		for i in self.Networks:
			arr=np.around(np.asarray(i),decimals=2)
			outstr=''.join(['Network',str(j),'\t'])
			j+=1
			m=1
			for k in arr:
				outstr=outstr+str(k)
				if m<len(arr):
					outstr=outstr+','
				else:
					outstr=outstr+"\n"
				m+=1
			out_file.write(outstr)

	# Step 4: Match Database

	def matchDb (self,inetaDb,ambig,near,match,topo,hitSc,covSc):
		""" Match the found INETA networks to the database entries.

		Args:
			inetaDb (dict): The INETA database in json format.
			ambig (float): Amibiguity tolerance (Removes all database entries with ambiguity higher than this tolerance).
			near (float): Tolerance for distance in the 13C dimension between unknown peak and a match database peak.
			match (int): Tolerance for number of matches in a single network.
			topo (float): Topology tolerance (How far the point can be in all directions from the match peak point).
			hitSc (float): Hit score threshold.
			covSc (float): Coverage score threshold.
		"""

		import json

		self.NetTag=[]
		self.NetMatch=[]
		q=0
		ct_match=0
		for Pvals in self.Networks:
			q=q+1
			(Ptags,unknownConn,FinalPairs)=matching.prepUnknowns(Pvals,self.Pairs)
			with open(inetaDb, 'r') as json_file:
				json_data = json.load(json_file)
			out=matching.matchDatabase(json_data,Pvals,Ptags,unknownConn,ambig,near,match,topo,hitSc,covSc,q,FinalPairs)
			self.NetMatch.append(out)
			self.NetTag.append(Ptags)
			if (len(out)>5):
				ct_match+=1
				print("Step4.1==> Match found for Network",q)
			else:
				print("Step4.1==> No matches found for Network",q)
		print("Matches found for ", ct_match, "Networks out of ",len(self.Networks))

	def writeMatches (self,outFolder,match_file):
		"""Write the matches report to a tab-delimited file.

		Args:
			outFolder (str): Path to the output folder.
			match_file (str): Output filename to write the matches.
		"""

		out_file = open(outFolder+"/"+match_file, 'w')
		header="#NetworkNum\tID\tMatchName\tSolvent\tAmbiguityScore\tHitscore\tCoverageScore\tMatchedConnections\tUnmatchedConnections\n"
		out_file.write(header)
		for i in self.NetMatch:
			outstr=str()
			if len(i)>5:
				for j in range(5,len(i)):
					name=i[j][0]
					name_parts=name.split('::')
					outline=i[0]+'\t'+name_parts[1]+'\t'+name_parts[2]+'\t'+name_parts[4]+'\t'+str(i[j][4])+'\t'+str(i[j][5])+'\t'+str(i[j][6])+"\t"
					for k in i[j][2]:
						outline+=k+'->'+str(i[j][2][k][0])+'-'+str(i[j][2][k][1])+','
					outline+='\t'
					for k in i[j][3]:
						outline+=k+'->'+str(i[j][3][k][0])+'-'+str(i[j][3][k][1])+','
					outline+='\n'
					outstr+=outline
			out_file.write(outstr)

	def summarize (self,outFolder,filename):
		"""Write a summary file with the main reports for the entire INETA run.

		Args:
			outFolder (str): Path to the output folder.
			filename (str): Output filename to write the report.
		"""

		out_file = open(outFolder+"/"+filename, 'w')
		now = datetime.now() # datetime object containing current date and time
		outstr="Summary for PyINETA run on "+str(now)+" :\n"
		out_file.write(outstr)
		try:

			# How many peaks picked
			ct_peak1=sum([len(self.Pts[x]) for x in self.Pts if isinstance(self.Pts[x], list)])
			outstr="# of picked peaks in all steps: "+str(ct_peak1)+"\n"
			out_file.write(outstr)

			# How many total peaks after clustering
			ct_peak2=np.vstack(list(self.clusteredPts.values())).shape[0]
			outstr="# of peaks after clustering: "+str(ct_peak2)+"\n"
			out_file.write(outstr)

			# How many clustered points after clustering close pts from steps
			outstr="# of peaks retained after merging steps: "+str(self.mergedPts.shape[0])+"\n"
			out_file.write(outstr)

			# How many networks
			outstr="# of networks found: "+str(len(self.Networks))+"\n"
			out_file.write(outstr)
		
			# How many matches
			ct_match=0
			for i in self.NetMatch:
				if (len(i)>5):
					ct_match+=1
			outstr="# of matches found: "+str(ct_match)+"\n"
			out_file.write(outstr)
		
		except:
			pass

def readConfig (configFile):
	"""Read the parameters and options from the config file.

	Args:
		configFile (str): The config filename. 

	Returns:
		dict : A dict mapping the parameters to their values.
	"""

	import configparser
	config = configparser.ConfigParser()
	config.read("config.ini")
	param=dict()
	try:
		param["Ft_File"]= config.get("PeakPick", "Ft_File")
		param["Data_Matrix_File"]= config.get("PeakPick", "Data_Matrix_File")
		param["13C_Ppm_File"]= config.get("PeakPick", "13C_Ppm_File")
		param["Double_Quantum_File"]= config.get("PeakPick", "Double_Quantum_File")
		param["Xrange_min"]= config.getint("PeakPick", "Xrange_min")
		param["Xrange_max"]= config.getint("PeakPick", "Xrange_max")
		param["Yrange_min"]= config.getint("PeakPick", "Yrange_min")
		param["Yrange_max"]= config.getint("PeakPick", "Yrange_max")
		param["OutImage_pick_separate"]= config.get("PeakPick", "OutImage_pick_separate")
		param["OutImage_pick_complete"]= config.get("PeakPick", "OutImage_pick_complete")
		param["Shift"]=config.get("PeakPick","Shift")
		param["Direction"]=config.get("PeakPick","Direction")
		param["Shift13C"]=config.getint("PeakPick","Shift13C")
		param["Full13C"]=config.getint("PeakPick","Full13C")
		param["FullDQ"]=config.getint("PeakPick","FullDQ")
		param["PPmin"]= config.getfloat("PeakPick", "PPmin")
		param["PPmax"]= config.getfloat("PeakPick", "PPmax")
		param["steps"]= config.getint("PeakPick", "steps")
		
		param["PPCS"]= config.getfloat("ClusterPoints", "PPCS")
		param["PPDQ"]= config.getfloat("ClusterPoints", "PPDQ")
		param["OutImage_cluster_separate"]= config.get("ClusterPoints", "OutImage_cluster_separate")
		param["OutImage_cluster_complete"]= config.get("ClusterPoints", "OutImage_cluster_complete")

		param["Select"]= config.get("FindNetwork", "Select")
		param["LevelPointsDistance"]= config.getfloat("FindNetwork", "LevelPointsDistance")
		param["DQT"]= config.getfloat("FindNetwork", "DQT")
		param["SumXY"]= config.getfloat("FindNetwork", "SumXY")
		param["SDT"]= config.getfloat("FindNetwork", "SDT")
		param["CST"]= config.getfloat("FindNetwork", "CST")
		param["Network_output_file"]= config.get("FindNetwork", "Network_output_file")
		param["OutImage_network_AllNets"]= config.get("FindNetwork", "OutImage_network_AllNets")

		param["Database_file"]= config.get("MatchDatabase", "Database_file")
		param["Ambiguity"]= config.getfloat("MatchDatabase", "Ambiguity")
		param["CSMT"]= config.getfloat("MatchDatabase", "CSMT")
		param["Match_tolerance"]= config.getint("MatchDatabase", "Match_tolerance")
		param["Topology_tolerance"]= config.getfloat("MatchDatabase", "Topology_tolerance")
		param["Hit_Score_threshold"]= config.getfloat("MatchDatabase", "Hit_Score_threshold")
		param["Coverage_Score_threshold"]= config.getfloat("MatchDatabase", "Coverage_Score_threshold")
		param["Matches_list_output_file"]= config.get("MatchDatabase", "Matches_list_output_file")
		param["Summary_file"]= config.get("MatchDatabase", "Summary_file")

	except:
		print("\nERROR: Encountered problems with the config file.")
		print("\tDo you have an older version?")
		print("\tIt should match template config file provided.")
		print("### ERROR MESSAGE ###")
		raise
	return(param)

def stepError (errmsg):
	"""Report error message.

	Args:
		errmsg (str): Error message reported by python interpreter.
	"""

	print("ERROR:  %s" %(errmsg))
	print("\tCound not find all required attributes in saved pickle file.")
	print("\tMake sure you've run all previous steps.")
	print("\tPyineta run not finished.")
	sys.exit(0)