import sys
sys.dont_write_bytecode = True
import numpy as np
np.set_printoptions(precision=3)

from datetime import datetime
import pyineta.picking as picking
import pyineta.filtering as filtering
import pyineta.finding as finding
import pyineta.matching as matching

# Read in the config file and store the parameters

# Step 5: Finalize plotting

class Pyineta:
	
	# Read the file into : 
	#	- the intensity matrix, 
	#	- 13C ppm vector and 
	#	- double quantum ppm vector

	def __init__(self, spectrum):	# Read the ft file using nmrglue
		(self.In,self.Cppm,self.DQppm)=picking.readFt(spectrum)

	@classmethod
	def readMat (cls, spectrum, xmat, ymat):	# Read in numpy arrays
		In=spectrum
		Cppm=xmat
		DQppm=ymat
		return cls(In,Cppm,DQppm)
	
	# Step 1: Peak Picking
	
	def pickPeak (self,PPmin,PPmax,steps,shift=None):
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
		(self.Pts,self.Xlist,self.Ylist,self.P)=picking.pick(self.In,self.Cppm,self.DQppm,PPmin,PPmax,steps)

	# Step 2: Filter points

	def filterPoints (self,PPcs,PPdq):
		self.filteredPts={}
		for k, P in list(self.Pts.items()):
			sortedP=filtering.sortX(P,float(PPcs),0)
			ysortedP=filtering.sortY(sortedP,float(PPdq),1)
			self.filteredPts[k]=filtering.centerMass(ysortedP,"median")

	# Step 3: Find Networks

	def findNetwork (self,levdist,dqt,sumXY,sdt,cst,sel='all'):
		q=0
		
		## Find horizontally aligned peaks
		self.horzPts={}
		if sel=="all":
			self.mergedPts=finding.mergeLevels(self.filteredPts,levdist)
		elif sel=="last":
			self.mergedPts=self.filteredPts[len(self.filteredPts)-1]
		else:
			exit("ERROR: Unknown value for sel. Use either all or last.")
		(self.horzPts)=finding.horzAlign(self.mergedPts,dqt,sumXY,sdt)
		
		## Find vertically aligned peaks from the set of horizontally aligned peaks
		horzMerged = [self.horzPts[k][0] for k in self.horzPts]
		horzMerged.extend([self.horzPts[k][1] for k in self.horzPts])
		self.vertPts=filtering.sortX(horzMerged,cst,0) # 0 for x-axis, 1 of y-axis

		## Build network from aligned peaks
		self.Networks=finding.buildNetwork(self.horzPts,self.vertPts)
		print("Step3.2==> Built ",len(self.Networks)," networks.")
		# print("BLABLABLA=+++",len(self.Networks))
		# print(Netfull)
		# print(self.Networks)
		## Generate a list of all connected pairs of peaks.
		self.Pairs=finding.listPairs(self.horzPts,self.vertPts)

	def writeNetwork (self,net_file):

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
		import json

		self.NetTag=[]
		self.NetMatch=[]
		TagDict={}
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
		out_file = open(outFolder+"/"+match_file, 'w')
		header="#NetworkNum\tID\tMatchName\tSolvent\tHitscore\tCoverageScore\tMatchedConnections\tUnmatchedConnections\n"
		out_file.write(header)
		for i in self.NetMatch:
			outstr=str()
			if len(i)>5:
				for j in range(5,len(i)):
					name=i[j][0]
					name_parts=name.split('::')
					outline=i[0]+'\t'+name_parts[1]+'\t'+name_parts[2]+'\t'+name_parts[4]+'\t'+str(i[j][4])+'\t'+str(i[j][5])+"\t"
					for k in i[j][2]:
						outline+=k+'->'+str(i[j][2][k][0])+'-'+str(i[j][2][k][1])+','
					outline+='\t'
					for k in i[j][3]:
						outline+=k+'->'+str(i[j][3][k][0])+'-'+str(i[j][3][k][1])+','
					outline+='\n'
					outstr+=outline
			out_file.write(outstr)

	def summarize (self,outFolder,filename):
		out_file = open(outFolder+"/"+filename, 'w')
		now = datetime.now() # datetime object containing current date and time
		outstr="Summary for PyINETA run on "+str(now)+" :\n"
		out_file.write(outstr)
		try:

			# How many peaks picked
			ct_peak1=sum([len(self.Pts[x]) for x in self.Pts if isinstance(self.Pts[x], list)])
			outstr="# of picked peaks in all steps: "+str(ct_peak1)+"\n"
			out_file.write(outstr)

			# How many total peaks after filtering
			ct_peak2=np.vstack(list(self.filteredPts.values())).shape[0]
			outstr="# of peaks after filtering: "+str(ct_peak2)+"\n"
			out_file.write(outstr)

			# How many filtered points after clustering close pts from steps
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
	import configparser
	config = configparser.ConfigParser()
	config.read("config.ini")
	param=dict()
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
	
	param["PPCS"]= config.getfloat("FilterPoints", "PPCS")
	param["PPDQ"]= config.getfloat("FilterPoints", "PPDQ")
	param["OutImage_filter_separate"]= config.get("FilterPoints", "OutImage_filter_separate")
	param["OutImage_filter_complete"]= config.get("FilterPoints", "OutImage_filter_complete")

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
	param["DQMT"]= config.getfloat("MatchDatabase", "DQMT")
	param["Topology_tolerance"]= config.getfloat("MatchDatabase", "Topology_tolerance")
	param["Hit_Score_threshold"]= config.getfloat("MatchDatabase", "Hit_Score_threshold")
	param["Coverage_Score_threshold"]= config.getfloat("MatchDatabase", "Coverage_Score_threshold")
	param["Matches_list_output_file"]= config.get("MatchDatabase", "Matches_list_output_file")
	param["Summary_file"]= config.get("MatchDatabase", "Summary_file")

	# config.optionxform=str
	# config.read(configFile)
	# general=dict(config.items('General'))
	# peakPick=dict(config.items('PeakPick'))
	# filterPoints=dict(config.items('FilterPoints'))
	# findNetwork=dict(config.items('FindNetwork'))
	# matchDatabase=dict(config.items('MatchDatabase'))
	# return(general,peakPick,filterPoints,findNetwork,matchDatabase)
	return(param)

def stepError ():
	print("ERROR:  Cound not find all required attributes in saved pickle file.")
	print("\tMake sure you've run all previous steps.")
	print("\tPyineta run not finished.")
	sys.exit(0)
# def writeFiles ():


