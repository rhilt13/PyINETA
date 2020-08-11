import sys
sys.dont_write_bytecode = True

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
	def readMat (cls, spectrum, xmat, ymat):	# Read in python lists/matrices
		In=spectrum
		Cppm=xmat
		DQppm=ymat
		return cls(In,Cppm,DQppm)
	
	# Step 1: Peak Picking
	
	def pickPeak (self,PPmin,PPmax,steps,shift=None):
		if shift is not None:
			print(shift,type(shift),len(shift))
			print("MUAHAHAHAHHAHAHA")
			if type(shift) is list and len(shift)==4:
				print("MUAHAHAHAHHAHAHA")
				self.In=picking.shifting(self.In,*shift)
			else:
				exit("ERROR: Argument shift needs to be a list with 4 items:padding units, total size of X axis, total size of Y axis and direction (either pos or neg). Eg: [20,4096,8192,'pos']")
		(self.Pts,self.Xlist,self.Ylist,self.P)=picking.pick(self.In,self.Cppm,self.DQppm,float(PPmin),float(PPmax),int(steps))
		return(self.Pts,self.Xlist,self.Ylist,self.P)
	
	# Step 2: Filter points

	def filterPoints (self,PPcs,PPdq):
		self.filteredPts={}
		for k, P in list(self.Pts.items()):
			sortedP=filtering.Sort(P,float(PPcs),0)
			ysortedP=filtering.YSort(sortedP,float(PPdq),1)
			self.filteredPts[k]=filtering.CenterMass(ysortedP,"median")

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
		self.vertPts=filtering.Sort(horzMerged,cst,0) # 0 for x-axis, 1 of y-axis
		print("Step3.3==> Identified vertically aligned peaks.")

		## Build network from aligned peaks
		(self.Networks,Netfull,Netct)=finding.BuildNetwork(self.horzPts,self.vertPts)
		print("Step3.4==> Built ",Netct," networks.")
		
		## Generate a list of all connected pairs of peaks.
		self.Pairs=finding.listPairs(self.horzPts,self.vertPts)

	def writeNetwork (self,net_file):
		import numpy as np

		out_file = open(net_file, 'w')
		j=1
		for i in self.Networks:
			outstr=''.join(['Network',str(j),' = ',str(np.around(np.asarray(i),decimals=2).tolist()),'\n'])
			j+=1
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
			# print(unknownConn)
		#	print "---MATCHING---"
			# print Ptags
			with open(inetaDb, 'r') as json_file:
				json_data = json.load(json_file)
			out=matching.MatchDatabase(json_data,Pvals,Ptags,unknownConn,ambig,near,match,topo,hitSc,covSc,q,FinalPairs)
			# print "out==>",len(out),out
			self.NetMatch.append(out)
			self.NetTag.append(Ptags)
			if (len(out)>5):
				ct_match+=1;
				print("Step4.2==> Match found for Network",q)
			else:
				print("Step4.2==> No matches found for Network",q)
		print("Matches found for ", ct_match, "Networks out of ",len(self.Networks))

	def writeMatches (self,outFolder,match_file):
		out_file = open(outFolder+"/"+match_file, 'w')

		for i in self.NetMatch:
			outstr=str()
			if len(i)>5:
				for j in range(5,len(i)):
					name=i[j][0]
					name_parts=name.split('::')
					outline=i[0]+','+name_parts[1]+','+name_parts[2]+','+name_parts[4]+',HitScore:'+str(i[j][4])+',CoverageScore:'+str(i[j][5])
					outline+="\n\tMatchedConnections: "
					for k in i[j][2]:
						outline+=k+'->'+str(i[j][2][k][0])+'-'+str(i[j][2][k][1])+','
					outline+='\n\tUnMatchedConnections: '
					for k in i[j][3]:
						outline+=k+'->'+str(i[j][3][k][0])+'-'+str(i[j][3][k][1])+','
					outline+='\n'
					outstr+=outline
			out_file.write(outstr)

def readConfig (configFile):
	import configparser
	config = configparser.ConfigParser()
	config.read("config.ini")
	param=dict()
	param["Scripts_location"]=config.get("General","Scripts_location")
	param["Generate_figure"]= config.get("General", "Generate_figure")
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
	param["Shift13C"]=config.getfloat("PeakPick","Shift13C")
	param["Full13C"]=config.getfloat("PeakPick","Full13C")
	param["FullDQ"]=config.getfloat("PeakPick","FullDQ")
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


	# config.optionxform=str
	# config.read(configFile)
	# general=dict(config.items('General'))
	# peakPick=dict(config.items('PeakPick'))
	# filterPoints=dict(config.items('FilterPoints'))
	# findNetwork=dict(config.items('FindNetwork'))
	# matchDatabase=dict(config.items('MatchDatabase'))
	# return(general,peakPick,filterPoints,findNetwork,matchDatabase)
	return(param)

# def writeFiles ():


