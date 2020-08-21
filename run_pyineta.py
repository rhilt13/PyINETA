#!/usr/bin/python

# Setup a download object
# If its there okay, otherwise download it
# If the downloaded file is more than 1 month old, re download it.
# Compress and make it available in place But uncompress it to run.

#import all dependencies

# create output folder

# check if config file exists and is properly formatted

#Open spectra and run pickPeak

import argparse
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pyineta.pyineta as pyineta
import pyineta.plotting as plotting
import pyineta.picking as picking


def main(args):

	## Read configuration file

	print("Step0.1==> Reading Config file...", end=" ")
	param=pyineta.readConfig(args.configfile)
	Xrng=(param['Xrange_min'],param['Xrange_max'])
	Yrng=(param['Yrange_min'],param['Yrange_max'])
	print("..Done.")
		
	## Read the NMR Ft or data matrix file
	if args.steps.lower() in {'all','load','load+'}:
		print("Step0.2==> Loading NMR file...", end=' ')
		try:
			spec=pyineta.Pyineta(param['Ft_File'])
		except FileNotFoundError:
			print("\nTrying to load from data matrix file.")
			try:
				Intsy = np.loadtxt(param['Data_Matrix_File'],dtype=np.float32)
				Xs = np.loadtxt(param['13C_Ppm_File'],dtype=np.float32)
				Ys = np.loadtxt(param['Double_Quantum_File'], dtype=np.float32)
				spec=pyineta.Pyineta(Intsy,Xs,Ys)
			except:
				print("ERROR: No input files found as specified in the config file.")
				print("\tPlease check the [Ft_File] or the [Data_Matrix_File, 13C_Ppm_File and Double_Quantum_File] options in the config file.")
				exit(0)
		else:
			pickle.dump(spec, open('ptf_pyINETAObj.tmp', 'wb'))
			print("..Done.")

	## Check for the presence of the pyineta object pickle file to load data from for subsequent steps.
	
	if args.steps.lower() not in {'all','load','load+'}:
		if not os.path.isfile("ptf_pyINETAObj.tmp"):
			print("ERROR: PyINETA object file not found. Please run the load module with -s load option first.")
			exit(0)
	
	## Step 1: Peak Picking

	PPmin=param['PPmin']
	PPmax=param['PPmax']
	steps=param['steps']

	if args.steps.lower() in {'all','pick','load+','pick+'}:

		if args.steps.lower() in {'pick','pick+'}:
			with open('ptf_pyINETAObj.tmp', 'rb') as handle:
				spec = pickle.loads(handle.read())

		if param["Shift"].lower() == "yes":
			padUnits=[param["Shift13C"],param["Full13C"],param["FullDQ"],param["Direction"]]
		else:
			padUnits=None
		print("Step1==> Peak picking with cutoffs min =","{:.2e}".format(PPmin),"and max =","{:.2e}".format(PPmax), "with ",steps," steps...", end=' ')
		try:
			spec.pickPeak(PPmin,PPmax,steps,padUnits)
		except AttributeError:
			pyineta.stepError()
		# print(spec.Pts)
		pickle.dump(spec, open('ptf_pyINETAObj.tmp', 'wb'))
		print("..Done.")

	## Step 2: Filter points

	PPcs=param['PPCS']
	PPdq=param['PPDQ']

	if args.steps.lower() in {'all','filter','load+','pick+','filter+'}:

		if args.steps.lower() in {'filter','filter+'}:
			with open('ptf_pyINETAObj.tmp', 'rb') as handle:
				spec = pickle.loads(handle.read())

		
		print("Step2==> Filtering points using PPCS=",PPcs,"and PPDQ=",PPdq,"and finding center of mass for clustered points...", end=' ')
		try:
			spec.filterPoints(PPcs,PPdq)
		except AttributeError:
			pyineta.stepError()
		pickle.dump(spec, open('ptf_pyINETAObj.tmp', 'wb'))
		print("..Done.")

	## Step 3: Find Networks
	
	if args.steps.lower() in {'all','find','load+','pick+','filter+','find+'}:

		if args.steps.lower() in {'find','find+'}:
			with open('ptf_pyINETAObj.tmp', 'rb') as handle:
				spec = pickle.loads(handle.read())

		levdist=param['LevelPointsDistance']
		dqt=param['DQT']
		sumXY=param['SumXY']
		sdt=param['SDT']
		cst=param['CST']
		print("Step3==> Finding networks...")
		try:
			spec.findNetwork(levdist,dqt,sumXY,sdt,cst)
		except AttributeError:
			pyineta.stepError()
		#
		net_file=param['Network_output_file']
		spec.writeNetwork(net_file)
		pickle.dump(spec, open('ptf_pyINETAObj.tmp', 'wb'))
		print("...Done.")

	## Step 4: Match Database

	inetaDb=param['Database_file']
	
	if args.steps.lower() in {'all','match','load+','pick+','filter+','find+','match+'}:

		if args.steps.lower() in {'match','match+'}:
			with open('ptf_pyINETAObj.tmp', 'rb') as handle:
				spec = pickle.loads(handle.read())

		ambig=param['Ambiguity']
		near=param['CSMT']
		match=param['Match_tolerance']
		topo=param['Topology_tolerance']
		hitSc=param['Hit_Score_threshold']
		covSc=param['Coverage_Score_threshold']
		print("Step4==> Searching database",inetaDb,"for matches...")
		try:
			spec.matchDb(inetaDb,ambig,near,match,topo,hitSc,covSc)
		except AttributeError:
			pyineta.stepError()
		#
		match_file= param['Matches_list_output_file']
		spec.writeMatches(args.output,match_file)
		pickle.dump(spec, open('ptf_pyINETAObj.tmp', 'wb'))
		print("...Done.")

	if args.steps.lower() in {'all','summary','load+','pick+','filter+','find+','match+'}:
		
		if args.steps.lower() in {'summary'}:
			with open('ptf_pyINETAObj.tmp', 'rb') as handle:
				spec = pickle.loads(handle.read())

		summary_file=param['Summary_file']
		spec.summarize(args.output,summary_file)
	#plotting for all

	out_sep= param['OutImage_pick_separate']
	out_comp= param['OutImage_pick_complete']
	out_sep2= param['OutImage_filter_separate']
	out_comp2= param['OutImage_filter_complete']
	out_nets=param['OutImage_network_AllNets']

	if args.steps.lower() not in {'load','summary'}:
		print("Generate_figure option set to:",args.figure)
	
		if args.figure.lower() == "yes":
			print("Step5==> Generating Plots...", end=' ')
		
			# Plot 1
			if (args.steps.lower() in {'all','pick','filter','load+','pick+','filter+'}):
				label=[]
				for i in picking.frange(PPmin,PPmax,steps):
					# print(i,PPmax,PPmin)
					label.append(i)
				if (len(spec.Xlist)>1):
					(figSep, axsSep)=plotting.plotFigSep(spec,Xrng,Yrng,label,out_sep) # Separate
					if (axsSep.ndim<2):
						axsSep = np.reshape(axsSep, (-1, 2))
				(figCom,axsCom)=plotting.plotSingle(spec,Xrng,Yrng,'Complete plot',out_comp) # Complete

			# Plot 2
			if (args.steps.lower() in {'all','filter','load+','pick+','filter+'}):
				plotting.plotFilteredPoints(figSep,axsSep,figCom,axsCom,spec,PPcs,PPdq,out_sep2,out_comp2)
			
			# Plot 3
			if (args.steps.lower() in {'all','find','load+','pick+','filter+','find+'}):
				plotting.plotNetwork(spec,Xrng,Yrng,PPcs,PPdq,out_nets)

			#Plot 4
			if (args.steps.lower() in {'all','match','load+','pick+','filter+','find+','match+'}):
				plotting.plotMatches(spec,args.output,inetaDb,Xrng,Yrng)
			print("...Done.")
		
		else:
			print("Skipping Plot Generation")

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="Script to run the INETA pipeline.")
	parser.add_argument('-c', '--configfile', required=True,
		help='Required: A config file with all the options and parameters required for the INETA run.')
	parser.add_argument('-o', '--output', default=os.getcwd(),
		help='Optional: Full path to the output folder.(Default: Current folder)')
	parser.add_argument('-s', '--steps', default="all",
		help='Optional: Specify which steps you want to run. Can be one of {all,load,pick,filter,find,match,summary,load+,pick+,filter+,find+,match+}. Adding a + to the end of option runs all steps after as well.')
	parser.add_argument('-f', '--figure', default="yes",
		help='Optional: Generate figures- yes or no. Default: Yes')
	args = parser.parse_args()
	main(args)