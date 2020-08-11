#!/usr/bin/python

#import all dependencies

# create output folder

# check if config file exists and is properly formatted

#Open spectra and run pickPeak

import pyineta.pyineta as pyineta
import pyineta.plotting as plotting
import pyineta.picking as picking
# import pyineta.ineta_functions_v1 as F
import matplotlib.pyplot as plt

def main(args):
	print(args.output)
	print("Step0.1==> Reading Config file...", end=" ")
	param=pyineta.readConfig(args.configfile)
	Xrng=(param['Xrange_min'],param['Xrange_max'])
	Yrng=(param['Yrange_min'],param['Yrange_max'])
	print("..Done.")

	print("Step0.2==> Loading NMR file...", end=' ')
	try:
		spec=pyineta.Pyineta(param['Ft_File'])
	except FileNotFoundError:
		print("\tLoading from data matrix file.")
		Intsy = np.loadtxt(param['Data_Matrix_File'],dtype=np.float32)
		Xs = np.loadtxt(param['13C_Ppm_File'],dtype=np.float32)
		Ys = np.loadtxt(param['Double_Quantum_File'], dtype=np.float32)
		spec=pyineta.Pyineta(Intsy,Xs,Ys)
	print("..Done.")

	# Load previously saved variables if running select steps
	# Use pickle to save the entire class
	# self.Pts
	# self.Clist

	# Step 1: Peak Picking

	if param["Shift"].lower() == "yes":
		padUnits=[param["Shift13C"],param["Full13C"],param["FullDQ"],param["Direction"]]
	else:
		padUnits=None
	PPmin=param['PPmin']
	PPmax=param['PPmax']
	steps=param['steps']
	print("Step1==> Peak picking with cutoffs min =",PPmin,"and max =",PPmax, "with ",steps," steps...", end=' ')
	spec.pickPeak(PPmin,PPmax,steps,padUnits)
	print("..Done.")
	# print(spec.P)

	#Saving variables

	# Step 2: Filter points

	PPcs=param['PPCS']
	PPdq=param['PPDQ']
	print("Step2.1==> Filtering points using PPCS=",PPcs,"and PPDQ=",PPdq,"and finding center of mass for clustered points...", end=' ')
	spec.filterPoints(PPcs,PPdq)
	print("..Done.")
	# print(spec.Clist)

	#plotting
	#Saving variables

	# Step 3: Find Networks

	levdist=param['LevelPointsDistance']
	dqt=param['DQT']
	sumXY=param['SumXY']
	sdt=param['SDT']
	cst=param['CST']
	print("Step3==>BLABALBALBALB")
	spec.findNetwork(levdist,dqt,sumXY,sdt,cst)

	#
	net_file=param['Network_output_file']
	spec.writeNetwork(net_file)
	
	# No plotting
	# Saving variables

	# Step 4: Match Database

	inetaDb=param['Database_file']
	ambig=param['Ambiguity']
	near=param['CSMT']
	match=param['Match_tolerance']
	topo=param['Topology_tolerance']
	hitSc=param['Hit_Score_threshold']
	covSc=param['Coverage_Score_threshold']
	spec.matchDb(inetaDb,ambig,near,match,topo,hitSc,covSc)

	# Initialize colors and output folders.
	# Plotting
	# Saving variables
	match_file= param['Matches_list_output_file']
	spec.writeMatches(args.output,match_file)
	
	#plotting for all

	out_sep= param['OutImage_pick_separate']
	out_comp= param['OutImage_pick_complete']
	out_sep2= param['OutImage_filter_separate']
	out_comp2= param['OutImage_filter_complete']
	out_nets=param['OutImage_network_AllNets']

	gen_fig=param['Generate_figure']
	print("Generate_figure option set to:",gen_fig)
	if (gen_fig == "Yes"):
		print("Step5.0==> Generating Plots...")
		
		# Plot 1
		label=[]
		for i in picking.frange(PPmin,PPmax,steps):
			# print(i,PPmax,PPmin)
			label.append(i)
		
		if (len(spec.Xlist)>1):
			(figSep, axsSep)=plotting.plotFigSep(spec,Xrng,Yrng,label,out_sep) # Separate
			print(type(axsSep),axsSep)
			if (axsSep.ndim<2):
				axsSep = np.reshape(axsSep, (-1, 2))

		(figCom,axsCom)=plotting.plotSingle(spec,Xrng,Yrng,'Complete plot',out_comp) # Complete

		# Plot 2
		plotting.plotFilteredPoints(figSep,axsSep,figCom,axsCom,spec,PPcs,PPdq,out_sep2,out_comp2)
		
		# Plot 3
		plotting.plotNetwork(spec,Xrng,Yrng,PPcs,PPdq,out_nets)

	# 	#Plot 4
		plotting.plotMatches(spec,args.output,inetaDb,Xrng,Yrng)

	else:
		print("Skipping Plot Generation")

if __name__ == '__main__':

	import argparse
	import os

	parser = argparse.ArgumentParser(description="Script to run the INETA pipeline.")
	parser.add_argument('-c', '--configfile', required=True,
		help='Required: tab delimited Template table listing the details of the new functions. Refer to template_sheet.tsv')
	parser.add_argument('-o', '--output', default=os.getcwd(),
		help='Optional: Full path to the output folder.(Default: Current folder)')
	# parser.add_argument('-w', '--wiki', help='Optional: Full path to the Edison lab metabolomics wiki folder. Used to check if a wiki page already exists before writing.', required=False)

	# parser.add_argument('-q', dest='pos', action='store',
	# 	 type=str,
	# 	 help='Query node to find neighbors')
	# parser.add_argument('-g', dest='pttrn', action='store',
	# 	 type=FileType('r'),
	# 	 help='List of patterns to generate network from/Alternative network')
	args = parser.parse_args()
	main(args)


