#! /usr/bin/env python

import os
import re
import json
from random import *
import pynmrstar
from argparse import ArgumentParser
from argparse import FileType

def list_files(directory, extension):
    filelist=list(f for f in os.listdir(directory) if f.endswith('.' + extension))
    filelist.sort()
    return(filelist)

class pyinetadb:
	
	def __init__(self,file):	# Initialize an entry using NMR star file
		self.pynmrEntry = pynmrstar.Entry.from_file(file)
	
	def get_entry_info(self):
		entry_list=self.pynmrEntry.get_saveframes_by_category('entry_information')
		self.entryID=self.pynmrEntry.get_tag('_Entry.ID')[0]
		try:
			name=self.pynmrEntry.get_tag('_Entry.BMRB_internal_directory_name')[0] # _Entry.Title _Chem_comp.Name _Entry.BMRB_internal_directory_name
		except:
			name=self.pynmrEntry.get_tag('_Entry.Title')[0].strip()
		name=re.sub('\([()RS0-9+-,]*\)','', name)
		name=re.sub("[()'\[\]]","", name)
		name=re.sub('[ ,-]+','_', name)
		name=re.sub('^_+','', name)
		name=re.sub('^[\+\/-]+','', name)
		self.entryTitle=name
		self.entryVersion=self.pynmrEntry.get_tag('_Entry.NMR_STAR_version')[0].strip()
		sample_results = []
		for sample_loop in self.pynmrEntry.get_loops_by_category("Sample_component"):
			sample_results.append(sample_loop.get_tag(['Type','Mol_common_name','Concentration_val','Concentration_val_units','Concentration_val_err']))
		self.entrySolvent=[]
		for sample_list in sample_results:
			for sample in sample_list:
				if sample[0]=="Solvent" and sample[1] not in self.entrySolvent:
					self.entrySolvent.append(sample[1])

	def get_shifts(self):
		cs_result_sets = []
		self.atomIDMap = []
		for chemical_shift_loop in self.pynmrEntry.get_loops_by_category("Atom_chem_shift"):
			cs_result_sets.append(chemical_shift_loop.get_tag(['Atom_type','Atom_ID','Auth_atom_ID','Val','Val_err','Ambiguity_code']))	
		self.shifts=[]
		for cs_list in cs_result_sets:
			shifts_list={}
			atomIDMap_list = {}
			for cs in cs_list:
				if cs[0] == "C":
					keyC=cs[1]
					if cs[2].startswith("C"):
						mainC=cs[2]
					else:
						mainC=cs[1]						
					if mainC not in shifts_list:
						shifts_list[mainC] = []
					shifts_list[mainC].append(float(cs[3]))
					if keyC not in atomIDMap_list:
						atomIDMap_list[keyC] = mainC
			self.atomIDMap.append(atomIDMap_list)
			self.shifts.append(shifts_list)

	def iter_bonds(self,idx,bd_list):
		bonds_list=[]
		for bd in bd_list:
			if (bd[0] in self.atomIDMap[idx]) and (bd[1] in self.atomIDMap[idx]):
				bonds_list.append([self.atomIDMap[idx][bd[0]],self.atomIDMap[idx][bd[1]]])
		return(bonds_list)

	def get_bonds(self):
		bond_result_sets = []
		for bond_loop in self.pynmrEntry.get_loops_by_category("Chem_comp_bond"):
			bond_result_sets.append(bond_loop.get_tag(['Atom_ID_1','Atom_ID_2']))
		self.bonds=[]
		if len(bond_result_sets)==len(self.shifts):
			for i,bd_list in enumerate(bond_result_sets):
				bonds_list=self.iter_bonds(i,bd_list)
				self.bonds.append(bonds_list)
		elif len(bond_result_sets)==1:
			for i in range(0,len(self.shifts)):
				bonds_list=self.iter_bonds(i,bond_result_sets[0])
				self.bonds.append(bonds_list)
		else:
			print("Bond results don't match chemical shifts!")

	def build_network(self):
		self.networks=[]
		for i,shifts_list in enumerate(self.shifts):
			accessed=[]
			networks_list=[]
			for key in sorted(shifts_list.keys()):
				if len(accessed) < len(self.bonds[i]):
					for idx, pair in enumerate(self.bonds[i]):
						val2=[]
						if idx not in accessed:
							if key in pair:
								val1=[]
								for v1 in shifts_list[pair[0]]:
									for v2 in shifts_list[pair[1]]:
										dq=round(v1+v2,3)
										val1.append([v1,dq])
								val2.append([pair[0],val1])
								val1=[]
								for v1 in shifts_list[pair[1]]:
									for v2 in shifts_list[pair[0]]:
										dq=round(v1+v2,3)
										val1.append([v1,dq])
								val2.append([pair[1],val1])
								accessed.append(idx)
								networks_list.append(val2)
			self.networks.append(networks_list)

	def get_ambig(self):
		self.ambiguity=[]
		for amb in self.shifts:
			ct_c=0
			ct_am=0
			for key,val in amb.items():
				ct_c+=1
				if(len(val)>1):
					ct_am+=1
			self.ambiguity.append(round(float(ct_am/ct_c),3))

def compileEntry(dbObj,num,ct,bugreport):
	# entryID,entryTitle,entryVersion,entrySolvent,shifts,bonds,networks,ambiguity
	uniqID=str(ct)+"::"+dbObj.entryID+"::"+dbObj.entryTitle+"::"+dbObj.entryVersion+"::"+dbObj.entrySolvent[num]
	dbEntry={}
	dbEntry["BMRBName"]=dbObj.entryTitle
	dbEntry["InternalID"]=uniqID
	dbEntry["Version"]=dbObj.entryVersion
	dbEntry["Solvent"]=dbObj.entrySolvent[num]
	dbEntry["ChemicalShifts"]=dbObj.shifts[num]
	dbEntry["Bonds"]=dbObj.bonds[num]
	dbEntry["Networks"]=dbObj.networks[num]
	dbEntry["Ambiguity"]=dbObj.ambiguity[num]
	if (bugreport):
		print(dbEntry["InternalID"])
	
	return(uniqID,dbEntry)

def main(args):	
	files=list_files(args.filedir,"str")
	inetadb={}
	ct=0
	for file in files:
		ct+=1
		if (args.debug):
			print("#########################")
			print("###   Entry details   ###")
			print("#########################")
			print("##Entry number:",ct," || File name:",file)
		fileloc=args.filedir+"/"+file
		try:
			entry=pyinetadb(fileloc)
			entry.get_entry_info()
			if (args.debug):
				print("#Entry info:",entry.entryID,entry.entryTitle,entry.entryVersion)
				print("#Solvents:",len(entry.entrySolvent),"==>",entry.entrySolvent)
			entry.get_shifts()
			if (args.debug):
				print("#Shifts:",len(entry.shifts),"==>",entry.shifts)
			entry.get_bonds()
			if (args.debug):
				print("#Bonds:",len(entry.bonds),"==>",entry.bonds)
			entry.build_network()
			if (args.debug):
				print("#Networks:",len(entry.networks),"==>",entry.networks)
				for i,val in enumerate(entry.networks):
					print("\t",i,"==========")
					print("\t",entry.shifts[i])
					print("\t",entry.bonds[i])
					print("\t",entry.networks[i])
			entry.get_ambig()
			if (args.debug):
				print("#Ambiguity:",len(entry.ambiguity),"==>",entry.ambiguity)
				print("############   Compiling entry information   ############")
			for i in range(0,len(entry.entrySolvent)):
				if (i>0 and i<len(entry.entrySolvent)):
					ct+=1 
				(key,dbitems)=compileEntry(entry,i,ct,args.debug)
				inetadb[key]=dbitems

		except:
			print("Skipped %s !!!" % (file))

	with open(args.outfile, 'w') as outf1:
		json.dump(inetadb, outf1)

if __name__ == '__main__':

	parser = ArgumentParser(description='Script to generate the INETA database from the BMRB star files.')
	parser.add_argument('-f', dest='filedir', required=True,
		 help='Folder with all the NMR STAR (*.str) files downloaded from BMRB')
	parser.add_argument('-o', dest='outfile', required=True,
		help='Output filename for the json formatted INETA DB.')
	parser.add_argument('-d', dest='debug', action="store_true", default=False, required=False,
		help='Use this option to run the script in debug mode. Generates output for every step to track errors.')
	args = parser.parse_args()
	main(args)