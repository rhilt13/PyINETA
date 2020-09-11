#! /usr/bin/env python

import os
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
		self.entryTitle=self.pynmrEntry.get_tag('_Entry.Title')[0]
		self.entryVersion=self.pynmrEntry.get_tag('_Entry.NMR_STAR_version')[0].strip()
		sample_results = []
		for sample_loop in self.pynmrEntry.get_loops_by_category("Sample_component"):
			sample_results.append(sample_loop.get_tag(['Type','Mol_common_name','Concentration_val','Concentration_val_units','Concentration_val_err']))
		print(sample_results)
		self.entrySolvent=[]
		for sample_list in sample_results:
			for sample in sample_list:
				if sample[0]=="Solvent":
					print(sample[1])
					self.entrySolvent.append(sample[1])

	def get_shifts(self):
		cs_result_sets = []
		self.atomIDMap = []
		for chemical_shift_loop in self.pynmrEntry.get_loops_by_category("Atom_chem_shift"):
			cs_result_sets.append(chemical_shift_loop.get_tag(['Atom_type','Atom_ID','Auth_atom_ID','Val','Val_err','Ambiguity_code']))	
		self.shifts=[]
		print(cs_result_sets)
		for cs_list in cs_result_sets:
			shifts_list={}
			atomIDMap_list = {}
			for cs in cs_list:
				if cs[0] == "C":
					if cs[2] not in shifts_list:
						shifts_list[cs[2]] = []
					shifts_list[cs[2]].append(float(cs[3]))
					if cs[1] not in atomIDMap_list:
						atomIDMap_list[cs[1]] = cs[2]
			self.atomIDMap.append(atomIDMap_list)
			self.shifts.append(shifts_list)


	def get_bonds(self):
		bond_result_sets = []
		for bond_loop in self.pynmrEntry.get_loops_by_category("Chem_comp_bond"):
			bond_result_sets.append(bond_loop.get_tag(['Atom_ID_1','Atom_ID_2']))
		self.bonds=[]
		for idx,bd_list in enumerate(bond_result_sets):
			bonds_list=[]
			for bd in bd_list:
				if (bd[0] in self.atomIDMap[idx]) and (bd[1] in self.atomIDMap[idx]):
					bonds_list.append([self.atomIDMap[idx][bd[0]],self.atomIDMap[idx][bd[1]]])
			self.bonds.append(bonds_list)

	def build_network(self):
		accessed=[]
		self.networks=[]
		for shifts_list in self.shifts:
			for key in shifts_list:
				if len(accessed) < len(self.bonds):
					for idx, pair in enumerate(self.bonds):
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

	def get_ambig(self):
		ct_c=0
		ct_am=0
		for key,val in self.shifts.items():
			ct_c+=1
			if(len(val)>1):
				ct_am+=1
		self.ambiguity=round(float(ct_am/ct_c),3)

def compileEntry(dbObj):
	# entryID,entryTitle,entryVersion,entrySolvent,shifts,bonds,nets,ambig):
	uniqID=str(randint(1,100))+"::"+dbObj.entryID+"::"+dbObj.entryTitle+"::"+dbObj.entryVersion+"::"+dbObj.entrySolvent
	dbEntry={}
	dbEntry["BMRBName"]=dbObj.entryTitle
	dbEntry["InternalID"]=uniqID
	dbEntry["Version"]=dbObj.entryVersion
	dbEntry["Solvent"]=dbObj.entrySolvent
	dbEntry["ChemicalShifts"]=dbObj.shifts
	dbEntry["Bonds"]=dbObj.bonds
	dbEntry["Networks"]=dbObj.networks
	dbEntry["Ambiguity"]=dbObj.ambiguity
	
	return(uniqID,dbEntry)

def main(args):	
	files=list_files(args.filedir,"str")
	inetadb={}
	for file in files:
		print(file)
		fileloc=args.filedir+"/"+file
		try:
			entry=pyinetadb(fileloc)
			print(entry)
			entry.get_entry_info()
			if (args.debug):
				print(entry.entryID,entry.entryTitle,entry.entryVersion,entry.entrySolvent)
			entry.get_shifts()
			if (args.debug):
				print(entry.shifts)
			entry.get_bonds()
			if (args.debug):
				print(entry.bonds)
			entry.build_network()
			if (args.debug):
				print(entry.networks)
			entry.get_ambig()
			if (args.debug):
				print(entry.ambiguity)
			(key,dbitems)=compileEntry(entry)
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