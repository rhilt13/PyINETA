"""Functions for matching INETA networks to the database.

This includes all the functions used by PyINETA to match the identified networks in the query spectrum to the INADEQUATE spectra of known metabolites in the INETA database.
Includes the following functions:
	* prepUnknowns - Add tags and connection names to the unknown networks.
	* matchTag - Match the tags from matched peaks to the list of unknown connections.
	* matchDatabase - Match unknown networks to networks in the InETA database.
"""

import itertools

def prepUnknowns (Pvals,Pairs):
	"""Add tags and connection names to the unknown networks.

	This function adds CX tags to the peaks that wil be used when matching them to teh database.

	Args:
		Pvals (list): List of all points in a single network.
		Pairs (list): A list of all horizontally connected pairs of points included in a network.
	
	Returns:
		list : List of Tags for naming the unknown peaks in a network.
		list : List of lists with pairs of connected tags.
		list : List of 13C and DQ ppm values corresponding to the list of connected tags.
	"""

	TagDict={}
	Ptags=[]
	for i in range(0,len(Pvals)):
		j=str('CX'+str(i+1))
		TagDict[Pvals[i]]=j
		Ptags.append('CX'+str(i+1))
	newP=[]
	for i in Pvals:
		new1=[]
		for j in Pairs:
			j=tuple(map(tuple, j))
			if i in j:
				new1.append(tuple(j))
		newP.extend(new1)
	newP.sort()
	FinalPairs=list(newP for newP,_ in itertools.groupby(newP))
	unknownConn= [[TagDict[j] for j in x] for x in FinalPairs]
	return(Ptags,unknownConn,FinalPairs)

def matchTag (match_peaks,unknown_conn):
	"""Match the tags from matched peaks to the list of unknown connections.

	This function is used to map known peak assignments to the matched peaks in an unknown network.

	Args:
		match_peaks (dict): Dict mapping the unknown peaks to matched peaks from database.
		unknown_conn (list): List of all unknown connections of a network.

	Returns:
		dict : Dict mapping pairs of peaks with unknown connections to known connections from database networks.
		dict : Dict mapping unmatched pairs of conenctions to '?'.
		int : Count of hte number of matched connections.
	"""

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
	return (final_match,no_match,hits)
	
def matchDatabase (json_db,Pval,tag,unknown_conn,amb_tol,near_tol,match_tol,top_tol,hit_tol,cov_tol,NetNum,Pts):
	"""Match unknown networks to networks in the InETA database.

	Args:
		json_db (dict): The INETA database in json format.
		Pval (list): List of all points in a single network.
		tag (list): List of tags for the unknown peaks.
		unknown_conn (list): List of unknown connections.
		amb_tol (float): Amibiguity tolerance (Removes all database entries with ambiguity higher than this tolerance).
		near_tol (float): Tolerance for distance in the 13C dimension between unknown peak and a match database peak.
		match_tol (int): Tolerance for number of matches in a single network.
		top_tol (float): Topology tolerance (How far the point can be in all directions from the match peak point).
		hit_tol (float): Hit score threshold.
		cov_tol (float): Coverage score threshold.
		NetNum (int): Network number.
		Pts (list): List of pairs of horizontally connected peaks.

	Returns:
		list : A list of networks with their corresponding hits and hit details.
	"""

	X = [float(i[0]) for i in Pval]
	Y = [float(i[1]) for i in Pval]
	hitList=[]
	hitList.append("Network"+ str(NetNum))
	hitList.append(X)
	hitList.append(Y)
	hitList.append(tag)
	hitList.append(Pts)
	for i in json_db:
		matchCnt=0
		matchCS=[]
		if (len(json_db[i]['Networks'] ) >=1):
			if (json_db[i]['Ambiguity'] <= float(amb_tol)):
				ambScore=json_db[i]['Ambiguity']
				for xval in X:
					for CS in json_db[i]['ChemicalShifts']:
						if CS not in matchCS:
							CSlist=json_db[i]['ChemicalShifts'][CS]
							meanCS=sum(CSlist) / float(len(CSlist))
							if (abs(xval-meanCS)<= float(near_tol)):
								matchCS.append(CS)
								matchCnt+=1
				if (matchCnt>= float(match_tol)):
					found={}
					matchPeaks={}
					for j in json_db[i]['Networks']:
						for k in j:
							break_outer_loop=False
							for l in k[1]:
								for ind,xval in enumerate(X):
									if (((xval-l[0])**2 + (Y[ind]-l[1])**2)**(1/2.0) <= float(top_tol)):
										matchPeaks[tag[ind]]=k[0]
										found[tag[ind]]=1
										break_outer_loop=True
										break
								if break_outer_loop: break
					(fin_match,no_match,hitCount)=matchTag(matchPeaks,unknown_conn)
					hitScore=float("{0:.3f}".format(float(hitCount)/float(len(json_db[i]['Networks']))))
					CovScore=float("{0:.3f}".format(float(len(found))/float(len(tag))))
					if (hitScore >= float(hit_tol) and CovScore >= float(cov_tol)):
						hitOut=[]
						hitOut.append(i)
						hitOut.append([json_db[i]['Networks']])
						hitOut.append(fin_match)
						hitOut.append(no_match)
						hitOut.append(ambScore)
						hitOut.append(hitScore)
						hitOut.append(CovScore)
						hitList.append(hitOut)
	return(hitList)