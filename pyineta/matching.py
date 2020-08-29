import itertools

def prepUnknowns (Pvals,Pairs):
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
	# print FinalPairs
	#print TagDict
	unknownConn= [[TagDict[j] for j in x] for x in FinalPairs]
	return(Ptags,unknownConn,FinalPairs)

def matchTag (match_peaks,candidate_peaks,unknown_conn):
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
#MatchTag
	X = [float(i[0]) for i in Pval]
	Y = [float(i[1]) for i in Pval]
	hitList=[]
	hitList.append("Network"+ str(NetNum))
	hitList.append(X)
	hitList.append(Y)
	hitList.append(tag)
	hitList.append(Pts)
	# print X
	# print Y
	# print Pval
	for i in json_db:
		# print i
		matchCnt=0
		matchCS=[]
		if (len(json_db[i]['Networks'] ) >=1):
			# print json_db[i]['Networks'], len(json_db[i]['Networks'] )
			if (json_db[i]['Ambiguity'] <= float(amb_tol)):
				ambScore=json_db[i]['Ambiguity']
				for xval in X:
					for CS in json_db[i]['ChemicalShifts']:
						if CS not in matchCS:
							CSlist=json_db[i]['ChemicalShifts'][CS]
							meanCS=sum(CSlist) / float(len(CSlist))
							# print i, xval, CS,"////", meanCS, matchCS
							if (abs(xval-meanCS)<= float(near_tol)):
								# print i, xval, meanCS, "==", abs(xval-meanCS), "//", near_tol
								matchCS.append(CS)
								matchCnt+=1
				if (matchCnt>= float(match_tol)):
					# print i, matchCS
					found={}
					candPeaks={}
					matchPeaks={}
					for j in json_db[i]['Networks']:
						# print "j==>",j
						candPeaks.setdefault(j[0][0], set()).add(j[1][0])
						for k in j:
							# print "k==>",k
							break_outer_loop=False
							for l in k[1]:
								# print "l==>",l
								for ind,xval in enumerate(X):
									mp=[]
									# print ind,xval,Y[ind],tag[ind]
									# print i,"=",k[0],l,"//",tag[ind],xval,Y[ind],"//",(xval-l[0])**2,(Y[ind]-l[1])**2,((xval-l[0])**2 + (Y[ind]-l[1])**2)**(1/2.0),top_tol
									if (((xval-l[0])**2 + (Y[ind]-l[1])**2)**(1/2.0) <= float(top_tol)):
										matchPeaks[tag[ind]]=k[0]
										# print "matchPeaks",len(tag),tag[ind],matchPeaks
										# print i,"HEE",l,xval,Y[ind],(xval-l[0])**2,(Y[ind]-l[1])**2,((xval-l[0])**2 + (Y[ind]-l[1])**2)**1/2,top_tol
										found[tag[ind]]=1
										break_outer_loop=True
										break
								if break_outer_loop: break
					# print(i)
					# print("matchPeaks==>",matchPeaks)
					# print("candPeaks==>",candPeaks)
					# print("unknown_conn==>",unknown_conn)
					(fin_match,no_match,hitCount)=matchTag(matchPeaks,candPeaks,unknown_conn)
					hitScore=float("{0:.3f}".format(float(hitCount)/float(len(json_db[i]['Networks']))))
					CovScore=float("{0:.3f}".format(float(len(found))/float(len(tag))))
					# print len(found),found,len(tag),tag 
					# print "fin_match==>",fin_match,
					# print "no_match==>",no_match,
					# print "hitCount=",hitCount," Total=",len(json_db[i]['Networks']),"=> hitScore=",hitScore
					if (hitScore >= float(hit_tol) and CovScore >= float(cov_tol)):
						hitOut=[]
						hitOut.append(i)
#						print json_db[i]['Networks']
						hitOut.append([json_db[i]['Networks']])
						# IndDict=ParseDbNetwork(i, json_db[i]['Networks'])
						# hitOut.append(IndDict)
						hitOut.append(fin_match)
						hitOut.append(no_match)
						hitOut.append(ambScore)
						hitOut.append(hitScore)
						hitOut.append(CovScore)
						hitList.append(hitOut)
						# print hitList
	return(hitList)