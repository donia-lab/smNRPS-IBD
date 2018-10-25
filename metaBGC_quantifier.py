#!/bin/env python

"""
Created on Fri Feb 2 14:52:23 2018
@author: francinecamacho
"""

import collections
import pandas as pd
import os 
import logging


class BGC:
	def __init__(self, name, length, hitCounter=0):
		self.name = name
		self.length = length
		self.hitCounter = hitCounter

"""Function populates a postion counter collection for the length of a sequences reads 
that has mapped to a given BGC."""
def populateCounter (startnum, stopnum,qseqID, counterdic): 
	for i in range (startnum, stopnum+1): 
		counterdic[qseqID][i]+=1
	return counterdic # defaultdict(<class 'collections.Counter'>, {'BGC1': Counter({1920: 1, 1921: 1, 1922: 1, 1923: 1, 1924: 1, etc...})

"""Function creates a pandas dataframe from BLAST tabular output"""
def makeDataFrame (PATH, colNames): 

	dataframe = pd.read_csv(PATH, names=colNames, header=None, delim_whitespace=True)
	
	logger.debug('dataframe:\n%s', dataframe)

	uniqueQID= dataframe['qseqid'].unique() 
	return uniqueQID,dataframe

"""Function initializes a collection counter based on the total unique BGC."""
def initializeCounter(qIDArray):
	counter = collections.defaultdict(collections.Counter)
	for j in range(0, len(qIDArray)):
		qID = qIDArray[j]
		counter[qID]
	return counter # defaultdict(<class 'collections.Counter'>, {'BGC1': Counter()})

"""Function filters sequences reads based on read coverage by position of read mapped to a BGC."""
def filterReadByPos(smin, smax, qmin, qmax, slength, qlength, percentReadCoveredByCluster):
	internal = None # if read is in the middle  
	filterPass = False

	if (qmin - 1 >= smin - 1) and (qlength - qmax >= slength - smax): # check for middle 
		internal = True
	else:
		internal = False

	if (internal == True):  # read pos == middle
		if (percentReadCoveredByCluster >= 90):
			filterPass = True
	elif (internal == False):  # read pos == edge

		if (percentReadCoveredByCluster >= 50):
			filterPass = True

	return filterPass

"""Function to populate each counter by iterating through dataframe rows
and calculate coverage and abundance"""
def calculateBreathAndDepth(counterObject, df, sampleName, sequenceReads, perc_identity):

	queryLenDict ={}
	completeCounter = {}

	for index, row in df.iterrows():
		logger.debug('index: %s', index)
		logger.debug('row:\n%s', row)
		bgcName = df.at[index,'qseqid']
		smin = min(df.at[index, 'sstart'], df.at[index, 'send'])
		smax = max(df.at[index, 'sstart'], df.at[index, 'send'])
		qmin = min(df.at[index, 'qstart'], df.at[index, 'qend'])
		qmax = max(df.at[index, 'qstart'], df.at[index, 'qend'])
		slength = df.at[index, 'slen']
		qlength = df.at[index, 'qlen']

		if bgcName not in queryLenDict:
			bgcClass = BGC(bgcName, qlength)
			queryLenDict[bgcName] = bgcClass 
			
		#Filter reads by percent identity cutoff (95%) and edge/internal read coverage cutoff, 50% and 90% respectively. 
		if (df.at[index, 'pident'] >= perc_identity):
			
			percentReadCoveredByCluster = (((smax - smin) +1) / slength) * 100
			readPass = filterReadByPos(smin, smax, qmin, qmax, slength, qlength, percentReadCoveredByCluster)

			if (readPass == True): 
				completeCounter = populateCounter(qmin, qmax, bgcName, counterObject)
				queryLenDict[bgcName].hitCounter+=1 
			else:
				continue 
		else:
			continue 
	
	resultsDict = {} # dictionary to display our results {qID: [COVERAGE, ABUNDANCE], qID: [COVERAGE, ABUNDANCE], etc..}

	for nID in completeCounter: # for each of the nID (keys) in our completeCounter find the total # of positions not zero
		logger.debug('nID: %s', nID)
		
		valuesNotZero = 0

		newDict = completeCounter[nID] # {qID:queryLength, qID: queryLength, etc..}

		try : 

			for key in newDict:
				value = newDict[key] # counts of each positions (used to calculate relativeAbundance depth)

				if value > 0: # if values are greater than zero add one to our counter of valuesNotZero and add all the counts (sumOfValues)
					valuesNotZero+=1

			seqReadsPerMillion = sequenceReads/1000000
			totalHitsPerBGC = queryLenDict[nID].hitCounter
			bgcPerKiloBase = queryLenDict[nID].length/1000 

			# relativeAbundance = totalHitsPerBGC *(seqReadsPerMillion)*(bgcPerKiloBase) # RPKM
			relativeAbundance = totalHitsPerBGC / (bgcPerKiloBase * seqReadsPerMillion)
		except ZeroDivisionError as relativeAbundance: # if counter object is  zero 
			relativeAbundance = 0 

		try: 
			# calculate percent coverage by total positions not zero divided by length of query
			coverage = (float(valuesNotZero) / queryLenDict[nID].length) * 100.0

		except ZeroDivisionError as coverage: # if counter object is zero 
			coverage =  0.0
		
		resultsDict[nID] = [coverage, relativeAbundance,sampleName] #{qID: [COVERAGE, ABUNDANCE], qID: [COVERAGE, ABUNDANCE], etc..}

	return resultsDict # {qID: [COVERAGE, ABUNDANCE, SampleID], ...}


def main(sampleID, sample_reads_total, outfile, tabularFilePath, perc_identity,tabular_colnames):

	outputDF = pd.DataFrame()

	df_colnames = tabular_colnames.split()

	with open(sample_reads_total) as f:
		totalReads = int(f.read().rstrip('\n'))
		logger.debug('totalReads: %d', totalReads)
	
	statinfo = os.stat(tabularFilePath)
	if statinfo.st_size != 0: #if tabular file is not an empty continue else skip 
			
		uniqueQIDArray, dfObject = makeDataFrame(tabularFilePath,df_colnames)
		counterStructure = initializeCounter(uniqueQIDArray)
		resultsDF = calculateBreathAndDepth(counterStructure, dfObject, sampleID, totalReads, perc_identity)

		logger.debug('resultsDF:\n%s', resultsDF)
		df2 = [(k, v[0], v[1], v[2]) for k, v in resultsDF.items()] # convert dictionary to list   
		if df2: # if there are no blast hit that made it through the cutoffs 
			logger.debug('df2:\n%s', df2)
			outputDF = outputDF.append(df2) # append list to our initalized dataframe 
		else: logger.info('No reported blast hits using cutoff: %s', sampleID)

	else:
		logger.info('Tabular File is empty for: %s', sampleID)

	
	#outputDF.columns = ["QueryID", "Breadth", "relativeAbundance", "Sample"] # rename column names
	outputDF = outputDF[(outputDF != 0).all(1)] # filter rows with zero (1) == any columns 
	outputDF.to_csv(outfile, index=False, sep='\t', header=False) # remove header names 


if __name__ == '__main__':
	
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--sampleID', required=True, type=str) #sample directory folder name 
	parser.add_argument('--sample_reads_total', required=True) # path to the read file name 
	parser.add_argument('--outfile', required=True, type= str) 
	parser.add_argument('--perc_identity', required=False, type=int, default=95)
	parser.add_argument('--tabular_file', required=True)
	parser.add_argument('--tabular_colnames', nargs='+', required=False, default = "sseqid slen sstart send qseqid qlen qstart qend pident Evalue") 
	parser.add_argument('-v', '--verbose', action='store_true')

	args = parser.parse_args()

	logger = logging.getLogger('metaBGC_quantifier')
	if args.verbose:
		logLevel = getattr(logging, 'DEBUG')
	else:
		logLevel = getattr(logging, 'INFO')

	logger.setLevel(logLevel)
	ch = logging.StreamHandler()
	ch.setLevel(logLevel)

	formatter = logging.Formatter('%(asctime)s - %(funcName)s - %(lineno)d - %(levelname)s - %(message)s', datefmt='%Y-%m-%d-%H%M%S')
	ch.setFormatter(formatter)
	logger.addHandler(ch)


	logger.info('sampleID: %s', args.sampleID)
	logger.info('sample_reads_total: %s', args.sample_reads_total)
	logger.info('outfile: %s', args.outfile)
	logger.info('tabular file: %s', args.tabular_file)
	logger.info('perc_identity: %d' % args.perc_identity)
	logger.info('tabular_colnames: %s', args.tabular_colnames)
	
	main(args.sampleID, args.sample_reads_total, args.outfile, args.tabular_file, args.perc_identity, args.tabular_colnames)

