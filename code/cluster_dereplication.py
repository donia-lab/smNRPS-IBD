#!/bin/env python

"""
Created on Tues Aug 14 8:59:11 2018
@author: francinecamacho
"""

"""This script takes in a BLAST tabular output file: a fasta file used to BLAST against itself to parse out the
 duplicated clusters, proteins or genes from the fasta file. Assumption is that the BLAST tabular file is run with 
 --perc_identity parameter in BLAST. This is necessary for BLAST to calculate the query coverage on hits that passed
 with the desired percent identity. The output is the de-replicated fasta file, and a text file with the unique cluster, 
 proteins or gene names."""
 
from Bio import SeqIO
import pandas as pd 
import networkx as nx 
import os 

#Function to make a panda data frame to iterate our BLAST tabular file output. 
def makeDataFrame(PATH, tabularHeader, coverageCutOff, perc_identity):
	names = tabularHeader.split(" ")
	dataframe = pd.read_csv(PATH, sep="\t", names=names, header=None)
	# filter dataframe to just have hits that meet qcovs threshold and pident threshold and qseqid are not identical to sseqid
	filteredDF = dataframe[(dataframe.sseqid != dataframe.qseqid) & (dataframe.qcovs >= coverageCutOff) &
						   (dataframe.pident >= perc_identity)]

	uniqueQID = dataframe['qseqid'].unique() # retrieve list of unique queries
	print("Number of inputted BGCs:", len(uniqueQID))
	return uniqueQID, filteredDF

def makeLengthDict(df):
	resultDict = {}
	for index, row in df.iterrows():
		queryName = df.at[index, 'qseqid']
		subjectName = df.at[index, 'sseqid']
		qlength = df.at[index, 'qlen']
		subjectLen = df.at[index, 'slen']
		
		if (queryName not in resultDict):
			resultDict[queryName] = qlength
		if (subjectName not in resultDict):
			resultDict[subjectName] = subjectLen
	return(resultDict)

def createNetworkGraph(df):
	graphtype = nx.Graph()
	bgc_network = nx.from_pandas_edgelist(df, source="qseqid", target="sseqid", create_using=graphtype)
	subgraphs = list(nx.connected_component_subgraphs(bgc_network, copy=False))

	print("Number of hubs:", len(subgraphs))
	return (subgraphs)

# Function to find the longest match for a query in terms of sequences length
def findMaxBGC(matchesArray, qlenDict):
	maxBGC = None
	for match in matchesArray:
		if (maxBGC == None):
			maxBGC = match
		if (qlenDict[maxBGC] < qlenDict[match]):
			maxBGC = match
	return(maxBGC)

def findUniqueBGCs(subgraphs, bgcList, outdir, outfile, allBGCs):

	os.chdir(outdir)
	uniqueBGCs = []
	all_matches_BGCs = []
	for subgraph in subgraphs:
		print("======================================================")
		print ('networkx subgraph:', subgraph.nodes())
		all_matches_BGCs = list(subgraph.nodes()) + all_matches_BGCs
		longest_bgc = findMaxBGC(subgraph.nodes(), bgcList)
		uniqueBGCs.append(longest_bgc)
		print ('Longest BGC in matches:', longest_bgc)
		print("======================================================")

	bgcs_no_matches = set(allBGCs) - set(all_matches_BGCs) # bgcs that did not find a match based on conditions 

	bgcs_no_matches_list = list(bgcs_no_matches) # convert to list 

	combined_unique_bgcs = bgcs_no_matches_list + uniqueBGCs

	if len(combined_unique_bgcs) == 0:
		print("Error: Could not identify any duplicates from input file")
	else:
		uniqueKeys_series = pd.Series(combined_unique_bgcs)
		uniqueKeys_DF = uniqueKeys_series.to_frame("bgcName")
		bgcNameFile = outfile + "_uniqueBGCNames.txt"
		uniqueKeys_DF.to_csv(bgcNameFile, index=False,
							 sep='\t')  # write dataframe to csv format (text file) of unique BGCs name
		print ('Total unique BGCs:', len(combined_unique_bgcs))

	# return unique cluster list to map list to our BGC master list to make a fasta file
	return(combined_unique_bgcs)
		
"""Function to create a fasta file using the list of unique BGCs and map those to the original master bgc fasta file """
def createFastaFile(uniqueBGCList, bgcMasterList, outdir, outfile):

	os.chdir(outdir)
	fastafileName = outfile + "_uniqueBGCs.fa"
	with open(fastafileName, 'w') as uniqueFile:
		for seq_record in SeqIO.parse(bgcMasterList, 'fasta'):
			if seq_record.id in uniqueBGCList:
				seq_record.description = ""
				SeqIO.write(seq_record, uniqueFile, "fasta")
			else:
				continue
	uniqueFile.close()

def main(tabular_file, outdir, outfile, perc_identity, coverage_cutoff, bgc_master_file,
		 tabular_file_header):

	statinfo = os.stat(tabular_file)
	if statinfo.st_size != 0:  #check size of tabular file 

		unique_qid_array, dfObject = makeDataFrame(tabular_file, tabular_file_header, coverage_cutoff, perc_identity)

		query_len_dict = makeLengthDict(dfObject)
		query_network_graph = createNetworkGraph(dfObject)
		unique_bgc_list = findUniqueBGCs(query_network_graph, query_len_dict, outdir, outfile, unique_qid_array)
		createFastaFile(unique_bgc_list, bgc_master_file, outdir, outfile)

	else:
		print("ERROR: Inputted tabular file is empty.")

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('--tabular_file', type=str, required=True, help='blast tabular file results')
	parser.add_argument('--outdir', type=str,
						help='directory to output files')
	parser.add_argument('--outfile', type=str, required=True, help='name of cohort or project')
	parser.add_argument('--perc_identity', type=int, required=False, default=95, help='default is 95')
	parser.add_argument('--coverage_cutoff', type=int, required=False, default=95, help='default is 95')
	parser.add_argument('--bgc_master_file', type=str, required=True,
						help='fasta file used for BLAST and de-replication')
	parser.add_argument('--tabular_file_header', type=str, required=False,
						default="sseqid qseqid slen qlen qcovs pident Evalue qstart qend", help='sseqid qseqid slen qlen qcovs pident Evalue qstart qend')

	args = parser.parse_args()

	main(args.tabular_file, args.outdir, args.outfile, args.perc_identity, args.coverage_cutoff,
		 args.bgc_master_file, args.tabular_file_header)		