#!/bin/env python

"""
Created on Thurs Mar 3 20:01:31 2016
@author: francinecamacho
"""

from Bio import SeqIO
import pandas as pd
import os 

"""This script will take the tabular file as input to detect BGCs based on percent identity (95%) 
and query coverage (95%) criteria to find the BGC taxa producer based on ref_seq NCBI database. The output
is a tabular result file with BGC name, percent identity, and coverage and taxa name.
"""

#Function to make a panda data frame tabular file 
def makeDataFrame (PATH, perc_ident_cutoff, coverage_cutoff): 
	names = ['seqid', 'stitle', 'sacc', 'qseqid', 'qlen', 'qcovs','pident', 'Evalue', 'qstart','qend']
	dataframe = pd.read_csv(PATH, sep = "\t", names=names, header=None)
	# filter dataframe by qcovs and percent
	filter_dataframe = dataframe[(dataframe.qcovs >= coverage_cutoff) & (dataframe.pident >= perc_ident_cutoff)]
	# unique_bgcs= dataframe['qseqid'].unique() 
	return filter_dataframe


def initializeDict(bgc_list):
	inital_dict = {}
	for i in range(0, len(bgc_list)):
		inital_dict[bgc_list[i]] = 'N/A'
	return inital_dict 

def mapTaxa(df, species_dict):
	for index, row in df.iterrows():
		bgc_name = df.at[index,'qseqid']
		species = df.at[index, 'stitle']
		coverage = df.at[index, 'qcovs']
		accession_id = df.at[index, 'sacc']
		percent_identity = df.at[index, 'pident']
		species_title = species.replace(" ", "_")
		species_dict[bgc_name] = [species_title, accession_id, percent_identity,coverage]
	for bgc in species_dict: 
		if species_dict[bgc] == 'N/A': 
			species_dict[bgc] = ['N/A', 'N/A','0.0', '0.0']
	return(species_dict)

# def updateFASTA(species_dict, fasta_file, outdir, outfile):
# 	os.chdir(outdir)
# 	bgc_fasta_file = SeqIO.parse(open(fasta_file),'fasta')
# 	with open(outfile, 'w') as updated_fasta:
# 		for seq_record in bgc_fasta_file:
# 			seq_id = seq_record.id
# 			if species_dict[seq_id][0]!= 'N/A':
# 				seq_record.id = seq_id + "__" + species_dict[seq_id][0]
# 				seq_record.description = "" 
# 				SeqIO.write(seq_record, updated_fasta, "fasta")
# 			else: 
# 				seq_record.description = "" 
# 				SeqIO.write(seq_record, updated_fasta, "fasta")

def main(tabular_file, bgc_fasta_file, outdir, outfile, perc_ident_cutoff, coverage_cutoff):

	output_df = pd.DataFrame() # initialize Data Frame 

	blast_df = makeDataFrame(tabular_file, perc_ident_cutoff, coverage_cutoff)
	bgc_fasta_file_seqIO = SeqIO.parse(bgc_fasta_file, "fasta")
	bgc_list = [record.id for record in bgc_fasta_file_seqIO ]

	bgc_dict = initializeDict(bgc_list)
	print(bgc_dict)
	species_result_dict = mapTaxa(blast_df, bgc_dict)  
	df2 = [(k, v[0], v[1], v[2], v[3]) for k, v in list(species_result_dict.items())]
	output_df = output_df.append(df2) # append list to our initalized dataframe 


	os.chdir(outdir)
	output_df.columns = ["BGC_NAME", "TAXA_NAME", "ACC_ID", "PERC_IDENT", "COVERAGE"] # rename column names 
	output_df.to_csv(outfile, index=False, sep='\t') # write dataframe to csv format (text file)

	# updateFASTA(species_result_dict, bgc_fasta_file, outdir, outfile)


if __name__ == '__main__':

	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--tabular_file', type =str, required=True)
	parser.add_argument('--bgc_fasta_file', type =str, required=True)
	parser.add_argument('--outdir', type= str) 
	parser.add_argument('--outfile', type= str, required=True)
	parser.add_argument('--perc_ident_cutoff', type= int, required=False, default=95)
	parser.add_argument('--coverage_cutoff', type= int, required=False, default=95)



	args = parser.parse_args()


	main(args.tabular_file, args.bgc_fasta_file, args.outdir, args.outfile, args.perc_ident_cutoff, args.coverage_cutoff)