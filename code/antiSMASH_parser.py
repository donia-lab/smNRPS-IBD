#!/bin/env python
"""
Created on Wed Jun  13 12:56:23 2016
@author: francinecamacho
"""
from Bio import SeqIO
import os 

"""Script parses antiSMASH genbank file to create a master fasta file for detected BGCs in scaffold >5KB in sample."""
def parseAntismashGBK(gbk_ext, sample_id, gbk_path, outfile, outdir): 

	os.chdir(outdir)

	with open(outfile, 'w') as parsed_file:

		for file_name in os.listdir(gbk_path):
			
			os.chdir(gbk_path)
			if file_name.endswith(gbk_ext): 

				for seq_record in SeqIO.parse(file_name, 'genbank'):
					scaffold_id = seq_record.description
					features = seq_record.features
					# scaffold_len = len(str(seq_record.seq))

					for gb_feature in features:  # iterate through features 

						# feature is the start of a cluster 
						if gb_feature.type == 'cluster': 
							cluster_start = gb_feature.location.nofuzzy_start
							cluster_end = gb_feature.location.nofuzzy_end
							cluster_len =  len(gb_feature.location)
							bgc_type = ''.join(gb_feature.qualifiers.get('product'))

							if bgc_type == '': # if no bgc type information 
								bgc_type = "unknown"

							cluster_seq = str(gb_feature.extract(seq_record.seq)) 
							fasta_info = ">" + sample_id + '__' + scaffold_id + '__' + str(cluster_len) + '__' + bgc_type + '__' + 'ANTISMASH'+ '__' + str(cluster_start) + '_' + str(cluster_end) + '\n' + cluster_seq +'\n'

							os.chdir(outdir)
							parsed_file.write(fasta_info)

						else: # if no cluster in feature is detected
							continue
			else: # if no antismash output file is detected 
				continue 

def main(sample_id, gbk_path, outdir, outfile):

	gbk_ext = ".final.gbk"
	
	parseAntismashGBK(gbk_ext, sample_id, gbk_path, outfile, outdir)

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--sample_id', type=str, required=True, help='sample name/id')
	parser.add_argument('--gbk_path', type=str, required=True, help="directory where sample's antiSMASH genbank is located")
	parser.add_argument('--outdir', type= str, required=True, help='directory where to ouput parsed fasta file') 
	parser.add_argument('--outfile', type= str, required=True, help='name of fasta file with parsed BGCs') 


	args = parser.parse_args()


	main(args.sample_id, args.gbk_path, args.outdir, args.outfile)