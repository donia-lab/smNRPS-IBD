#!/bin/env python

"""
Created on Tuesday Jan 26 14:32:16 2016
@author: francinecamacho
"""

from Bio import SeqIO

def scaffoldCounter(assemblyFile, scaffold_len_cutoff):

        for ct, scaffold_record in enumerate(SeqIO.parse(assemblyFile, "fasta")):
                if int(len(scaffold_record.seq)) < scaffold_len_cutoff:
                        break
        print (ct) 


def main(assemblyFile, scaffold_len_cutoff):

        limitNumber = scaffoldCounter(assemblyFile,scaffold_len_cutoff)


if __name__ == '__main__':

        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('--assemblyFile', type=str)
        parser.add_argument('--scaffold_len_cutoff', type=int, default=5000)
        args = parser.parse_args()

        main(args.assemblyFile, args.scaffold_len_cutoff)
