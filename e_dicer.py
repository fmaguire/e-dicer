#!/usr/bin/env python

from Bio import SeqIO
import os.path

def parse_fasta(fasta_file_name):
    '''
    Function to a parse a fasta file into a list of Seq object using SeqIO
    input: fasta_file_name
    output: iterator for SeqRecord objects
    '''

    if not os.path.isfile(fasta_file_name):
        raise IOError

    seq_generator = SeqIO.parse(fasta_file_name, 'fasta')

    return seq_generator
