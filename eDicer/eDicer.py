#!/usr/bin/env python

import Bio
from Bio import SeqIO
import os.path
import subprocess
import warnings
import argparse


def get_parser():

    parser = argparse.ArgumentParser(description='Generate dicer fragments'
                                                 'from a multi-line fasta')
    parser.add_argument('--fasta', '-f',
                        action='store',
                        type=str,
                        dest='input_file',
                        required=True,
                        help='Input fasta file for dicer fragment generation '
                             '(REQUIRED)')

    parser.add_argument('--output', '-o',
                        action='store',
                        type=str,
                        dest='output_file',
                        default='output.fas',
                        help='Filename for fragment output (default=%(default)s)')

    parser.add_argument('--frag_size', '-k',
                        action='store_true',
                        dest='k',
                        default=23,
                        help="Fragment size to use (default=%(default)s)")
    return parser


def parse_fasta(fasta_file_name):
    '''
    Function to a parse a fasta file into a list of Seq object using SeqIO
    input: fasta_file_name
    output: iterator for SeqRecord objects
    '''
    if not os.path.isfile(fasta_file_name):
        raise IOError("{0} can't be found".format(fasta_file_name))

    seq_generator = Bio.SeqIO.parse(fasta_file_name, "fasta")
    return seq_generator

def create_bt2_index(reference_sequences_fn, r_seed=7):
    '''
    Function which takes in a fasta file name generates its bowtie2 index
    and returns the list of bowtie2 index filenames
    input: reference_sequences_fn file name of fasta file to use as reference
    output: list of filenames that make up index
    '''

    if not os.path.exists(reference_sequences_fn):
        raise ValueError('File does not exist: {0}'.format(reference_sequences_fn))

    # remove extension if it exists
    bt2_index_base = "eDicer_" + os.path.splitext(reference_sequences_fn)[0]
    devnull = open(os.devnull, 'w')
    bowtie2_build_command = 'bowtie2-build '\
                            '--seed {0} -f {1} {2}'.format(r_seed,
                                                           reference_sequences_fn,
                                                           bt2_index_base)

    # suppress stdout but maintain stderr defaul in case it fails
    ret_code = subprocess.call(bowtie2_build_command.split(),
                               stdout=devnull)

    if ret_code is not 0:
        raise OSError('Bowtie2 failed: {0}'.format(bowtie2_build_command))

    devnull.close()

    output_files = glob.glob(bt2_index_base + '.*')
    output_files.sort()

    return output_files


def generate_fragments(seqrec, k=23):
    '''
    Function to generate all possible contiguous k-length fragments of a
    specific seq
    input:   k - fragment size default is 23bp
           seq - seqrecord containing seq for fragmentation
    output: seq_fragments - list of seqrecords containing all the fragments
    '''

    #validate input
    if type(k) is not int or not int(k) > 0:
        raise ValueError("K must be an int > 0, "\
                         "'{0}' is thus invalid".format(str(k)))

    if type(seqrec) is not Bio.SeqRecord.SeqRecord:
        raise ValueError("seqrec must be a single SeqRecord obj, "\
                         "'{0}' is thus invalid".format(str(seqrec)))

    seq_len = len(seqrec.seq)

    if k > seq_len:
        # raise warning here if sequence is shorter than k
        warnings.warn("Sequence (seqrec={0}) shorter than k (k={1}) "\
                      "therefore discarding".format(seqrec.id, k),
                      UserWarning)
        return None

    elif seq_len >= k:
        #generate all fragments
        fragments = [Bio.SeqRecord.SeqRecord(seq=seqrec.seq[x: x + k],
                                             name=seqrec.name,
                                             id=seqrec.id,
                                             description=seqrec.description,
                                             dbxrefs=[]) \
                     for x in range(seq_len) if not x+k > seq_len]

        #add identifer information to fragment
        #this specified fragment size used and which fragment within the
        #set generated the current fragment is
        nfragments = len(fragments)
        labelled_fragments = [Bio.SeqRecord.SeqRecord(\
                                seq=fragments[x].seq,
                                id=seqrec.id,
                                name=fragments[x].name,
                                description=\
                                "{0} [{1}bp fragment {2} of {3}]"\
                                "".format(fragments[x].description,
                                          k,
                                          x+1,
                                          nfragments)) \
                              for x in range(nfragments)]

        return labelled_fragments

def write_fasta(seq_list, output_file):
    '''
    Output a list of SeqRecords to a file
    input:  seq_list
            output_file
    output: exit_status
    '''
    output_file = os.path.abspath(output_file)
    if not os.access(os.path.dirname(output_file), os.W_OK):
        os.utime(output_file, None)

    if not os.access(os.path.dirname(output_file), os.W_OK):
        raise IOError('{0} is not writeable'.format(output_file))

    with open(output_file, 'a') as out_fh:
        Bio.SeqIO.write(seq_list, out_fh, 'fasta')

def main(input_fasta, output_file, k=23):
    '''
    Main function to run the fragment generation if required
    '''

    for seq in parse_fasta(input_fasta):
        seq_fragments = generate_fragments(seq, k)
        write_fasta(seq_fragments, output_file)
