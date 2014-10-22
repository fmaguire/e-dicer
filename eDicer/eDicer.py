#!/usr/bin/env python

import Bio
from Bio import SeqIO
import os.path
import warnings

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

def generate_fragments(seqrec, k=21):
    '''
    Function to generate all possible contiguous k-length fragments of a
    specific seq
    input:   k - fragment size default is 21bp
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
