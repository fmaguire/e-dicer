#!/usr/bin/env python

import Bio
from Bio import SeqIO
import os.path
import warnings
import glob
import sys
import os
import shutil
from joblib import Parallel, delayed
import multiprocessing
import pickle

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

def write_fasta(seq_list, output_file):
    '''
    Output a list of SeqRecords to a file
    input:  seq_list
            output_file
    output: exit_status
    '''

    if not os.access(os.path.dirname(output_file), os.W_OK):
        raise IOError('{0} is not writeable'.format(output_file))

    with open(output_file, 'a') as out_fh:
        Bio.SeqIO.write(seq_list, out_fh, 'fasta')



def edicer(input_fasta, output_file, k=21):
    '''
    Main function to run the fragment generation if required
    '''

    for seq in parse_fasta(input_fasta):
        seq_fragments = generate_fragments(seq, k)
        write_fasta(seq_fragments, output_file)


def collision_search(search_params):

    mismatches = search_params[0]
    name = search_params[1]
    folder = search_params[2]
    fasta = search_params[3]
    fragments = search_params[4]
    db = search_params[5]

    bowtie_index = fasta + '.index'

    if len(glob.glob(bowtie_index + '*')) == 0:
        cmd = "bowtie-build {} {} > /dev/null 2>&1".format(fasta, bowtie_index)
        print("Building index {}".format(fasta))
        exit = os.system(cmd)
        if exit != 0:
            print("Failure: {}".format(cmd))
            sys.exit(1)

    collisions = (os.path.join(folder, name + '_collisions'),
                  os.path.join(folder, name + '_summary'))

    if not os.path.isfile(collisions[0]) or not os.path.isfile(collisions[1]):
        # need to compute -v for 95% match
        print("Mapping {}".format(fasta))
        cmd = "bowtie -f -v {} -a {} {} > {} 2> {}".format(mismatches,
                                                           bowtie_index,
                                                           fragments,
                                                           collisions[0],
                                                           collisions[1])
        exit = os.system(cmd)
        if exit != 0:
            print("Failure: {}".format(cmd))
            sys.exit(1)

        cmd = 'grep -v "^>" {} | wc -m >> {}'.format(db, collisions[1])
        exit = os.system(cmd)

    return collisions[1]

def summarise_results(results, run_dir):
    '''
    Summarise collision results
    '''
    data = {"source": [], "collision_count": [], 'size':[]}
    for fp in results:
        query = os.path.split(fp)[-1].split("_summary")[0]
        with open(fp) as fh:
            output = [x.strip() for x in fh.readlines()]
            alignments = output[-2]
            size = int(output[-1])

            if alignments.startswith("No alignments"):
                count = 0
            else:
                count = int(alignments.split()[1])
        data['source'].append(query)
        data['collision_count'].append(count)
        data['size'].append(size)

    with open(os.path.join(run_dir, "collision_summary.pkl"), 'wb') as fh:
        pickle.dump(data, fh)


def main(args):
    '''
    Main runner function to dice and evaluate collisions against a db
    '''

    # check arguments

    if len(args) != 5:
        print("Requires 4 arguments : run.py QUERY_FASTA DB_FOLDER K RUN_NAME")
        sys.exit(1)

    input_fp = os.path.abspath(args[1])
    db_folder = os.path.abspath(args[2])
    k = int(args[3])
    run_name = args[4]

    if not os.path.isfile(input_fp):
        print("{} does not exist".format(input_fp))
        sys.exit(1)

    dbfiles = glob.glob(os.path.join(db_folder, '*.fas'))
    if len(dbfiles) == 0:
        print("Cannot find any db files in {} (need .fas) suffix".format(db_folder))
        sys.exit(1)


    run_dir = '{}_k{}_output'.format(run_name, k)

    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)


    fragments = os.path.join(run_dir,
                             os.path.split(input_fp)[-1] + "_fragments")

    if not os.path.isfile(fragments):
        print("eDicing {}".format(input_fp))
        edicer(input_fp, fragments, k)

    db_dir = os.path.join(run_dir, "db")
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)

    mismatches = int(k * 0.05)

    search_params = []

    for db in dbfiles:
        name = os.path.splitext(os.path.split(db)[-1])[0]

        folder = os.path.join(db_dir, name)

        if not os.path.isdir(folder):
            os.mkdir(folder)

        fasta = os.path.join(folder, name)

        if not os.path.isfile(fasta):
            shutil.copyfile(db, fasta)

        search_params.append((mismatches, name, folder, fasta, fragments, db))

    num_cores = multiprocessing.cpu_count()

    results = Parallel(n_jobs=num_cores)(delayed(collision_search)(x) for x in search_params)

    summarise_results(results, run_dir)


