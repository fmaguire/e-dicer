#!/usr/bin/env python

__version__ = '1.0.0'

import Bio
from Bio import SeqIO
import os.path
import glob
import logging
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
        raise IOError(f"{fasta_file_name} can't be found")

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
        logging.error(f"K must be an int > 0, '{k}' is thus invalid")
        sys.exit(1)

    if type(seqrec) is not Bio.SeqRecord.SeqRecord:
        logging.error("seqrec must be a single SeqRecord obj, '{seqrec}' is thus invalid")
        sys.exit(1)

    seq_len = len(seqrec.seq)

    if k > seq_len:
        # raise warning here if sequence is shorter than k
        logging.warn(f"Sequence (seqrec={seqrec.id}) shorter than k (k={k}) therefore discarding")
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
                                description=f"{fragments[x].description} [{k}bp fragment {x+1} of {nfragments}]") \
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
        logging.error(f"{output_file} is not writeable")
        sys.exit(1)

    with open(output_file, 'a') as out_fh:
        Bio.SeqIO.write(seq_list, out_fh, 'fasta')



def edicer(input_fasta, output_file, k=21):
    '''
    Main function to run the fragment generation if required
    '''

    for seq in parse_fasta(input_fasta):
        seq_fragments = generate_fragments(seq, k)
        if seq_fragments is not None:
            write_fasta(seq_fragments, output_file)


def collision_search(search_params):

    mismatches = search_params[0]
    name = search_params[1]
    folder = search_params[2]
    fasta = search_params[3]
    fragments = search_params[4]
    db = search_params[5]
    k = search_params[6]

    bowtie_index = fasta + '.index'

    if len(glob.glob(bowtie_index + '*')) == 0:
        cmd = "bowtie-build {fasta} {bowtie_index} > /dev/null 2>&1"
        logging.info(f"Building index {fasta}")
        exit = os.system(cmd)
        if exit != 0:
            logging.error(f"Failure: {cmd}")
            sys.exit(1)

    collisions = os.path.join(folder, name + '_collisions')
    collisions_summary = os.path.join(folder, name + '_summary')

    if not os.path.isfile(collisions) or not os.path.isfile(collisions_summary):
        # need to compute -v for 95% match
        logging.info(f"Mapping {fasta}")
        cmd = f"bowtie -f -v {mismatches} -a {bowtie_index} {fragments} > {collisions} 2> {collisions_summary}"

        exit = os.system(cmd)
        if exit != 0:
            logging.error(f"Failure: {cmd}")
            sys.exit(1)

        # get the kmer count statistics for query and db
        query_hashes = os.path.join(folder, name + '_query.jf')
        db_hashes = os.path.join(folder, name + '_db.jf')
        cmds = [f"jellyfish count -m {k} -s 50000 -o {query_hashes} {fragments}",
                f"echo 'query {fragments} {k}' >> {collisions_summary}",
                f"jellyfish stats {query_hashes} >> {collisions_summary}",
                f"jellyfish count -m {k} -s 50000 -o {db_hashes} {db}",
                'echo "\ndb {}" >> {}'.format(bowtie_index, collisions_summary),
                f"jellyfish stats {db_hashes} >> {collisions_summary}",
                f"rm {collisions} {query_hashes} {db_hashes}"]

        for cmd in cmds:
            exit = os.system(cmd)
            if exit != 0:
                logging.error(f"Failure: {cmd}")
                sys.exit(1)

    return collisions_summary

def summarise_results(results, run_dir):
    '''
    Summarise collision results
    '''
    data = {"db": [], "collision_count": [], 'db_kmer_count': [], 'query': [], 'k': [],
            'query_kmer_count': [], 'db_kmer_count_distinct': [], 'query_kmer_count_distinct': []}


    for fp in results:
        data['db'].append(os.path.split(fp)[-1].split("_summary")[0])


        with open(fp) as fh:
            for line in fh:
                if line.startswith('No alignments'):
                    data['collision_count'].append(0)
                elif line.startswith('Reported'):
                    data['collision_count'].append(line.split()[1])
                elif line.startswith('query'):
                    data['query'].append(line.split()[1])
                    data['k'].append(line.strip().split()[2])

                    # skip unique kmers
                    next(fh)
                    data['query_kmer_count_distinct'].append(int(next(fh).strip().split()[-1]))
                    data['query_kmer_count'].append(int(next(fh).strip().split()[-1]))
                elif line.startswith('db'):
                    next(fh)
                    data['db_kmer_count_distinct'].append(int(next(fh).strip().split()[-1]))
                    data['db_kmer_count'].append(int(next(fh).strip().split()[-1]))

    with open(os.path.join(run_dir, "collision_summary.pkl"), 'wb') as fh:
        pickle.dump(data, fh)

def check_if_valid_fasta(fp):
    """
    Checks if a file is a valid fasta file
    """
    try:
        parse = SeqIO.parse(fp, "fasta")
    except b:
        pass

def check_dependencies():
    """
    Check all dependencies exist and work
    """
    missing=False
    for program in ['bowtie --version', 'bowtie-build']:
        try:
            output = subprocess.run(program, shell=True, check=True,
                    stdout=subprocess.PIPE, encoding='utf-8')
            logging.debug(f"Tool {program} is installed: {output.stdout}")
        except:
            logging.error(f"Tool {program} is not installed")
            missing = True

    if missing:
        logging.error("One or more dependencies are missing please install")
        sys.exit(1)
    else:
        logging.debug("All dependencies found")


def run(args):
    '''
    Main runner function to dice and evaluate collisions against a db
    '''
    check_dependencies()

    # start logging
    if args.debug or args.verbose:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.DEBUG,
                            handlers=[logging.FileHandler(f"{run_name}.log"),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.INFO,
                            handlers=[logging.FileHandler(f"{run_name}.log"),
                                      logging.StreamHandler()])

    logging.info(f"Started ETD '{run_name}' with input '{args.input_genome}'")








    input_fp = os.path.abspath(args.query)
    db_folder = os.path.abspath(args.database_dir)


    # check
    dbfiles = list[glob.glob(os.path.join(db_folder, '*.fas'))]
    dbfiles += list[glob.glob(os.path.join(db_folder, "*.fasta"))]
    dbfiles += list[glob.glob(os.path.join(db_folder, "*.fna"))]
    dbfiles += list[glob.glob(os.path.join(db_folder, "*.fa"))]
    if len(dbfiles) == 0:
        logging.error(f"Cannot find any fasta files in {args.database_dir}")
        sys.exit(1)


    run_dir = f"{args.run_name}_k{args.frag_size}_output"

    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)


    fragments = os.path.join(run_dir,
                             os.path.split(input_fp)[-1] + "_fragments")

    if not os.path.isfile(fragments):
        print(f"eDicing {input_fp}")
        edicer(input_fp, fragments, args.frag_size)

    db_dir = os.path.join(run_dir, "db")
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)

    mismatches = int(args.frag_size * args.mismatches)

    search_params = []

    for db in dbfiles:
        name = os.path.splitext(os.path.split(db)[-1])[0]

        folder = os.path.join(db_dir, name)

        if not os.path.isdir(folder):
            os.mkdir(folder)

        fasta = os.path.join(folder, name)

        if not os.path.isfile(fasta):
            shutil.copyfile(db, fasta)

        search_params.append((mismatches, name, folder, fasta, fragments, db, k))

    results = Parallel(n_jobs=args.num_threads)(delayed(collision_search)(x) for x in search_params)

    summarise_results(results, run_dir)


