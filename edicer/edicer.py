#!/usr/bin/env python

__version__ = '1.0.0'

from Bio import SeqIO, SeqRecord
import glob, sys, os, logging, json, shutil
from joblib import Parallel, delayed
import pandas as pd

from edicer import utils


def generate_fragments(seqrec, k=21):
    '''
    Generate all possible contiguous k-length fragments of a
    specific seq
    input:
        - k: fragment size (default is 21bp)
        - seq: seqrecord containing seq for fragmentation
    output:
        - seq_fragments: list of seqrecords for all fragments from this sequence
    '''

    #validate input
    if type(k) is not int or not int(k) > 0:
        logging.error(f"K must be an int > 0, '{k}' is thus invalid")
        sys.exit(1)

    if type(seqrec) is not SeqRecord.SeqRecord:
        logging.error("seqrec must be a single SeqRecord obj, '{seqrec}' is thus invalid")
        sys.exit(1)

    seq_len = len(seqrec.seq)

    if k > seq_len:
        # raise warning here if sequence is shorter than k
        logging.warning(f"Sequence (seqrec={seqrec.id}) shorter than k (k={k}) therefore discarding")
        return None

    elif seq_len >= k:
        #generate all fragments
        fragments = [SeqRecord.SeqRecord(seq=seqrec.seq[x: x + k],
                                             name=seqrec.name,
                                             id=seqrec.id,
                                             description=seqrec.description,
                                             dbxrefs=[]) \
                     for x in range(seq_len) if not x+k > seq_len]

        #add identifer information to fragment
        #this specified fragment size used and which fragment within the
        #set generated the current fragment is
        nfragments = len(fragments)
        labelled_fragments = [SeqRecord.SeqRecord(\
                                seq=fragments[x].seq,
                                id=seqrec.id,
                                name=fragments[x].name,
                                description=f"{fragments[x].description} [{k}bp fragment {x+1} of {nfragments}]") \
                              for x in range(nfragments)]

        return labelled_fragments


def edicer(input_fasta, output_file, k=21):
    '''
    Dice all sequences in the query fasta into their k-size sequence fragments
    input:
        - input_fasta: path to the query fasta file
        - output_fasta: path to output the fasta file of k-size fragments
        - k: fragment size to use (default k=21)
    '''

    for seq in SeqIO.parse(input_fasta, "fasta"):
        seq_fragments = generate_fragments(seq, k)
        if seq_fragments is not None:
            with open(output_file, 'a') as out_fh:
                SeqIO.write(seq_fragments, out_fh, 'fasta')


def collision_search(search_params):
    """
    Use bowtie to search for shared k-size sequences between the generated
    query fragments and each of the database fasta files

    input:
        - search_params: dictionary containing all the necessary search params
    output:
        - collisions_summary: path to collision stats file
    """

    # check if index already exists or not for this sequence
    # if not generate it using bowtie
    bowtie_index = search_params['db_fasta'] + '.index'
    if len(glob.glob(bowtie_index + '*')) == 0:
        cmd = f"bowtie-build {search_params['db_fasta']} {bowtie_index}"
        logging.debug(f"Building index {search_params['db_fasta']}")
        logging.debug(f"Running command: {cmd}")
        output = subprocess.run(cmd, shell=True, check=True,
                        stdout=subprocess.PIPE, encoding='utf-8')
    else:
        logging.debug(f"Database index already exists, ensure it was built "
            f"using same version of bowtie, using: {bowtie_index}")


    # check if the collisions output files exist if not run bowtie to generate
    # them and parse output to build stats file
    collisions = os.path.join(search_params['db_output_folder'],
                              search_params['database_name'] + '_collisions.sam')
    collision_stats = os.path.join(search_params['db_output_folder'],
                                   search_params['database_name'] + \
                                           '_summary.json')

    if not os.path.isfile(collisions) or not os.path.isfile(collision_stats):
        # need to compute -v for 95% match
        logging.debug(f"Mapping {search_params['query_fragment_fp']} against "
                     f"{search_params['db_fasta']} with "
                     f"{search_params['frag_size']}bp fragments and "
                     f"{search_params['mismatches']} mismatches")
        cmd = f"bowtie -f -v {search_params['mismatches']} -a {bowtie_index} \
                {search_params['query_fragment_fp']} > {collisions}"
        logging.debug(f"Running command: {cmd}")
        output = subprocess.run(cmd, shell=True, check=True,
                        stderr=subprocess.PIPE, encoding='utf-8')
        output = output.stderr.split('\n')

        # parse collision summary
        logging.debug(f"Parsing collision summary statistics: {output}")
        collision_stats_data = {"db_fasta": search_params['db_fasta'],
                                'db_name': search_params['database_name'],
                                'fragment_size': search_params['frag_size'],
                                'query_fragments': search_params['query_fragment_fp'],
                                'mismatches': search_params['mismatches']}

        for line in output:
            if line.startswith('No alignments'):
                collision_stats_data['collision_count'] = 0
            elif line.startswith('Reported'):
                collision_stats_data['collision_count'] = line.split()[1]
            elif line.startswith('query'):
                collision_stats_data['query'] = line.split()[1]
                collision_stats_data['fragment_size'] = line.strip().split()[2]

        # get k-mer statistics for database
        logging.info("Calculating database fragment statistics: "
                     f"{search_params['db_fasta']}")
        db_hashes = os.path.join(search_params['db_output_folder'],
                                 search_params['database_name'] + '_db.jf')
        count_cmd = f"jellyfish count -m {search_params['frag_size']} -s 50000 -o {db_hashes} {search_params['db_fasta']}"
        logging.debug(f"Running command: {count_cmd}")
        output = subprocess.run(count_cmd, shell=True, check=True,
                                stdout=subprocess.PIPE, encoding='utf-8')

        stats_cmd = f"jellyfish stats {db_hashes}"
        logging.debug(f"Running command: {stats_cmd}")
        output = subprocess.run(stats_cmd, shell=True, check=True,
                                stdout=subprocess.PIPE, encoding='utf-8')
        output = output.stdout.split('\n')
        logging.debug(f"Parsing output from '{stats_cmd}': {output}")
        collision_stats_data['db_kmer_count_distinct'] = int(output[1].strip().split()[-1])
        collision_stats_data['db_kmer_count'] = int(output[2].strip().split()[-1])



        logging.debug(f"Writing to {collision_stats}: {collision_stats_data}")
        with open(collision_stats, 'w') as out_fh:
            json.dump(collision_stats_data, out_fh)
    else:
        logging.debug("Collision data exists using: "
                      f"{collisions} and {collision_stats}")
        with open(collision_stats) as fh:
            collision_stats_data = json.load(fh)

    return collision_stats_data


def run(args):
    '''
    Main runner function to dice and evaluate collisions against a db
    '''
    query_name = os.path.basename(os.path.splitext(args.query)[0])
    if not args.run_name:
        run_name = f"{query_name}_k{args.frag_size}_output"
    else:
        run_name = f"{args.run_name}_k{args.frag_size}_output"

    # start logging
    if args.verbose:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.DEBUG,
                            handlers=[
                                logging.FileHandler(f"{run_name}.log"),
                                logging.StreamHandler()])
    elif args.quiet:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.ERROR,
                            handlers=[
                                logging.FileHandler(f"{run_name}.log"),
                                logging.StreamHandler()])
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s',
                            level=logging.INFO,
                            handlers=[
                                logging.FileHandler(f"{run_name}.log"),
                                logging.StreamHandler()])


    logging.info(f"Started edicer: '{run_name}'")

    # check all depencies and installed and working
    utils.check_dependencies()

    # check input fasta is valid
    utils.is_valid_fasta(args.query)

    # check reference fasta files exist and are valid fasta files
    allowed_fasta_extensions = ["*.fas", "*.fasta", "*.fna", "*.fa"]
    dbfiles = []
    logging.debug(f"Searching '{args.database_dir}' for files ending in "
                  f"'{allowed_fasta_extensions}'")
    for ext in allowed_fasta_extensions:
        dbfiles.extend(glob.glob(os.path.join(args.database_dir, ext)))

    if len(dbfiles) == 0:
        logging.error(f"Cannot find any fasta files in {args.database_dir}: "
                      f" {allow_fasta_extensions}")
        sys.exit(1)
    for db_fasta in dbfiles:
        utils.is_valid_fasta(db_fasta)

    logging.info(f"Detecting {(1-args.mismatches) * 100}% {args.frag_size}bp "
            f"matches between '{args.query}' and reference sequences: "
            f"'{dbfiles}' in '{run_name}'")

    # make the output dir/or continue on from a previous analysis
    logging.debug(f"Making output directory: {run_name}")
    if os.path.isdir(run_name):
        if args.overwrite:
            logging.debug(f"Removing and remaking (--overwrite used): {run_name}")
            shutil.rmtree(run_name)
            os.mkdir(run_name)
        else:
            logging.warning(f"Previous analysis found, resuming: '{run_name}' "
                             "(use debugging for details) ")
    else:
        os.mkdir(run_name)

    # cut the input query sequences into k-size fragments for comparison with
    # databases
    logging.info(f"Dicing sequences in query: {args.query}")
    fragment_file = os.path.join(run_name, query_name + "_fragments.fna")
    if not os.path.isfile(fragment_file):
        logging.info(f"Splitting {args.query} in {args.frag_size}bp fragments")
        edicer(args.query, fragment_file, args.frag_size)
    else:
        logging.debug(f"Query fragments exists, using: {fragment_file}")

    # getting k-mer count statistics on query fragments
    logging.info(f"Calculating k-mer statistics on {args.query}")
    query_statistics_file = os.path.join(run_name,
                                       query_name + "_stats.json")

    # if stats file doesn't exist use jellyfish to calculate statistics
    if not os.path.exists(query_statistics_file):
        query_hashes = os.path.join(run_name, query_name  + "_fragment_hashes.jf")
        count_cmd = f"jellyfish count -m {args.frag_size} -s 50000 -o {query_hashes} {args.query}"
        logging.debug(f"Running command: {count_cmd}")
        output = subprocess.run(count_cmd, shell=True, check=True,
                        stdout=subprocess.PIPE, encoding='utf-8')

        stats_cmd = f"jellyfish stats {query_hashes}"
        logging.debug(f"Running command: {stats_cmd}")
        output = subprocess.run(stats_cmd, shell=True, check=True,
                        stdout=subprocess.PIPE, encoding='utf-8')

        logging.debug(f"Parsing output from {stats_cmd}: {output.stdout}")
        output = output.stdout.split('\n')
        query_stats = {'query': args.query,
                       'fragment_file': fragment_file,
                       'fragment_size': args.frag_size,
                       'query_kmer_count_distinct': int(output[1].strip().split()[-1]),
                       'query_kmer_count': int(output[2].strip().split()[-1])}
        logging.debug(f"Writing to {query_statistics_file}: {query_stats}")
        with open(query_statistics_file, 'w') as out_fh:
            json.dump(query_stats, out_fh)

        #if not args.keep_temp:
        #    logging.debug(f"Deleting {query_hashes}")
        #    os.remove(query_hashes)
    else:
        logging.debug(f"Query statistics exists, using: {query_statistics_file}")
        with open(query_statistics_file) as fh:
            query_stats = json.load(fh)


    logging.debug("Making directory for per database outputs")
    db_dir = os.path.join(run_name, "database_collisions")
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)
    else:
        logging.debug(f"Database output directory exists, using: {db_dir}")

    # number of base mismatches permitted is rounded down to nearest integer
    # i.e. 21 fragments with 5% permitted mismatches = 1 permitted non-matching
    # base
    mismatches = int(args.frag_size * args.mismatches)

    search_params = {}
    logging.info("Finding collisions between query and each database")

    # to make this more efficient for comparing against a large number of
    # reference genomes we want to run it the alignments in parallel
    all_search_params = []
    for db in dbfiles:
        db_search_params = {'mismatches': mismatches,
                            'frag_size': args.frag_size,
                            'query_fragment_fp': fragment_file}

        db_search_params['database_name'] = os.path.splitext(\
                                                    os.path.basename(db))[0]

        db_search_params['db_output_folder'] = os.path.join(db_dir,
                                            db_search_params['database_name'])


        if not os.path.isdir(db_search_params['db_output_folder']):
            os.mkdir(db_search_params['db_output_folder'])
        else:
            logging.debug(f"Individual database output directory exists, "
                            f"using: {db_search_params['db_output_folder']}")


        db_search_params['db_fasta'] = os.path.join(\
                                        db_search_params['db_output_folder'],
                                        db_search_params['database_name'])

        if not os.path.isfile(db_search_params['db_fasta']):
            logging.debug(f"Copying database fasta {db} to {db_search_params['db_fasta']})")
            shutil.copyfile(db, db_search_params['db_fasta'])
        else:
            logging.debug(f"Database fasta exists, using: {db_search_params['db_fasta']})")


        all_search_params.append(db_search_params)

    results = Parallel(n_jobs=args.num_threads)(delayed(collision_search)\
                (db_search_params) for db_search_params in all_search_params)

    run_summary = os.path.join(run_name, "run_summary.json")
    logging.info(f"Summarising results to {run_summary}")
    if not os.path.exists(run_summary):
        summary_results = {'query': [], 'fragment_size': [], 'mismatch_proportion': [],
                           'query_kmer_count': [], 'query_kmer_count_distinct': [],
                           'collision_count': [], 'database': [], 'db_kmer_count': [],
                           'db_kmer_count_distinct': []}
        for result in results:
            summary_results['query'].append(query_name)
            summary_results['fragment_size'].append(args.frag_size)
            summary_results['mismatch_proportion'].append(args.mismatches)
            summary_results['query_kmer_count'].append(query_stats['query_kmer_count'])
            summary_results['query_kmer_count_distinct'].append(query_stats['query_kmer_count_distinct'])
            summary_results['collision_count'].append(result['collision_count'])
            summary_results['database'].append(result['db_name'])
            summary_results['db_kmer_count'].append(result['db_kmer_count'])
            summary_results['db_kmer_count_distinct'].append(result['db_kmer_count_distinct'])
        logging.debug(f"Dumping collated results: {summary_results}")
        with open(run_summary, 'w') as fh:
            json.dump(summary_results, fh)
    else:
        logging.debug(f"Overall results already exists: {run_summary}")
        with open(run_summary) as fh:
            summary_results = json.load(fh)

    run_summary_table = os.path.join(run_name, "run_summary.tsv")
    logging.info(f"Generating results table {run_summary_table}")
    pd.DataFrame(summary_results).to_csv(run_summary_table, sep='\t')
