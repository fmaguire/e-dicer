eDicer
=======

[![Build Status](https://travis-ci.org/fmaguire/eDicer.svg)](https://travis-ci.org/fmaguire/eDicer)

Detect potential sRNA collisions between a query and folder of reference sequences. This simulates the generation of Dicer fragments from a transcriptome of a given k-mer size (e.g. 21-mers) and searches for all collisions in a set of references with a certain amount of mismatch permitted.

Installation
------------

Dependencies can be installed using conda:

`conda env create -f env.yml`

Or manually installed (tested with following versions although older may be compatible) with external dependecies:

```
  - kmer-jellyfish=2.3.0
  - bowtie=1.3.0
```

and python libraries:
```
  - biopython=1.77
  - pandas=1.1.0
  - joblib=0.16.
```

Then just run `pip install .` within the directory to install the library.

Usage
-------

```
> python edicer.py -h

usage: eDicer [-h] [-v] -q QUERY -d DATABASE_DIR -k FRAG_SIZE [-n RUN_NAME]
              [-m MISMATCHES] [-j NUM_THREADS] [--verbose] [--overwrite]
              [--quiet] [--keep_tmp]

Detect potential sRNA collisions between a query and folder of reference
sequences

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -q QUERY, --query QUERY
                        Path to fasta containing query DNA ORF sequences
  -d DATABASE_DIR, --database_dir DATABASE_DIR
                        Folder containing one or more fasta files of DNA to
                        use as the reference databases
  -k FRAG_SIZE, --frag_size FRAG_SIZE
                        sRNA fragment size to use e.g. k=21
  -n RUN_NAME, --run_name RUN_NAME
                        Name/path to store edicer output
  -m MISMATCHES, --mismatches MISMATCHES
                        Proportion of mismatches allowed
  -j NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads to use
  --verbose             Run with verbose output
  --overwrite           Overwrite previous output directory
  --quiet               Run with minimal output
  --keep_tmp            Keep all intermediate outputs
```

If your file layout is:

```
├── sample_cds_seqs.fasta
└── test_references
    ├── test1.fasta
    └── test2.fasta
```

eDicer is invoked:

```python edicer.py -q test/sample_cds_seqs.fasta -d test/test_references -k 21```

This will generate an automatically named output folder:

```
sample_cds_seqs_k21_output.log
sample_cds_seqs_k21_output
├── database_collisions
│   ├── test1
│   │   ├── test1                      # copy of reference fasta
│   │   ├── test1_collisions.sam       # all collisions between query and test1 reference
│   │   └── test1_summary.json         # json summary of collisions and k-mer counts (unique, distinct, total) of test 1 reference
│   └── test2                          
│       ├── test2
│       ├── test2_collisions.sam
│       └── test2_summary.json
├── run_summary.json                   # combined summary of run in json
├── run_summary.tsv                    # and tsv format
├── sample_cds_seqs_fragments.fna      # generated query sequence fragments/k-mers
└── sample_cds_seqs_stats.json         # summary of query k-mer counts (unique, distinct, total)
```


To combine results from several queries into one table:

```
> python edicer-summarise.py -h                                                                                                                                                                                                                                            

usage: edicer-summarise [-h] [-v] [-i SUMMARY_JSONS [SUMMARY_JSONS ...]]
                        [-o OUTPUT]

Collate outputs from multiple eDicer runs into one file

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i SUMMARY_JSONS [SUMMARY_JSONS ...], --summary_jsons SUMMARY_JSONS [SUMMARY_JSONS ...]
                        List of run_summary.json files
  -o OUTPUT, --output OUTPUT
                        Combined summary output location
```

So if we ran the same query with different kmer sizes in our test data:

```python edicer-summarise.py -i sample_cds_seqs_k21_output/run_summary.json sample_cds_seqs_k22_output/run_summary.json```

This script will generate a single combined `combined_output.tsv` and `combined_output.json`

