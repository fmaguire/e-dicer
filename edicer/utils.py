#!/usr/bin/env python

import subprocess
import logging
import sys
import os
from Bio import SeqIO

def check_dependencies():
    """
    Check all dependencies exist and work
    """
    missing=False
    for program in ['bowtie --version', 'bowtie-build -h',
                    'jellyfish --version']:
        try:
            output = subprocess.run(program, shell=True, check=True,
                        stdout=subprocess.PIPE, encoding='utf-8')
        except:
            logging.error(f"Tool {program} is not installed/working")
            missing = True
    if missing:
        logging.error("One or more dependencies are missing please install")
        sys.exit(1)
    else:
        logging.debug("All dependencies found")

def is_valid_fasta(fp):
    """
    Checks if a file is a valid fasta file
    """
    possible_fasta = SeqIO.parse(fp, "fasta")
    # if no valid records this returns False
    if not any(possible_fasta):
        logging.error(f"Fasta file contains no readable sequence records: {fp}")
        sys.exit(1)

def is_valid_file(parser, arg):
    """
    Checks whether input file exists

    Parameters:
        arg (str):The filename to be checked
    Returns:
        arg (str):The filename if it exists
    """
    if not os.path.exists(arg):
        parser.error(f"The path {arg} does not exist")
    else:
        return arg




