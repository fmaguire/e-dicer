#!/usr/bin/env python
#-*- coding: utf-8 -*-

from edicer import edicer
import argparse
import os

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Detect potential sRNA '
                                                 'collisions between a query '
                                                 'and folder of reference '
                                                 'sequences',
                                     prog='eDicer')


    parser.add_argument('-v', '--version', action='version',
                        version=f"%(prog)s {edicer.__version__}")

    parser.add_argument('-q', '--query',
                        type=lambda x: is_valid_file(parser, x),
                        required=True,
                        help="Path to fasta containing query DNA ORF sequences")

    parser.add_argument('-d', '--database_dir',
                        type=lambda x: is_valid_file(parser, x),
                        required=True,
                        help="Folder containing one or more fasta files of DNA "
                            " to use as the reference databases")

    parser.add_argument('-k', '--frag_size',
                        type=int,
                        required=True,
                        help="sRNA fragment size to use e.g. k=21")

    parser.add_argument('-n', '--run_name',
                        default="eDicer_results",
                        help="Name for output")

    parser.add_argument('-m', '--mismatches',
                        type=float,
                        default=0.05,
                        help="Proportion of mismatches allowed")

    parser.add_argument('-j', '--num_threads', default=1, type=int,
                        help="Number of threads to use")

    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Run with verbose output")

    args = parser.parse_args()

    edicer.run(args)
