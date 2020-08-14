#!/usr/bin/env pythom
#-*- coding: utf-8 -*-

from edicer import edicer, utils
import argparse
import os


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Detect potential sRNA '
                                                 'collisions between a query '
                                                 'and folder of reference '
                                                 'sequences',
                                     prog='eDicer')


    parser.add_argument('-v', '--version', action='version',
                        version=f"%(prog)s {edicer.__version__}")

    parser.add_argument('-q', '--query',
                        type=lambda x: utils.is_valid_file(parser, x),
                        required=True,
                        help="Path to fasta containing query DNA ORF sequences")

    parser.add_argument('-d', '--database_dir',
                        type=lambda x: utils.is_valid_file(parser, x),
                        required=True,
                        help="Folder containing one or more fasta files of DNA "
                            " to use as the reference databases")

    parser.add_argument('-k', '--frag_size',
                        type=int,
                        required=True,
                        help="sRNA fragment size to use e.g. k=21")

    parser.add_argument('-n', '--run_name',
                        default=False,
                        help="Name/path to store edicer output")

    parser.add_argument('-m', '--mismatches',
                        type=float,
                        default=0.05,
                        help="Proportion of mismatches allowed")

    parser.add_argument('-j', '--num_threads', default=1, type=int,
                        help="Number of threads to use")

    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Run with verbose output")

    parser.add_argument('--overwrite', action='store_true', default=False,
                        help="Overwrite previous output directory")

    parser.add_argument('--quiet', action='store_true', default=False,
                        help="Run with minimal output")

    parser.add_argument('--keep_tmp', action='store_true', default=False,
                        help="Keep all intermediate outputs")

    args = parser.parse_args()

    edicer.run(args)
