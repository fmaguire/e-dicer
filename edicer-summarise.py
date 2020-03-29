#!/usr/bin/env python
#-*- coding: utf-8 -*-

from edicer import edicer, utils
import argparse
import os


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Collate outputs from multiple '
                                                 'eDicer runs into one file',
                                                 prog='edicer-summarise')


    parser.add_argument('-v', '--version', action='version',
                        version=f"%(prog)s {edicer.__version__}")

    parser.add_argument('-i', '--summary_jsons', #type=utils.is_valid_file,
                        nargs='+', help='List of run_summary.json files')

    parser.add_argument('-o', '--output', type=str, default="combined_output",
                         help="Combined summary output location")


    args = parser.parse_args()

    edicer.summarise(args)
