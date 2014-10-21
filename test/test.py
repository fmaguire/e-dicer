#!/usr/bin/env python

import unittest
import e_dicer

import os
import sys

import types
import warnings
import random

from StringIO import StringIO
from Bio import Alphabet
from Bio import SeqRecord
from Bio import Seq
from Bio import SeqIO

class TestFastaParsing(unittest.TestCase):
    '''
    Class to test program is correctly parsing fasta files
    '''

    def setUp(self):
        '''
        Initialise data for tests by finding
        '''
        self.sample_fasta_file = u'test/sample_cds_seqs.fasta'
        self.seq_number = 100

    def test_exception_raising(self):
        '''
        Test if fasta parser correctly raises IO error for non-existant fasta file
        '''
        self.assertRaises(IOError, e_dicer.parse_fasta, u'NOT_A_REAL_FASTA.fasta')

    def test_generator_creation(self):
        '''
        Test for correct creation of SeqIO generator
        '''
        self.seq_gen = e_dicer.parse_fasta(self.sample_fasta_file)
        self.assertIsInstance(self.seq_gen, types.GeneratorType)

    def test_parsed_sequences(self):
        '''
        Test if generator contains correct number of sequences and
        Test if
        '''
        self.seqs = [seq for seq in \
                     e_dicer.parse_fasta(self.sample_fasta_file)]

        self.assertIsInstance(self.seqs, list)
        self.assertEqual(len(self.seqs), 100)

        self.random_seq = random.sample(self.seqs, 1)[0]
        self.assertIs(type(self.random_seq), SeqRecord.SeqRecord)

    #todo: implement test (and feature) to ensure parsed input is nucl
    #      and not prot

if __name__=='__main__':

    unittest.main()

