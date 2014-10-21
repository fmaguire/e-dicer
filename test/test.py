#!/usr/bin/env python

import unittest
import e_dicer

import os
import sys

import types
import warnings
import random

import Bio
from Bio import SeqIO
import StringIO

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
        self.assertRaisesRegexp(IOError ,"NOT_A_REAL_FASTA.fasta can't be found",
                                eDicer.parse_fasta, u'NOT_A_REAL_FASTA.fasta')

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
        self.assertIs(type(self.random_seq), Bio.SeqRecord.SeqRecord)

    #todo: implement test (and feature) to ensure parsed input is nucl
    #      and not prot


class TestDicerFunction(unittest.TestCase):
    '''
    Class to test the correct generation of all possible k-fragments
    from a single SeqRecord
    '''

    def setUp(self):
        '''
        Parse and grab a random seq
        '''
        self.raw_seq = """>m.32354 g.32354  ORF g.32354 m.32354 \
type:3prime_partial len:74 (+) \
comp10036_c0_seq1:394-618(+)_
ATGCAGATGCGTGTTCGAATGTACGTTGGTGGCG\
TGCTGCTCGCGCTTCTGCTTGCGA"""
        self.example_seq = Bio.SeqIO.parse(StringIO.StringIO(self.raw_seq), "fasta").next()

        self.example_21_frag_file = "test/sample_fragments.fasta"

        self.example_frags = [seqrecord for seqrecord in \
                             Bio.SeqIO.parse(self.example_21_frag_file, "fasta")]

        self.example_seq_len = len(self.example_seq.seq)

    def test_sequence_too_short(self):
        '''
        Test when len(seq) < k function correctly discards seq with a warning
        '''

        with warnings.catch_warnings(record=True) as w:
            self.short = eDicer.generate_fragments(self.example_seq,
                                                   self.example_seq_len + 1)

            self.assertIs(self.short, None)

            self.assertEqual(len(w), 1)
            self.assertIs(w[-1].category, UserWarning)
            self.assertEqual(str(w[-1].message), "Sequence (seqrec={0}) shorter than "\
                                            "k (k={1}) therefore "\
                                            "discarding".format(self.example_seq.id,
                                                                self.example_seq_len+1))

    def test_sequence_same_length(self):
        '''
        Test when k==len(seq) function just returns the same sequence
        '''
        self.same = eDicer.generate_fragments(self.example_seq,
                                              self.example_seq_len)

        self.assertIsNot(self.same, None)
        self.assertEqual(len(self.same), 1)
        self.assertEqual(str(self.same[0].seq), str(self.example_seq.seq))
        self.assertEqual(self.same[0].description, "{0} [{1}bp fragment" \
                                                   " 1 of 1]".format(self.example_seq.description,
                                                                     self.example_seq_len))


    def test_sequence_fragmentation(self):
        '''
        Test for correct fragmentation of sample seq with default k (21)
        '''

        test_frags = eDicer.generate_fragments(self.example_seq)
        for frag_index in range(len(self.example_frags)):
            self.assertEqual(str(self.example_frags[frag_index].seq),
                             str(test_frags[frag_index].seq))
            self.assertEqual(self.example_frags[frag_index].description,
                             test_frags[frag_index].description)


    def test_incorrect_input(self):
        '''
        Test function raises exception if input isn't a seqrec and positive int
        '''

        self.randint = random.randint(-1000, 0)

        self.assertRaisesRegexp(ValueError,
                                "K must be an int > 0, '{0}' is thus "\
                                "invalid".format(str(self.randint)),
                                eDicer.generate_fragments, self.example_seq, self.randint)

        self.assertRaisesRegexp(ValueError,
                                "K must be an int > 0, 'a' is thus invalid",
                                eDicer.generate_fragments, self.example_seq, 'a')

        self.assertRaisesRegexp(ValueError,
                                "seqrec must be a single SeqRecord obj, 'test'"\
                                " is thus invalid",
                                eDicer.generate_fragments, 'test', 54)

if __name__=='__main__':

    unittest.main()
