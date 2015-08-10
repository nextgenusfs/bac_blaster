#!/usr/bin/env python
 
import sys
 
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
 
import argparse
 
parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Fasta file')
parser.add_argument('qual', help='Qual file')
args = parser.parse_args()
 
records = PairedFastaQualIterator(
    open(args.fasta),
    open(args.qual)
)
for rec in records:
    sys.stdout.write(rec.format('fastq'))