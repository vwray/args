# -*- coding: utf-8 -*-
"""
A script to convert VCF (Variant Calling Format) files to the format 
required as input for RENT+, has a first row of positions followed by rows
of 0's and 1's.

Created on Tue Nov  7 10:32:18 2023

@author: Valerie Wray
"""

import csv
import sys

vcfInputFile = sys.argv[1]
outputFile = sys.argv[2]

with open(vcfInputFile, 'r') as fin, open(outputFile, 'w') as fout:
    tsv_reader = csv.reader(fin, delimiter="\t")
    for i in range(5):
        next(tsv_reader)
    header = next(tsv_reader)
    
    positions = []
    data = []
    
    for row in tsv_reader:
        positions.append(row[1])
        newRow = []
        for i in range(9,len(row)):
            newRow.extend(row[i].split('|'))
          
        data.append(newRow)
    fout.write(' '.join(positions)+'\n')
    for i in range(len(data[0])):
        for j in range(len(data)):
            fout.write(data[j][i])
        fout.write('\n')  