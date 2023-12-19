# -*- coding: utf-8 -*-
"""
Converts ARGweaver consensus trees created after running ARGweaver 
and then using arg-cons to determine the consensus tree topologies.

Created on Tue Dec 12 09:02:55 2023

@author: Valerie Wray
"""

import csv
import sys

numberOfSamples = sys.argv[1]
sequenceLength = sys.argv[2]
inputFile = sys.argv[3]
outputFile = sys.argv[4] 
positionsFile = sys.argv[5] 

with open(inputFile, 'r') as fin, open(outputFile, 'w') as fout, open(positionsFile, 'w') as pout:
    tsv_reader = csv.reader(fin, delimiter="\t")
    next(tsv_reader)
    next(tsv_reader)
    positions=[]
    
    for row in tsv_reader:
        positions.append(int(row[1]))
        newickString = row[3]
        fout.write(newickString + '\n')
    
    pout.write(str(positions))
