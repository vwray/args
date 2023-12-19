# -*- coding: utf-8 -*-
"""
A script to convert VCF (Variant Calling Format) files to the sites format 
required as input for ARGweaver.

Created on Sun Nov 26 10:03:49 2023

@author: Valerie Wray
"""

import csv
import sys

vcfInputFile = sys.argv[1]
outputFile = sys.argv[2]

with open(vcfInputFile, 'r') as fin, open(outputFile, 'w') as fout:
    tsv_reader = csv.reader(fin, delimiter="\t")
    sequenceLength='0'
    for i in range(5):
        if(i==3):
            line = next(tsv_reader)
            line=line[0].split('=')
            sequenceLength=line[len(line)-1][:-1]
        else:
            next(tsv_reader)
    header = next(tsv_reader)
    population = header[9:]
    
    positions = []
    data = []
    
    numberOfSamples = len(population) * 2
    chrom = '1'
    
    
    fout.write("NAMES\tn" + "\tn".join(map(str,range(numberOfSamples))))
    first=True
    
    for row in tsv_reader:
        if(first):
            chrom = row[0]
            fout.write("\nREGION\tchr\t"+chrom+"\t"+sequenceLength + "\n")
            first=False
        
        ancestralState=row[3]
        mutatedState=row[4]
        newRow = []
        for i in range(9,len(row)):
            newRow.extend(row[i].split('|'))

        states=[]
        for i in newRow:
            if i=='0':
                states.append(ancestralState)
            else:
                states.append(mutatedState)
        
        #write the position followed by the states
        fout.write(row[1] + "\t" + ''.join(states) + '\n')
    