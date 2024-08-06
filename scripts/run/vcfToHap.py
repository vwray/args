# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 09:08:27 2023

@author: valer
"""

import csv

vcfInputFile = "H:\\My Drive\\Fall2023\\research\\outputVCFFile_10samples_1thousSequenceLen.vcf" #input("Enter vcf file name:")
hapOutputFile = "H:\\My Drive\\Fall2023\\research\\hapOutput1.hap" #input("Enter output file name:")

with open(vcfInputFile, 'r') as fin, open(hapOutputFile, 'w') as fout:
    tsv_reader = csv.reader(fin, delimiter="\t")
    #help(tsv_reader)
    for i in range(5):
        next(tsv_reader)
    header = next(tsv_reader)
        
    #first row
    

    for row in tsv_reader:
        for i in range(10,len(row)):
            #fout.write('>' + header[i] + '_' +  + '\n')
            #fout.write("{} {}".format(header[i], i%2 + 1))
            print("{} {}".format(header[i], i%2 + 1))
        print(type(row))
        help(fout.write)

#        (name, color, ranking) = row
#        print(f"{name} is rank {ranking}")
        
    #contents = fin.read()
    #print(contents)
    
    