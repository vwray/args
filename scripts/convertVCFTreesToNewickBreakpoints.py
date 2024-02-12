# -*- coding: utf-8 -*-
"""
Creates breakpoints and Newicks file for existing VCF and trees file.

Created on Thu Feb  1 15:03:33 2024

@author: Valerie Wray
"""
import sys
import tskit

numberOfSamples = int(sys.argv[1])
inputVCFFile = sys.argv[2]
inputTreesFile = sys.argv[3]
outputNewickFile = sys.argv[4]
outputBreakpointsFile = sys.argv[5]

mts = tskit.load(inputTreesFile)

nodeLabels = {}
for i in range(4*numberOfSamples):
    nodeLabels[i] = str(i+1)

with open(outputNewickFile, 'w') as newickFile:
    for tree in mts.trees():
        newickFile.write(tree.as_newick(precision=3, node_labels=nodeLabels))
        newickFile.write('\n')
    
with open(outputBreakpointsFile, 'w') as breakpointsFile:
    breakpointsFile.write(str(list(mts.breakpoints())))
