# -*- coding: utf-8 -*-
"""
Converts Relate output to Newick strings and a list of positions.

Created on Tue Dec 5 17:59:07 2023

@author: Valerie Wray
"""

import tskit
import sys

numberOfSamples = sys.argv[1]
sequenceLength = sys.argv[2]
inputTreeFile = sys.argv[3]
outputNewickFile = sys.argv[4]
outputBreakpointsFile = sys.argv[5]


relateTs = tskit.load(inputTreeFile)

print(relateTs.first().draw_text())
print("Number of trees: ", len(relateTs.trees()))

nodeLabels = {}
for i in range(4*numberOfSamples):
    nodeLabels[i] = str(i+1)

with open(outputNewickFile, 'w') as newickFile:
    for tree in relateTs.trees():
        newickFile.write(tree.as_newick(precision=3, node_labels=nodeLabels))
        newickFile.write('\n')
    
with open(outputBreakpointsFile, 'w') as breakpointsFile:
    breakpointsFile.write(str(list(relateTs.breakpoints())))

print("Number of breakpoints: ", len(list(relateTs.breakpoints())))

