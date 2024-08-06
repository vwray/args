# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 17:59:07 2023

Converts Relate output to newick strings and breakpoints.

@author: Valerie Wray
"""

import tskit

numberOfSamples = 3
sequenceLength = '100_000'

relateTs = tskit.load("H:\\My Drive\\Fall2023\\research\\relateScripts\\" + sequenceLength + "\\relateOutput_2.trees")

print(relateTs.first().draw_text())
print("Number of trees: ", len(relateTs.trees()))

nodeLabels = {}
for i in range(4*numberOfSamples):
    nodeLabels[i] = str(i+1)

with open("H:\\My Drive\\Fall2023\\research\\relateScripts\\" + sequenceLength + "\\relateNewick.txt", 'w') as newickFile:
    for tree in relateTs.trees():
        #print(tree.as_newick(precision=3, node_labels=nodeLabels))
        #print(tree.draw_text())
        newickFile.write(tree.as_newick(precision=3, node_labels=nodeLabels))
        newickFile.write('\n')
    
with open("H:\\My Drive\\Fall2023\\research\\relateScripts\\" + sequenceLength + "\\relateBreakpoints.txt", 'w') as breakpointsFile:
    breakpointsFile.write(str(list(relateTs.breakpoints())))

print("Number of breakpoints: ", len(list(relateTs.breakpoints())))

