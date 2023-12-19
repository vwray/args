# -*- coding: utf-8 -*-
"""
A script to compare branch lengths and TMRCA of a true tree sequence created 
from msprime to estimated tree sequences created by RENT+, Relate, and ARGweaver.

Created on Wed Nov 22 11:55:14 2023

@author: valer
"""

import dendropy
import csv
from bisect import bisect_left
from compare_two_trees import compare_trees
from matplotlib import pyplot as plt
import sys

numberOfSamples = sys.argv[1]
sequenceLength = sys.argv[2]
trueTreesNewickFile = sys.argv[3]
trueTreesBreakpointsFile = sys.argv[4]
relateNewickFile = sys.argv[5]
relateBreakpointsFile = sys.argv[6]
rentTreesFile = sys.argv[7]
rentNewickFile = sys.arg[8]
argweaverNewickFile = sys.argv[9]
argweaverBreakpointsFile = sys.argv[10]
plotFile = sys.argv[11]

#True msprime trees
trueTreeList = dendropy.TreeList.get(path=trueTreesNewickFile, schema="newick")
print("number of true trees: ", len(trueTreeList))

#read in msprime breakpoints
with open(trueTreesBreakpointsFile, 'r') as breakpointsFile:
    breakpoints = breakpointsFile.readline()
breakpoints = breakpoints[1:len(breakpoints)-1].split(', ')
breakpoints = [int(float(x)) for x in breakpoints]

#Relate trees
relateTreeList = dendropy.TreeList.get(path=relateNewickFile,
                                     schema="newick",
                                     taxon_namespace=trueTreeList.taxon_namespace)
print("number of Relate trees: ", len(relateTreeList))

#read in Relate breakpoints
with open(relateBreakpointsFile, 'r') as breakpointsFile:
    relateBreakpoints = breakpointsFile.readline()
relateBreakpoints = relateBreakpoints[1:len(relateBreakpoints)-1].split(', ')
relateBreakpoints = [int(float(x)) for x in relateBreakpoints[1:]]

#convert rent trees file to Newick format
with open(rentTreesFile, 'r') as fin, open(rentNewickFile, 'w') as fout:
    tsv_reader = csv.reader(fin, delimiter="\t")
    positions = []
    for row in tsv_reader:
        positions.append(int(row[0]))
        fout.write(row[1] + ';\n')

estTreeList = dendropy.TreeList.get(
        path=rentNewickFile,
        schema="newick",
        taxon_namespace=trueTreeList.taxon_namespace)
print("number of RENT+ trees: ", len(estTreeList))

#import ARGweaver trees
argweaverTreeList = dendropy.TreeList.get(path=argweaverNewickFile, schema="newick")
print("number of ARGweaver trees: ", len(argweaverTreeList))                          

#convert taxon labels from 0-indexed to 1-indexed
for i in range((11), -1, -1):
    for tree in argweaverTreeList:
        node_to_change = tree.find_node_with_taxon_label(str(i))
        if(node_to_change!=None):
            node_to_change.taxon.label = str(i+1)

argweaverTreeList.migrate_taxon_namespace(trueTreeList.taxon_namespace, unify_taxa_by_label=True)
print("argweaver taxa:")
print(argweaverTreeList.poll_taxa())

with open(argweaverBreakpointsFile, 'r') as breakpointsFile:
    argweaverPositions = breakpointsFile.readline()
argweaverPositions = argweaverPositions[1:len(argweaverPositions)-1].split(', ')
argweaverPositions = [int(float(x)) for x in argweaverPositions]

combinedPositions = []
rentDistances = []
relateDistances = []

estPositionsIndex=0
        
xValues = (x for x in range(1,int(sequenceLength)) if x % 1 == 0)

newXValues = []

newXValues = []
trueBranchLengths = []
estBranchLengths = []
relateEstBranchLengths = []

#for each msprime position, get the TMRCA of the estimated RENT+ tree
for i in range(len(trueTreeList)):
    trueBranchLengths.append(trueTreeList[i].max_distance_from_root())


#for each RENT+ position, get the TMRCA of the estimated RENT+ tree
for i in range(len(estTreeList)):
    estBranchLengths.append(estTreeList[i].max_distance_from_root())

xValuesBranchLengths = []
relateXValues=[]
previousMsprimeIndex = 0
previousRelateIndex = 0
for i in xValues:
    msprimeIndex = bisect_left(breakpoints, i)
    relateIndex = bisect_left(relateBreakpoints, i)
    if(msprimeIndex < len(trueTreeList) and ((previousMsprimeIndex!=msprimeIndex and previousRelateIndex !=relateIndex) or i==1)
       and relateIndex < len(relateTreeList)):
        #compare Relate trees to true trees
        relateEstBranchLengths.append(relateTreeList[relateIndex].max_distance_from_root())
        xValuesBranchLengths.append(trueTreeList[msprimeIndex].max_distance_from_root())
        relateXValues.append(i)
        previousMsprimeIndex=msprimeIndex
        previousRelateIndex=relateIndex

plt.plot(breakpoints[1:], trueBranchLengths, label='true', marker='o')
#plt.plot(positions, estBranchLengths, label='RENT+', marker='o')
plt.plot(relateXValues, relateEstBranchLengths, label='Relate', marker='o')

#plt.plot(xValuesBranchLengths, relateEstBranchLengths, marker='o')
#plt.plot([0,1],[0,1], transform=plt.transAxes)

#print("true branch lengths: ", trueBranchLengths)
#print("est branch lengths: ", estBranchLengths)

#plt.plot(newXValues, trueBranchLengths)#, label='msprime')
#plt.plot(newXValues, estBranchLengths)#, label='RENT+')
plt.title("TMRCA of True Msprime and Relate Estimated Trees")
plt.xlabel("Position on genome")
plt.ylabel("TMRCA (in number of generations)")
plt.legend()
plt.tight_layout()
plt.savefig(plotFile)


