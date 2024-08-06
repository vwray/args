# -*- coding: utf-8 -*-
"""
A script to compare topologies of a true tree sequence created from msprime to 
estimated tree sequences created by RENT+, Relate, and ARGweaver.

Created on Wed Nov 15 19:51:52 2023

@author: Valerie Wray
"""

import dendropy
import csv
from bisect import bisect_left
from compare_two_trees import compare_trees
from matplotlib import pyplot as plt
import sys
import seaborn as sns
import tskit

#sns.set_theme()

numberOfSamples = sys.argv[1]
sequenceLength = sys.argv[2]
trueTreesNewickFile = sys.argv[3]
trueTreesBreakpointsFile = sys.argv[4]
relateNewickFile = sys.argv[5]
relateBreakpointsFile = sys.argv[6]
rentTreesFile = sys.argv[7]
rentNewickFile = sys.argv[8]
argweaverNewickFile = sys.argv[9]
argweaverBreakpointsFile = sys.argv[10]
plotFile = sys.argv[11]
trueTreesTreesFile = sys.argv[12]

#True msprime trees
trueTreesTskitTrees = tskit.load(trueTreesTreesFile)
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

rentDistances = []
relateDistances = []
argweaverDistances = []

rentFalsePositives = []
relateFalsePositives = []
argweaverFalsePositives = []

estPositionsIndex=0
        
xValues = (x for x in range(1,int(sequenceLength)) if x % 1 == 0)
newXValues = []

#for every position on the genome (or spaced apart by number specified above), 
#compare true msprime tree with estimated Relate trees
previousMsprimeIndex = 0
previousRelateIndex = 0
for i in xValues:
    msprimeIndex = bisect_left(breakpoints, i)
    relateIndex = bisect_left(relateBreakpoints, i)
    if(msprimeIndex < len(trueTreeList) and ((previousMsprimeIndex!=msprimeIndex and previousRelateIndex !=relateIndex) or i==1)
       and relateIndex < len(relateTreeList)):
        #compare Relate trees to true trees
        [nl, ei1, ei2, fn, fp, rf_2] = compare_trees(trueTreeList[msprimeIndex-1],relateTreeList[relateIndex])
        relateDistances.append(rf_2)
        relateFalsePositives.append(fp)
        newXValues.append(i)
        previousMsprimeIndex=msprimeIndex
        previousRelateIndex=relateIndex


#for each RENT+ position, compare the true msprime tree with the estimated RENT+ tree
for i in range(len(positions)):
    index1 = bisect_left(breakpoints, positions[i])
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],estTreeList[i])
    rentDistances.append(rf)
    rentFalsePositives.append(fp)
  
#for each ARGweaver position, compare the true msprime tree with the estimated ARGweaver tree
for i in range(len(argweaverPositions)):
    index1= bisect_left(breakpoints, argweaverPositions[i])
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],argweaverTreeList[i])
    argweaverDistances.append(rf)
    argweaverFalsePositives.append(fp)


trueDistances = []
'''
#compare adjacent true msprime trees
previousTree = trueTreeList[0]
for i in range(1, len(trueTreeList)):
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(previousTree,trueTreeList[i])
    trueDistances.append(rf)
    previousTree = trueTreeList[i]

#compare adjacent estimated RENT+ trees
estDistances = []
previousTree = estTreeList[0]
for i in range(1, len(estTreeList)):
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(previousTree,estTreeList[i])
    estDistances.append(rf)
    previousTree = estTreeList[i] 
    
#compare adjacent estimated Relate trees
relateEstDistances = []
previousTree = relateTreeList[0]
for i in range(1, len(relateTreeList)):
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(previousTree,relateTreeList[i])
    relateEstDistances.append(rf)
    previousTree = relateTreeList[i]
    
#compare adjacent estimated ARGweaver trees
argweaverEstDistances = []
previousTree = argweaverTreeList[0]
for i in range(1, len(argweaverTreeList)):
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(previousTree,argweaverTreeList[i])
    argweaverEstDistances.append(rf)
    previousTree = argweaverTreeList[i]
'''
trueXValues = breakpoints[1:len(breakpoints)-1]
estXValues = positions[:-1]
relateEstXValues = relateBreakpoints[:-1]
argweaverEstXValues = argweaverPositions[1:]

plt.figure(1)
ax1 = plt.subplot(211)
#Plot rf distances between true and est trees
ax1.plot(positions, rentDistances, label='RENT+')
ax1.plot(newXValues, relateDistances, label='Relate')
#ax1.plot(argweaverPositions, argweaverDistances, label='ARGweaver')
plt.title("RF Distances and False Positives between True and Estimated Trees")
plt.tick_params(labelbottom = False, bottom = False)
plt.ylabel('RF Distance')
plt.legend()
ax2 = plt.subplot(212, sharex=ax1)
#Plot false positives
ax2.plot(positions, rentFalsePositives, label='RENT+')
ax2.plot(newXValues, relateFalsePositives, label='Relate')
#ax2.plot(argweaverPositions, argweaverFalsePositives, label='ARGweaver')


#get migrating tracts (introgressed segments)
neanderthal_id = [p.id for p in trueTreesTskitTrees.populations() if p.metadata['name']=='B'][0]
migrating_tracts = []
# Get all tracts that migrated into the neanderthal population
for migration in trueTreesTskitTrees.migrations():
    if migration.dest == neanderthal_id:
        migrating_tracts.append((migration.left, migration.right))
        
for tract in migrating_tracts:
    print("begin tract: ", tract[0])
    print("end tract: ", tract[1])
    #ax1.scatter(tract[0], 0.5)
    #ax1.scatter(tract[1], 0.7)
    ax1.fill_between(tract, [1,1], facecolor='green', alpha=.6)
    ax2.fill_between(tract, [9,9], facecolor='green', alpha=.6)


#Plot adjacent tree distances
#plt.plot(trueXValues, trueDistances, label='msprime', marker='o')
#plt.plot(estXValues, estDistances, label='RENT+', marker='o')
#plt.plot(relateEstXValues, relateEstDistances, label='Relate', marker='o')
#plt.plot(newXValues,relateDistances, label='Relate2', marker='o')
#plt.plot(argweaverEstXValues, argweaverEstDistances, label='ARGweaver', marker='o')

#plt.title("Distances between Adjacent True and Estimated Trees")
#plt.title("False Positives in Topology of Estimated Trees")
plt.xlabel("Position on genome")
#plt.ylabel("RF Distance")
plt.ylabel("Number of False Positives")
plt.legend()
plt.savefig(plotFile)
