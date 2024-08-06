# -*- coding: utf-8 -*-
"""
A script to compare topologies of a true tree sequence created from msprime to 
estimated tree sequences created by RENT+, Relate, and ARGweaver, with both
ARGweaver consensus trees and ARGweaver trees from last MCMC iteration ARG.

Created on Fri Feb 16 10:38:21 2024

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

sns.set_theme()

numberOfSamples = sys.argv[1]
sequenceLength = sys.argv[2]
trueTreesNewickFile = sys.argv[3]
trueTreesBreakpointsFile = sys.argv[4]
relateNewickFile = sys.argv[5]
relateBreakpointsFile = sys.argv[6]
rentTreesFile = sys.argv[7]
rentNewickFile = sys.argv[8]
argweaverConsensusNewickFile = sys.argv[9]
argweaverConsensusBreakpointsFile = sys.argv[10]
argweaverLastARGTreesFile = sys.argv[11]
argweaverLastARGNewickFile = sys.argv[12]
plotFile = sys.argv[13]
trueTreesTreesFile = sys.argv[14]
averageErrorFile = sys.argv[15]

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

#import ARGweaver consensus trees
argweaverConsensusTreeList = dendropy.TreeList.get(path=argweaverConsensusNewickFile, schema="newick")
print("number of ARGweaver consensus trees: ", len(argweaverConsensusTreeList))                          

#convert taxon labels from 0-indexed to 1-indexed
for i in range((11), -1, -1):
    for tree in argweaverConsensusTreeList:
        node_to_change = tree.find_node_with_taxon_label(str(i))
        if(node_to_change!=None):
            node_to_change.taxon.label = str(i+1)

argweaverConsensusTreeList.migrate_taxon_namespace(trueTreeList.taxon_namespace, unify_taxa_by_label=True)
print("argweaver taxa:")
print(argweaverConsensusTreeList.poll_taxa())

with open(argweaverConsensusBreakpointsFile, 'r') as breakpointsFile:
    argweaverConsensusPositions = breakpointsFile.readline()
argweaverConsensusPositions = argweaverConsensusPositions[1:len(argweaverConsensusPositions)-1].split(', ')
argweaverConsensusPositions = [int(float(x)) for x in argweaverConsensusPositions]

# import ARGweaver last ARG trees
nodeLabels = {}
for i in range(4*int(numberOfSamples)):
    nodeLabels[i] = str(i+1)
argweaverLastARGTrees = tskit.load(argweaverLastARGTreesFile)
with open(argweaverLastARGNewickFile, 'w') as newickFile:
    for tree in argweaverLastARGTrees.trees():
        newickFile.write(tree.as_newick(precision=3, node_labels=nodeLabels))
        newickFile.write('\n')

print("ARGweaver last ARG first tree:")
print(argweaverLastARGTrees.first().draw_text())   

trueTreesTskitTrees = tskit.load(trueTreesTreesFile)

print("True trees first tree:")
print(trueTreesTskitTrees.first().draw_text()) 

argweaverLastARGTreeList = dendropy.TreeList.get(path=argweaverLastARGNewickFile, schema="newick")
argweaverLastARGTreeList.migrate_taxon_namespace(trueTreeList.taxon_namespace, unify_taxa_by_label=True)
print("number of argweaver last ARG trees: ", len(argweaverLastARGTreeList))
    
#with open(argweaverLastARGBreakpointsFile, 'w') as breakpointsFile:
#    breakpointsFile.write(str(list(argweaverLastARGTrees.breakpoints())))
argweaverLastARGPositions = list(argweaverLastARGTrees.breakpoints())
print("number of argweaver last ARG positions: ", len(argweaverLastARGPositions))

rentDistances = []
relateDistances = []
argweaverConsensusDistances = []
argweaverLastARGDistances = []

rentFalsePositives = []
relateFalsePositives = []
argweaverConsensusFalsePositives = []
argweaverLastARGFalsePositives = []

rentFalseNegatives = []
relateFalseNegatives = []
argweaverConsensusFalseNegatives = []
argweaverLastARGFalseNegatives = []

rentKCDistances = []
relateKCDistances = []
argweaverConsensusKCDistances = []
argweaverLastARGKCDistances = []

estPositionsIndex=0
        
xValues = []#(x for x in range(1,int(sequenceLength)) if x % 1 == 0)
for i in range(1,int(sequenceLength)):
    if i % 1 ==0:
        xValues.append(i)
newXValues = []
newARGweaverConsensusXValues = []
newARGweaverLastARGXValues = []

#for every position on the genome (or spaced apart by number specified above), 
#compare true msprime tree with estimated Relate trees
previousMsprimeIndex = 0
previousRelateIndex = 0
previousARGweaverConsensusIndex = 0
previousARGweaverLastARGIndex = 0
for i in xValues:
    msprimeIndex = bisect_left(breakpoints, i)
    relateIndex = bisect_left(relateBreakpoints, i)
    if(msprimeIndex < len(trueTreeList) and ((previousMsprimeIndex!=msprimeIndex or previousRelateIndex !=relateIndex) or i==1)
       and relateIndex < len(relateTreeList)):
        #compare Relate trees to true trees
        [nl, ei1, ei2, fn, fp, rf_2] = compare_trees(trueTreeList[msprimeIndex-1],relateTreeList[relateIndex])
        #kcDist = trueTreeList[msprimeIndex-1].kc_distance(relateTreeList[relateIndex])
        relateDistances.append(rf_2)
        relateFalsePositives.append(fp)
        relateFalseNegatives.append(fn)
        #relateKCDistances.append(kcDist)
        newXValues.append(i)
        previousMsprimeIndex=msprimeIndex
        previousRelateIndex=relateIndex
   
for j in xValues:
    msprimeIndex = bisect_left(breakpoints, j)
    argweaverConsensusIndex = bisect_left(argweaverConsensusPositions, j)
    if(msprimeIndex < len(trueTreeList) and ((previousMsprimeIndex!=msprimeIndex or previousARGweaverConsensusIndex !=argweaverConsensusIndex) or j==1)
       and argweaverConsensusIndex < len(argweaverConsensusTreeList)):
        #compare ARGweaver consensus trees to true trees
        [nl, ei1, ei2, fn, fp, rf_2] = compare_trees(trueTreeList[msprimeIndex-1],argweaverConsensusTreeList[argweaverConsensusIndex])
        #kcDist = trueTreeList[msprimeIndex-1].kc_distance(argweaverConsensusTreeList[argweaverConsensusIndex])
        argweaverConsensusDistances.append(rf_2)
        argweaverConsensusFalsePositives.append(fp)
        argweaverConsensusFalseNegatives.append(fn)
        #argweaverConsensusKCDistances.append(kcDist)
        newARGweaverConsensusXValues.append(j)
        previousMsprimeIndex=msprimeIndex
        previousARGweaverConsensusIndex=argweaverConsensusIndex
        
for i in xValues:
    msprimeIndex = bisect_left(breakpoints, i)
    argweaverLastARGIndex = bisect_left(argweaverLastARGPositions, i)
    if(msprimeIndex < len(trueTreeList) and ((previousMsprimeIndex!=msprimeIndex or previousARGweaverLastARGIndex !=argweaverLastARGIndex) or i==1)
       and argweaverLastARGIndex < len(argweaverLastARGTreeList)):
        #compare ARGweaver consensus trees to true trees
        [nl, ei1, ei2, fn, fp, rf_2] = compare_trees(trueTreeList[msprimeIndex-1],argweaverLastARGTreeList[argweaverLastARGIndex])
        argweaverLastARGDistances.append(rf_2)
        argweaverLastARGFalsePositives.append(fp)
        argweaverLastARGFalseNegatives.append(fn)
        newARGweaverLastARGXValues.append(i)
        previousMsprimeIndex=msprimeIndex
        previousARGweaverLastARGIndex=argweaverLastARGIndex


#take RF distance only at SNPs
relateSNPDists = []
argweaverConsensusSNPDists = []
argweaverLastARGSNPDists = []
#for each RENT+ position, compare the true msprime tree with the estimated RENT+ tree
for i in range(len(positions)):
    index1 = bisect_left(breakpoints, positions[i])
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],estTreeList[i])
    rentDistances.append(rf)
    rentFalsePositives.append(fp)
    rentFalseNegatives.append(fn)
    
    relateIndex = bisect_left(relateBreakpoints, positions[i])
    if (relateIndex < len(relateTreeList)):
        [nl, ei1, ei2, fn, fp, rf_2] = compare_trees(trueTreeList[index1-1],relateTreeList[relateIndex])
        relateSNPDists.append(rf_2)
        
    argweaverConsensusIndex = bisect_left(argweaverConsensusPositions, positions[i])
    if (argweaverConsensusIndex < len(argweaverConsensusTreeList)):
        [nl, ei1, ei2, fn, fp, rf_3] = compare_trees(trueTreeList[index1-1],argweaverConsensusTreeList[argweaverConsensusIndex])
        argweaverConsensusSNPDists.append(rf_3)
    
    argweaverLastARGIndex = bisect_left(argweaverLastARGPositions, positions[i])
    if (argweaverLastARGIndex < len(argweaverLastARGTreeList)):
        [nl, ei1, ei2, fn, fp, rf_4] = compare_trees(trueTreeList[index1-1],argweaverLastARGTreeList[argweaverLastARGIndex])
        argweaverLastARGSNPDists.append(rf_4)
    
  
'''
#for each ARGweaver consensus position, compare the true msprime tree with the estimated ARGweaver tree
for i in range(len(argweaverConsensusPositions)):
    index1= bisect_left(breakpoints, argweaverConsensusPositions[i])
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],argweaverConsensusTreeList[i])
    argweaverConsensusDistances.append(rf)
    argweaverConsensusFalsePositives.append(fp)
  
    
print("first dendropy trees:")
print("first true tree ")
trueTreeList[0].print_plot()
print("first ARGweaver last ARG tree ")
argweaverLastARGTreeList[0].print_plot()
    
#for each ARGweaver last ARG position, compare the true msprime tree with the estimated ARGweaver tree
for i in range(1,len(argweaverLastARGPositions)):
    index1= bisect_left(breakpoints, argweaverLastARGPositions[i])
    if i==0:
        print("comparing true tree ", index1)
        trueTreeList[index1-1].print_plot()
        print("to ARGweaver last ARG tree ", i)
        argweaverLastARGTreeList[i].print_plot()
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],argweaverLastARGTreeList[i-1])
    argweaverLastARGDistances.append(rf)
    argweaverLastARGFalsePositives.append(fp)
'''  



#get migrating tracts (introgressed segments)
neanderthal_id = [p.id for p in trueTreesTskitTrees.populations() if p.metadata['name']=='B'][0]
migrating_tracts = []
# Get all tracts that migrated into the neanderthal population
for migration in trueTreesTskitTrees.migrations():
    if migration.dest == neanderthal_id:
        migrating_tracts.append((migration.left, migration.right))


print("first dendropy trees:")
print("first true tree ")
trueTreeList[0].print_plot()
print("first ARGweaver last ARG tree ")
argweaverLastARGTreeList[0].print_plot()
print("first ARGweaver consensus tree ")
argweaverConsensusTreeList[0].print_plot()


#calculate average error of each method    
with open(averageErrorFile, 'w') as fout:
    fout.write(',Average RF Distance Over All Data Points,Average RF Distance at SNPs\n')
    fout.write('RENT+,' + str(sum(rentDistances)/len(rentDistances)) + ',' + str(sum(rentDistances)/len(rentDistances)) + '\n')
    fout.write('Relate,' + str(sum(relateDistances)/len(relateDistances)) + ',' + str(sum(relateSNPDists)/len(relateSNPDists)) + '\n')
    fout.write('ARGweaver (last ARG),' + str(sum(argweaverLastARGDistances)/len(argweaverLastARGDistances)) + ',' + str(sum(argweaverLastARGSNPDists)/len(argweaverLastARGSNPDists)) + '\n')
    fout.write('ARGweaver (consensus trees),' + str(sum(argweaverConsensusDistances)/len(argweaverConsensusDistances)) + ',' + str(sum(argweaverConsensusSNPDists)/len(argweaverConsensusSNPDists)) + '\n')
    

trueDistances = []

plt.figure(1)
ax1 = plt.subplot(311)
#Plot rf distances between true and est trees
ax1.plot(positions, rentDistances, label='RENT+')
ax1.plot(newXValues, relateDistances, label='Relate')
ax1.plot(newARGweaverConsensusXValues, argweaverConsensusDistances, label='ARGweaver consensus')
ax1.plot(newARGweaverLastARGXValues, argweaverLastARGDistances, label='ARGweaver last ARG')
plt.title("Tree Topology Comparison")
plt.tick_params(labelbottom = False, bottom = False)
plt.ylabel('RF Distance')
plt.legend()
ax2 = plt.subplot(312, sharex=ax1)
#Plot false positives
ax2.plot(positions, rentFalsePositives, label='RENT+')
ax2.plot(newXValues, relateFalsePositives, label='Relate')
ax2.plot(newARGweaverConsensusXValues, argweaverConsensusFalsePositives, label='ARGweaver consensus')
ax2.plot(newARGweaverLastARGXValues, argweaverLastARGFalsePositives, label='ARGweaver last ARG')

plt.ylabel("FP")

plt.tick_params(labelbottom = False, bottom = False)

ax3 = plt.subplot(313, sharex=ax1)
#Plot false positives
ax3.plot(positions, rentFalseNegatives, label='RENT+')
ax3.plot(newXValues, relateFalseNegatives, label='Relate')
ax3.plot(newARGweaverConsensusXValues, argweaverConsensusFalseNegatives, label='ARGweaver consensus')
ax3.plot(newARGweaverLastARGXValues, argweaverLastARGFalseNegatives, label='ARGweaver last ARG')

for tract in migrating_tracts:
    print("begin tract: ", tract[0])
    print("end tract: ", tract[1])
    #ax1.scatter(tract[0], 0.5)
    #ax1.scatter(tract[1], 0.7)
    ax1.fill_between(tract, [1,1], facecolor='green', alpha=.1)
    ax2.fill_between(tract, [1,1], facecolor='green', alpha=.1)
    ax3.fill_between(tract, [1,1], facecolor='green', alpha=.1)


'''
ax1 = plt.subplot(211)
#Plot rf distances between true and est trees
ax1.plot(positions, rentDistances, label='RENT+')
ax1.plot(newXValues, relateDistances, label='Relate')
ax1.plot(newARGweaverConsensusXValues, argweaverConsensusDistances, label='ARGweaver consensus')
ax1.plot(newARGweaverLastARGXValues, argweaverLastARGDistances, label='ARGweaver last ARG')
plt.title("RF Distances and False Positives between True and Estimated Trees")
plt.tick_params(labelbottom = False, bottom = False)
plt.ylabel('RF Distance')
plt.legend()
ax2 = plt.subplot(212, sharex=ax1)
#Plot false positives
ax2.plot(positions, rentFalsePositives, label='RENT+')
ax2.plot(newXValues, relateFalsePositives, label='Relate')
ax2.plot(newARGweaverConsensusXValues, argweaverConsensusFalsePositives, label='ARGweaver consensus')
ax2.plot(newARGweaverLastARGXValues, argweaverLastARGFalsePositives, label='ARGweaver last ARG')
'''



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
plt.ylabel("FN")
plt.savefig(plotFile)
