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
#from compare_two_trees import compare_trees
from matplotlib import pyplot as plt
import seaborn as sns

sns.set_theme()

sequenceLength = '100_000'

#True msprime trees
trueTreeList = dendropy.TreeList.get(path="H://My Drive//Fall2023//research//outputFiles//trueTrees//" + sequenceLength + "//newick.txt",
                                     schema="newick")
print("number of true trees: ", len(trueTreeList))
#print(trueTreeList)

#read in breakpoints
with open("H://My Drive//Fall2023//research//outputFiles//trueTrees//" + sequenceLength + "//breakpoints.txt", 'r') as breakpointsFile:
    breakpoints = breakpointsFile.readline()
#print(breakpoints)
breakpoints = breakpoints[1:len(breakpoints)-1].split(', ')
#print(breakpoints)
print("first breakpoint: ", breakpoints[0])
#breakpoints = list(map(int, breakpoints))
breakpoints = [int(float(x)) for x in breakpoints]

#Relate trees
relateTreeList = dendropy.TreeList.get(path="H://My Drive//Fall2023//research//relateScripts//" + sequenceLength + "//relateNewick.txt",
                                     schema="newick",
                                     taxon_namespace=trueTreeList.taxon_namespace)
print("number of Relate trees: ", len(relateTreeList))
#print(trueTreeList)

#read in breakpoints
with open("H://My Drive//Fall2023//research//relateScripts//" + sequenceLength + "//relateBreakpoints.txt", 'r') as breakpointsFile:
    relateBreakpoints = breakpointsFile.readline()
print(relateBreakpoints)
relateBreakpoints = relateBreakpoints[1:len(relateBreakpoints)-1].split(', ')
print("first Relate breakpoint: ", relateBreakpoints[0])
print(relateBreakpoints)
#breakpoints = list(map(int, breakpoints))
relateBreakpoints = [int(float(x)) for x in relateBreakpoints[1:]]
print("number of Relate breakpoints: ", len(relateBreakpoints))

#fix rent tree format
rentTreeFile="H://My Drive//Fall2023//research//rentOutput//" + sequenceLength + "//rentInputFile.txt.trees"
modifiedRentTreeFile = "H://My Drive//Fall2023//research//rentOutput//" + sequenceLength + "//rentNewick.txt"
with open(rentTreeFile, 'r') as fin, open(modifiedRentTreeFile, 'w') as fout:
    tsv_reader = csv.reader(fin, delimiter="\t")
    positions = []
    
    for row in tsv_reader:
        positions.append(int(row[0]))
        #newRow = row.
        fout.write(row[1] + ';\n')
#print(positions)
#print("first position: ", positions[0])
estTreeList = dendropy.TreeList.get(
        path=modifiedRentTreeFile,
        schema="newick",
        taxon_namespace=trueTreeList.taxon_namespace)

print("number of RENT+ trees: ", len(estTreeList))
print("number of RENT+ positions: ", len(positions))




'''
#import ARGweaver trees
argweaverTreeList = dendropy.TreeList.get(path="H://My Drive//Fall2023//research//argweaverOutput//" + sequenceLength + "//out_processedNewickTrees_5000iter_d20.txt",
                                     schema="newick")
                                     #taxon_namespace=trueTreeList.taxon_namespace)
print("first ARGweaver tree before migrating namespace")                          
print(argweaverTreeList[0])   


for i in range((11), -1, -1):
    for tree in argweaverTreeList:
        #tree.find_node_with_taxon_label("")
        node_to_change = tree.find_node_with_taxon_label(str(i))
        if(node_to_change!=None):
            print(f"node to change: {node_to_change}. Changing to: {i+1}")
            node_to_change.taxon.label = str(i+1)
            
argweaverTreeList.migrate_taxon_namespace(trueTreeList.taxon_namespace, unify_taxa_by_label=True)
print("ARGweaver taxon namespace: ")
print(argweaverTreeList.taxon_namespace)

print("true tree list poll taxa:")
print(trueTreeList.poll_taxa())
print("argweaver poll taxa:")
print(argweaverTreeList.poll_taxa())

with open("H://My Drive//Fall2023//research//argweaverOutput//" + sequenceLength + "//out_positions_5000iter_d20.txt", 'r') as breakpointsFile:
    argweaverPositions = breakpointsFile.readline()
#print(breakpoints)
argweaverPositions = argweaverPositions[1:len(argweaverPositions)-1].split(', ')
#print(breakpoints)
print("first breakpoint: ", argweaverPositions[0])
#breakpoints = list(map(int, breakpoints))
argweaverPositions = [int(float(x)) for x in argweaverPositions]
'''




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
    #index1 = bisect_left(breakpoints, positions[i])
    #[nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],estTreeList[i])
    #rentDistances.append(rf)
    trueBranchLengths.append(trueTreeList[i].max_distance_from_root())


#for each RENT+ position, get the TMRCA of the estimated RENT+ tree
for i in range(1, len(estTreeList)):
    #index1 = bisect_left(breakpoints, positions[i])
    #[nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],estTreeList[i])
    #rentDistances.append(rf)
    estBranchLengths.append(estTreeList[i].max_distance_from_root()*20000)
    
print(estBranchLengths)

xValuesBranchLengths = []
'''
#for each Relate position, get the TMRCA of the estimated Relate tree
for i in range(len(relateTreeList)):
    index1 = bisect_left(breakpoints, relateBreakpoints[i])
    #[nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],relateTreeList[i])
    #relateDistances.append(rf)
    relateEstBranchLengths.append(relateTreeList[i].max_distance_from_root())
    xValuesBranchLengths.append(trueTreeList[index1].max_distance_from_root())
'''
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

#for i in range(len(argweaverTreeList)):
#    argweaverEstBranchLengths.append(argweaverTreeList[i].max_distance_from_root())


with open("C://Users//valer//Downloads//vcfToSitesTest//argweaverTMRCA.txt", 'r') as file:
    tsv_reader = csv.reader(file, delimiter="\t")
    argweaverPositions = []
    argweaverTmrcaAvg = []
    argweaverTmrcaLow = []
    argweaverTmrcaHigh = []
    
    for row in tsv_reader:
        argweaverPositions.append(int(row[1]))
        argweaverTmrcaAvg.append(float(row[3]))
        argweaverTmrcaLow.append(float(row[4]))
        argweaverTmrcaHigh.append(float(row[5]))


'''

for i in xValues:
    index1 = bisect_left(breakpoints, i)
    index2 = bisect_left(positions, i)
    #print(f"Index 1: {index1}, Index 2: {index2}")
    if(index1 < len(trueTreeList) and index2 < len(estTreeList)):
        [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],estTreeList[index2])
        rentDistances.append(rf)
        newXValues.append(i)
        
        trueBranchLengths.append(trueTreeList[index1-1].max_distance_from_root())
        relateEstBranchLengths.append(relateTreeList[index2].max_distance_from_root())
    else:
        break
    
'''
        
  
#plt.style.use('seaborn-darkgrid')
plt.plot(breakpoints[1:], trueBranchLengths, label='true')#, marker='o')
plt.plot(positions[1:], estBranchLengths, label='RENT+')#, marker='o')
plt.plot(relateXValues, relateEstBranchLengths, label='Relate')#, marker='o')
plt.plot(argweaverPositions, argweaverTmrcaAvg, label='ARGweaver')
#plt.plot(argweaverPositions, argweaverTmrcaLow, label='ARGweaver low')
#plt.plot(argweaverPositions, argweaverTmrcaHigh, label='ARGweaver high')

#plt.plot(xValuesBranchLengths, relateEstBranchLengths, marker='o')
#plt.plot([0,1],[0,1], transform=plt.transAxes)

#print("true branch lengths: ", trueBranchLengths)
#print("est branch lengths: ", estBranchLengths)

#plt.plot(newXValues, trueBranchLengths)#, label='msprime')
#plt.plot(newXValues, estBranchLengths)#, label='RENT+')
plt.title("TMRCA of True Msprime and Estimated Trees")
plt.xlabel("Position on genome")
plt.ylabel("TMRCA (in number of generations)")
plt.legend()
plt.tight_layout()
plt.savefig('H://My Drive//Fall2023//research//outputFiles//plots//' + sequenceLength + '//TMRCA_true_est_withRentandARGweaver.png')


