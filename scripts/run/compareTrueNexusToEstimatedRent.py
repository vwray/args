# -*- coding: utf-8 -*-
"""

Created on Wed Nov 15 19:51:52 2023

A script to compare a Nexus true tree sequence (created from msprime) to 
estimated tree sequences (created by RENT+, Relate, and ARGweaver).

@author: valer
"""

import dendropy
import csv
from bisect import bisect_left
from compare_two_trees import compare_trees
#import httpimport
#with httpimport.remote_repo('https://github.com/ekmolloy/fastmulrfs/tree/master/python-tools'):
#    import compare_two_trees
from matplotlib import pyplot as plt

sequenceLength = '100_000'
numberOfSamples = 3



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
            



'''
taxaDictionary={}
for taxonA in trueTreeList.poll_taxa():
    for taxonB in estTreeList.poll_taxa():
        if int(taxonA.description().split("'")[-2])==int(taxonB.description().split("'")[-2])+1:
            taxaDictionary[taxonB]=taxonA
print(taxaDictionary) 
'''                              
                           
#dictionary = {}
#for i in range(4*numberOfSamples):
#    dictionary[str(i)] = str(i+1)
    
#print(dictionary)



argweaverTreeList.migrate_taxon_namespace(trueTreeList.taxon_namespace, unify_taxa_by_label=True)
print("ARGweaver taxon namespace: ")
print(argweaverTreeList.taxon_namespace)

print("true tree list poll taxa:")
print(trueTreeList.poll_taxa())
print("argweaver poll taxa:")
print(argweaverTreeList.poll_taxa())
#print(trueTreeList[0].poll_taxa()[0].)



    

with open("H://My Drive//Fall2023//research//argweaverOutput//" + sequenceLength + "//out_positions_5000iter_d20.txt", 'r') as breakpointsFile:
    argweaverPositions = breakpointsFile.readline()
#print(breakpoints)
argweaverPositions = argweaverPositions[1:len(argweaverPositions)-1].split(', ')
#print(breakpoints)
print("first breakpoint: ", argweaverPositions[0])
#breakpoints = list(map(int, breakpoints))
argweaverPositions = [int(float(x)) for x in argweaverPositions]



combinedPositions = []
rentDistances = []
relateDistances = []
argweaverDistances = []

rentFalsePositives = []
relateFalsePositives = []
argweaverFalsePositives = []

rentFalseNegatives = []
relateFalseNegatives = []
argweaverFalseNegatives = []

estPositionsIndex=0
        
xValues = (x for x in range(1,int(sequenceLength)) if x % 1 == 0)

newXValues = []

#for every position on the genome (or spaced apart by number specified above), 
#compare true msprime tree with estimated trees
previousIndex1 = 0
previousIndex3 = 0

for i in xValues:
    index1 = bisect_left(breakpoints, i)
    #index2 = bisect_left(positions, i)
    index3 = bisect_left(relateBreakpoints, i)
    if(index1 < len(trueTreeList) and ((previousIndex1!=index1 and previousIndex3 !=index3) or i==1) #and index2 < len(estTreeList) 
       and index3 < len(relateTreeList)):
        print(f"i: {i}")
        #[nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],estTreeList[index2])
        #rentDistances.append(rf)
        
        #use line below to compare Relate trees to true trees
        [nl, ei1, ei2, fn, fp, rf_2] = compare_trees(trueTreeList[index1-1],relateTreeList[index3])
        relateDistances.append(rf_2)
        relateFalsePositives.append(fp)
        relateFalseNegatives.append(fn)
        newXValues.append(i)
        previousIndex1=index1
        previousIndex3=index3
 #       if rf < 0.1:
 #           print(f"True tree at position {breakpoints[index1]} is similar to estimated tree at position {positions[index2]} with RF distance: {rf}")
 #       if rf > 0.9:
 #           print(f"True tree at position {breakpoints[index1]} is very different from estimated tree at position {positions[index2]} with RF distance: {rf}")

        






#for each RENT+ position, compare the true msprime tree with the estimated RENT+ tree
for i in range(len(positions)):
    index1 = bisect_left(breakpoints, positions[i])
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],estTreeList[i])
    rentDistances.append(rf)
    rentFalsePositives.append(fp)
    rentFalseNegatives.append(fn)
 
#for each Relate position, compare the true msprime tree with the estimated Relate tree
#for i in range(len(relateBreakpoints)):
#    index1 = bisect_left(breakpoints, relateBreakpoints[i])
#    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],relateTreeList[i])
#    relateDistances.append(rf)
#    relateFalsePositives.append(fp)

for i in range(len(argweaverPositions)):
    
    if(i==0):
        print("true tree:")
        print(trueTreeList[index1-1])
        print(trueTreeList[index1-1].taxon_namespace_scoped_copy)
        print("estimated tree:")
        print(argweaverTreeList[i])
    
    index1= bisect_left(breakpoints, argweaverPositions[i])
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[index1-1],argweaverTreeList[i])
    argweaverDistances.append(rf)
    argweaverFalsePositives.append(fp)
    argweaverFalseNegatives.append(fn)

#[nl, ei1, ei2, fn, fp, rf] = compare_trees(trueTreeList[0],argweaverTreeList[0])






trueDistances = []
#trueXValues = []
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
    
relateEstDistances = []
previousTree = relateTreeList[0]
for i in range(1, len(relateTreeList)):
    [nl, ei1, ei2, fn, fp, rf] = compare_trees(previousTree,relateTreeList[i])
    relateEstDistances.append(rf)
    previousTree = relateTreeList[i]
    
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
#argweaverEstXValues = argweaverPositions[1:]
print(len(trueXValues))
print(len(breakpoints))




#print(distances)

#plt.plot(newXValues, rentDistances, label='RENT+')
#plt.plot(newXValues, relateDistances, label='Relate')
'''
plt.figure(1)
ax1 = plt.subplot(211)
ax1.plot(positions, rentDistances, label='RENT+', marker='o')
ax1.plot(newXValues, relateDistances, label='Relate', marker='o')
ax1.plot(argweaverPositions, argweaverDistances, label='ARGweaver', marker='o')
plt.title("RF Distances and False Positives between True and Estimated Trees")
plt.tick_params(labelbottom = False, bottom = False)
plt.ylabel('RF Distance')
plt.legend()
ax2 = plt.subplot(212, sharex=ax1)
ax2.plot(positions, rentFalsePositives, label='RENT+', marker='o')
ax2.plot(newXValues, relateFalsePositives, label='Relate', marker='o')
ax2.plot(argweaverPositions, argweaverFalsePositives, label='ARGweaver', marker='o')
plt.ylabel('False Positives')
plt.xlabel('Position on genome')
plt.tight_layout()
'''

#Plot rf distances between true and est trees
#plt.plot(positions, rentDistances, label='RENT+', marker='o')
#plt.plot(newXValues, relateDistances, label='Relate', marker='o')
#plt.plot(argweaverPositions, argweaverDistances, label='ARGweaver', marker='o')

#Plot false positives
#plt.plot(positions, rentFalsePositives, label='RENT+', marker='o')
#plt.plot(newXValues, relateFalsePositives, label='Relate', marker='o')
#plt.plot(argweaverPositions, argweaverFalsePositives, label='ARGweaver', marker='o')

#Plot adjacent tree distances
#plt.plot(trueXValues, trueDistances, label='msprime', marker='o')
#plt.plot(estXValues, estDistances, label='RENT+', marker='o')
#plt.plot(relateEstXValues, relateEstDistances, label='Relate', marker='o')
#plt.plot(newXValues,relateDistances, label='Relate2', marker='o')
#plt.plot(argweaverEstXValues, argweaverEstDistances, label='ARGweaver', marker='o')

plt.figure(1)
ax1 = plt.subplot(211)
#ax1.plot(positions, rentDistances, label='RENT+', marker='o')
#ax1.plot(newXValues, relateDistances, label='Relate', marker='o')
#ax1.plot(argweaverPositions, argweaverDistances, label='ARGweaver', marker='o')
#plt.title("RF Distances and False Positives between True and Estimated Trees")
#plt.tick_params(labelbottom = False, bottom = False)
#plt.ylabel('RF Distance')
plt.legend()
ax2 = plt.subplot(211)
ax2.plot(positions, rentFalsePositives, label='RENT+', marker='o')
ax2.plot(newXValues, relateFalsePositives, label='Relate', marker='o')
ax2.plot(argweaverPositions, argweaverFalsePositives, label='ARGweaver', marker='o')
plt.tick_params(labelbottom = False, bottom = False)
plt.ylabel('False Positives')
plt.xlabel('Position on genome')
ax3 = plt.subplot(212, sharex=ax2)
ax3.plot(positions, rentFalseNegatives, label='RENT+', marker='o')
ax3.plot(newXValues, relateFalseNegatives, label='Relate', marker='o')
ax3.plot(argweaverPositions, argweaverFalseNegatives, label='ARGweaver', marker='o')
plt.ylabel('False Negatives')
plt.xlabel('Position on genome')
plt.tight_layout()




#plt.title("False Positives in Topology of Estimated Trees")
#plt.xlabel("Position on genome")
#plt.ylabel("RF Distance")
#plt.ylabel("Number of False Positives")
plt.legend()
plt.savefig('H://My Drive//Fall2023//research//outputFiles//plots//' + sequenceLength + '//rfdist_fp_combined_withArgweaver_fn.png')
