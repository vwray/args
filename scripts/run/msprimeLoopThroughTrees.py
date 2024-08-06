# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 17:55:33 2023

@author: Valerie Wray and Yulin Zhang
"""

import msprime
import demes

numberOfSamples = 3
sequenceLength = '100_000'

model = demes.load("C:\\Users\\valer\\Downloads\\simple_ooa_neanderthal5r19_simp.yaml")
demo = msprime.Demography.from_demes(model)
ts = msprime.sim_ancestry(
                samples={'A': numberOfSamples, 'YRI': numberOfSamples},
                sequence_length= int(sequenceLength), #10mbp #params.sequence_length,
                recombination_rate= 1e-8, #10^-8#params.recombination_rate,
                random_seed=83, #seed,
                demography=demo,
                record_migrations=True,
            )

mts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed= 83) #seed) #1.2e-8

print(model)
print(demo)

#with open("H://My Drive//Fall2023//research//outputFiles//demography.txt", 'w') as demoFile:
#    demoFile.write(demo.__str__())

print(ts.first().draw_text())
#print(mts)
#SVG(ts.draw_svg())

#tszip.compress(ts, "H://My Drive//Fall2023//research//outputFiles//simulation-trees_3samples_10_000SequenceLen.tsz")

# output vcf file
with open("H://My Drive//Fall2023//research//msprimeOutput//" + sequenceLength + "//outputVCFFile.vcf", 'w') as vcf:
    mts.write_vcf(vcf, contig_id=1)

#print(list(ts.breakpoints()))
#print(list(mts.breakpoints()))
'''
previousTree = ts.first()
for x in ts.breakpoints():
    print("breakpoint: ", x)

#print("size of breakpoints: ", len(ts.breakpoints))
for i in range(len(ts.mutations_site)):
    print("site: ", ts.mutations_site[i])
for tree in ts.trees():
    print(f"Tree {tree.index} covers {tree.interval}")
    #print(f"Tree {tree.index} at position {tree.sites_position} has a mutation.")
    print(tree.as_newick(precision = 0))
    #print("distance: ", tree.kc_distance(previousTree))
    previousTree = tree
    if tree.index >= 4:
        print("...")
        break
print(f"Tree {ts.last().index} covers {ts.last().interval}")
print(ts.first().draw_text())
print(ts.at_index(1).draw_text())
display(ts.draw_svg(x_lim=(0, 50)))
#help(ts.draw_svg)
'''
nodeLabels = {}
for i in range(4*numberOfSamples):
    nodeLabels[i] = str(i)
    
#print(nodeLabels)
print(mts.first().draw_text())
#print(mts.first().as_newick())
print(mts.first().as_newick(node_labels=nodeLabels))#, include_branch_lengths=False))
#help(mts.first().draw_text)

'''
with open("H://My Drive//Fall2023//research//msprimeOutput//" + sequenceLength + "//trueTreesNewick.txt", 'w') as newickFile:
    for tree in mts.trees():
        #print(tree.as_newick(precision=3, node_labels=nodeLabels))
        newickFile.write(tree.as_newick(precision=3, node_labels=nodeLabels))
        newickFile.write('\n')
    
with open("H://My Drive//Fall2023//research//msprimeOutput//" + sequenceLength + "//trueTreesNewickBreakpoints.txt", 'w') as breakpointsFile:
    breakpointsFile.write(str(list(mts.breakpoints())))
'''