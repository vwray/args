# -*- coding: utf-8 -*-
"""
Generates ancestry using msprime based on a demography file. Converts the tree
sequence to VCF and to a Newick file and a list of positions.

Created on Wed Nov 15 17:55:33 2023

@author: Valerie Wray
"""

import msprime
import demes
import sys

numberOfSamples = int(sys.argv[1])
sequenceLength = sys.argv[2]
inputDemographyFile = sys.argv[3]
outputVCFFile = sys.argv[4]
outputNewickFile = sys.argv[5]
outputBreakpointsFile = sys.argv[6]
outputTreesFile = sys.argv[7]
recombRate = float(sys.argv[8])

model = demes.load(inputDemographyFile)
demo = msprime.Demography.from_demes(model)
ts = msprime.sim_ancestry(
                samples={'A': numberOfSamples, 'YRI': numberOfSamples},
                sequence_length= int(sequenceLength),
                recombination_rate= recombRate,
                random_seed=83,
                demography=demo,
                record_migrations=True,
            )
mts = msprime.sim_mutations(ts, rate=1.2e-8, random_seed=83)

#print(model)
#print(demo)

# output vcf file
with open(outputVCFFile, 'w') as vcf:
    mts.write_vcf(vcf, contig_id=1)
    
mts.dump(outputTreesFile)

nodeLabels = {}
for i in range(4*numberOfSamples):
    nodeLabels[i] = str(i+1)

with open(outputNewickFile, 'w') as newickFile:
    for tree in mts.trees():
        newickFile.write(tree.as_newick(precision=3, node_labels=nodeLabels))
        newickFile.write('\n')
    
with open(outputBreakpointsFile, 'w') as breakpointsFile:
    breakpointsFile.write(str(list(mts.breakpoints())))