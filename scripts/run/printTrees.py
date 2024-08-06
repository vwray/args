# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 09:19:16 2023

Print out the full ARG in an interactive HTML file.

@author: valer
"""

#import tsconvert

#with open("H://My Drive//Fall2023//research//outputFiles//outputVCFFile_3samples_10_000SequenceLen_nodeLabels.vcf", 'r') as vcf:
#    mts = vcf.

import numpy as np
import msprime

parameters = {
    "samples": 3, # Three diploid individuals == six sample genomes
    "sequence_length": 1e4,
    "recombination_rate": 1e-7,
    "population_size": 1e3,
    "random_seed": 333,
}

ts_arg = msprime.sim_ancestry(**parameters, record_full_arg=True, discrete_genome=False)
# NB: the strict Hudson ARG needs unique crossover positions (i.e. a continuous genome)

print('Simulated a "full ARG" under the Hudson model:')
print(
    f" ARG stored in a tree sequence with {ts_arg.num_nodes} nodes and"
    f" {ts_arg.num_edges} edges (creating {ts_arg.num_trees} local trees)"
)

mu = 1e-7
ts_arg = msprime.sim_mutations(ts_arg, rate=mu, random_seed=888)
print("     Sample node:  " + "   ".join(str(u) for u in ts_arg.samples()))
for v in ts_arg.variants():
    print(f"Variable site {v.site.id}:", np.array(v.alleles)[v.genotypes])
    
import tskit_arg_visualizer
d3arg = tskit_arg_visualizer.D3ARG.from_ts(ts=ts_arg)
w, h = 450, 300  # width and height
d3arg.draw(w, h, edge_type="ortho", sample_order=[0, 2, 1, 3, 5, 4])