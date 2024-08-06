# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:42:20 2024

@author: valer
"""

import msprime
import tskit
import tsinfer
import tqdm

true_ts = tskit.load("msprimeTrees.trees")

with tsinfer.SampleData(
    path="tsinferInput.samples",
    sequence_length=true_ts.sequence_length,
    num_flush_threads=2
) as sim_sample_data:
    for var in tqdm.tqdm(true_ts.variants(), total=true_ts.num_sites):
        sim_sample_data.add_site(var.site.position, var.genotypes, var.alleles)



#just do this instead of above:
sample_data = tsinfer.SampleData.from_tree_sequence(true_ts, path="tsinferInput2.samples", num_flush_threads=2)


inferred_ts = tsinfer.infer(sim_sample_data).dump("tsinferOutput.trees")
inferred_ts = tsinfer.infer(sim_sample_data)
inferred_ts = inferred_ts.simplify(keep_unary=False)

import tsdate
dated_ts = tsdate.date(inferred_ts, Ne=10000, mutation_rate=1e-8)
dated_ts.dump("dated_ts.trees")