# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 17:27:59 2023

@author: valer
"""

#import argutils
import tskit
#import tszip

#with open("H:\My Drive\Fall2023\research\argweaverOutput\out.0.arg") as f:
#    ts = argutils.convert_argweaver(f)
#    ts.dump("H:\My Drive\Fall2023\research\argweaverOutput\argweaverTrees.trees")
    
msprimeTrees = tskit.load("C:/Users/valer/Downloads/vcfToSitesTest/2024-2-13/sim_l250kb_0.trees")

#msprimeTrees = tszip.decompress("C:/Users/valer/Downloads/vcfToSitesTest/2024-2-13/msprimeTrees.trees")

#msprimeTrees.first()

for tree in msprimeTrees.trees():
    print(tree.draw_text())
    
print("first true tree:")
print(msprimeTrees.first().draw_text())