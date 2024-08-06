# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 20:35:18 2023

A script to plot ARGweaver stats to monitor convergence.

@author: valer
"""

import csv
from matplotlib import pyplot as plt

sequenceLength = '100_000'
numberOfSamples = 3

inputFile = "H:\\My Drive\\Fall2023\\research\\argweaverOutput\\" + sequenceLength + "\\out.stats"

with open(inputFile, 'r') as fin:
    tsv_reader = csv.reader(fin, delimiter="\t")
    next(tsv_reader)
    next(tsv_reader)
    iteration=[]
    prior=[]
    likelihood=[]
    joint=[]
    recombs=[]
    noncompats=[]
    arglen=[]
    
    for row in tsv_reader:
        iteration.append(int(row[1]))
        prior.append(float(row[2]))
        likelihood.append(float(row[3]))
        joint.append(float(row[4]))
        recombs.append(int(row[5]))
        noncompats.append(int(row[6]))
        arglen.append(float(row[7]))
        #tree = tsconvert.from_newick(newickString, min_edge_length=0)
        #print(tree)
    
    #print(iteration)
    print(prior)

    import numpy as np
    
    t = np.arange(0, 2001, 1)
    '''
    ax1 = plt.subplot(211)
    ax1.plot(t, np.sin(2*np.pi*t))
    
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.plot(t, np.sin(4*np.pi*t))
    
    plt.show()
    
    
    '''
    plt.figure(1)
    ax1 = plt.subplot(321)
    ax1.plot(iteration,prior)
    plt.ylabel('Prior')
    plt.tick_params(labelbottom = False, bottom = False)
    ax2 = plt.subplot(322, sharex=ax1)
    ax2.plot(iteration,likelihood)
    plt.ylabel('Likelihood')
    plt.tick_params(labelbottom = False, bottom = False)
    ax3 = plt.subplot(323, sharex=ax1)
    ax3.plot(iteration,joint)
    plt.ylabel('Joint')
    plt.tick_params(labelbottom = False, bottom = False)
    ax4 = plt.subplot(324, sharex=ax1)
    ax4.plot(iteration,recombs)
    plt.ylabel('Recombs')
    plt.tick_params(labelbottom = False, bottom = False)
    ax5 = plt.subplot(325, sharex=ax1)
    ax5.plot(iteration,noncompats)
    plt.ylabel('Noncompats')
    plt.xlabel('MCMC Iteration')
    ax6 = plt.subplot(326, sharex=ax1)
    ax6.plot(iteration,arglen)
    plt.ylabel('Arglen')
    plt.xlabel('MCMC Iteration')
    plt.tight_layout()
    plt.savefig('H://My Drive//Fall2023//research//outputFiles//plots//' + sequenceLength + '//argweaverConvergence_5000iter.png')
