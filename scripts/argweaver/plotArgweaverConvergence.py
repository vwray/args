# -*- coding: utf-8 -*-
"""
A script to plot ARGweaver statistics to monitor convergence.

Created on Sat Dec 16 20:35:18 2023

@author: Valerie Wray
"""

import csv
from matplotlib import pyplot as plt
import sys

numberOfSamples = sys.argv[1]
sequenceLength = sys.argv[2]
inputFile = sys.argv[3]
plotFile = sys.argv[4]

with open(inputFile, 'r') as fin:
    tsv_reader = csv.reader(fin, delimiter="\t")
    next(tsv_reader)
    next(tsv_reader)
    iteration=[]
    prior=[]
    prior2=[]
    likelihood=[]
    joint=[]
    recombs=[]
    noncompats=[]
    arglen=[]
    
    for row in tsv_reader:
        iteration.append(int(row[1]))
        prior.append(float(row[2]))
        prior2.append(float(row[3]))
        likelihood.append(float(row[4]))
        joint.append(float(row[5]))
        recombs.append(int(row[6]))
        noncompats.append(int(row[7]))
        arglen.append(float(row[8]))

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
    plt.savefig(plotFile)
