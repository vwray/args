# -*- coding: utf-8 -*-
"""
Creates a genetic map file for Relate based on a sequence length and recombination rate.

Created on Thu Dec 28 13:58:21 2023

@author: Valerie Wray
"""

import sys

sequenceLength = int(sys.argv[1])
recombRate = float(sys.argv[2])
outputFile = sys.argv[3]

recombRateConvertedUnits = recombRate * 1e8
print(recombRateConvertedUnits)
geneticPosition = recombRateConvertedUnits*sequenceLength/1e6

with open(outputFile, 'w') as file:
    file.write("position COMBINED_rate.cM.Mb. Genetic_Map.cM.\n")
    file.write(f"0 {recombRateConvertedUnits} 0\n")
    file.write(f"{sequenceLength} {recombRateConvertedUnits} {geneticPosition}")