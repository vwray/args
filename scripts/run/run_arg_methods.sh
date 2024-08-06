module load conda
conda activate arg
numberOfSamples=3
sequenceLength="100_000"
recombRate=1e-8
mutationRate=1.2e-8
baseDirectory="/fs/cbcb-lab/ekmolloy/vwray/args"
scriptRepo="${baseDirectory}/gitRepoScripts/args/scripts"
outputDirectory="${baseDirectory}/output"

cd ${outputDirectory}

newDirectory=$(date +"%Y-%m-%dT%H:%M:%S")_${sequenceLength}

mkdir ${newDirectory}

outputDirectory=${outputDirectory}/${newDirectory}

#call the Python script to generate msprime VCF, Newick, and breakpoint files
python ${scriptRepo}/generateMsprimeInput.py ${numberOfSamples} \
                ${sequenceLength} \
                "${baseDirectory}/simulationData/simple_ooa_neanderthal5r19_simp.yaml" \
                "${outputDirectory}/msprimeVCF.vcf" \
                "${outputDirectory}/msprimeNewick.newick" \
                "${outputDirectory}/msprimeBreakpoints.txt" \
                ${recombRate}

#RENT+

#convert VCF to RENT+ input
python ${scriptRepo}/rentPlus/vcfToRentInput.py "${outputDirectory}/msprimeVCF.vcf" \
                 "${outputDirectory}/rentOutput.txt"

#run RENT+ with branch lengths
cd ${outputDirectory}
java -jar /fs/cbcb-lab/ekmolloy/group/software/RentPlus/RentPlus.jar -t rentOutput.txt
 
#Relate

#convert VCF to haps and sample for Relate input
cd ${outputDirectory}
/fs/cbcb-lab/ekmolloy/group/software/relate/bin/RelateFileFormats \
                 --mode ConvertFromVcf \
                 --haps relateHaps.haps \
                 --sample relateSample.sample \
                 -i msprimeVCF        

#create genetic map for Relate input
python ${scriptRepo}/relate/createGeneticMap.py ${sequenceLength} \
                 ${recombRate} \
                 "${outputDirectory}/relateGeneticMap.txt"

#run Relate
cd ${outputDirectory}
/fs/cbcb-lab/ekmolloy/group/software/relate/bin/Relate \
      --mode All \
      -m $mutationRate \
      -N 20000 \
      --haps relateHaps.haps \
      --sample relateSample.sample \
      --map relateGeneticMap.txt \
      --seed 48 \
      -o relateOutput

#convert Relate output to tree sequence
/fs/cbcb-lab/ekmolloy/group/software/relate/bin/RelateFileFormats \
                 --mode ConvertToTreeSequence \
                 -i relateOutput \
                 -o relateOutput

#convert Relate tree sequence to Newick and breakpoints
python ${scriptRepo}/relate/convertRelateOutput.py ${numberOfSamples} \
                  ${sequenceLength} \
                  relateOutput.trees \
                  relateNewick.txt \
                  relateBreakpoints.txt

#ARGweaver

#convert VCF to ARGweaver sites input
cd ${outputDirectory}
python ${scriptRepo}/argweaver/convertVcfToArgweaverInput.py msprimeVCF.vcf argweaverSites.sites

#run ARGweaver
conda activate py2
/fs/cbcb-lab/ekmolloy/group/software/argweaver/argweaver/bin/arg-sample -s argweaverSites.sites \
                 -N 10000 \
                 -r ${recombRate} \
                 -m ${mutationRate} \
                 --ntimes 20.0 \
                 --maxtime 200e4 \
                 -c 10 \
                 -n 2000 \
                 --randseed 632 \
                 -o argweaverOutput

#extract ARGweaver TMRCA
/fs/cbcb-lab/ekmolloy/group/software/argweaver/argweaver/bin/arg-extract-tmrca argweaverOutput.%d.smc.gz \
    > argweaverTMRCA.txt

#determine ARGweaver local consensus trees
/fs/cbcb-lab/ekmolloy/group/software/argweaver/argweaver/bin/arg-cons argweaverOutput.%d.smc.gz -s 1000 -d 20 > argweaverConsensusTrees.txt

#convert ARGweaver consensus trees to Newick and breakpoints
conda activate arg
python ${scriptRepo}/argweaver/convertArgweaverConsensusTrees.py ${numberOfSamples} \
                  ${sequenceLength} \
                  argweaverConsensusTrees.txt \
                  argweaverNewick.txt \
                  argweaverBreakpoints.txt

python ${scriptRepo}/argweaver/plotArgweaverConvergence.py ${numberOfSamples} \
                  ${sequenceLength} \
                  argweaverOutput.stats \
                  argweaverConvergencePlot.png 


# Compare trees from all 3 methods
python ${scriptRepo}/compareEstimatedToTrueTrees.py ${numberOfSamples} \
                  ${sequenceLength} \
                  msprimeNewick.newick \
                  msprimeBreakpoints.txt \
                  relateNewick.txt \
                  relateBreakpoints.txt \
                  rentOutput.txt.trees \
                  rentNewick.txt \
                  argweaverNewick.txt \
                  argweaverBreakpoints.txt \
                  treeComparison.png