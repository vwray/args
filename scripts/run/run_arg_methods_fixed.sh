module load conda
conda init bash
source ~/.bashrc
conda activate arg
numberOfSamples=3
diploidSamples=6
sequenceLength="1_000_000"
sequenceLenInt=1000000
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
                "${outputDirectory}/msprimeTrees.trees" \
                ${recombRate}

#RENT+

#convert VCF to RENT+ input
python ${scriptRepo}/rentPlus/vcfToRentInput.py "${outputDirectory}/msprimeVCF.vcf" \
                 "${outputDirectory}/rentOutput.txt"

#run RENT+ with branch lengths
cd ${outputDirectory}
java -Xmx24g -jar /fs/cbcb-lab/ekmolloy/group/software/RentPlus/RentPlus.jar -t rentOutput.txt

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
      --seed 94577 \
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
python /fs/cbcb-lab/ekmolloy/group/software/ARGsims/scripts/argweaver/2_vcf2sites.py msprimeVCF.vcf ${diploidSamples} ${sequenceLenInt} ${outputDirectory}
#run ARGweaver  can add: --smc-prime \
conda activate py2
/fs/cbcb-lab/ekmolloy/group/software/argweaver-fork/argweaver/bin/arg-sample -s msprimeVCF.sites \
                 -N 10000 \
                 -r ${recombRate} \
                 -m ${mutationRate} \
                 --ntimes 20.0 \
                 --maxtime 200e4 \
                 -c 10 \
                 -n 2000 \
                 --randseed 632 \
                 -o argweaverOutput/arg-sample

#extract ARGweaver TMRCA
/fs/cbcb-lab/ekmolloy/group/software/argweaver-fork/argweaver/bin/arg-extract-tmrca argweaverOutput/arg-sample.%d.smc.gz \
    > argweaverTMRCA.txt

#determine ARGweaver local consensus trees
/fs/cbcb-lab/ekmolloy/group/software/argweaver-fork/argweaver/bin/arg-cons argweaverOutput/arg-sample.%d.smc.gz -s 1980 -d 1 > argweaverConsensusTrees.txt

#convert ARGweaver consensus trees to Newick and breakpoints
conda activate arg
python ${scriptRepo}/argweaver/convertArgweaverConsensusTrees.py ${numberOfSamples} \
                  ${sequenceLength} \
                  argweaverConsensusTrees.txt \
                  argweaverNewick.txt \
                  argweaverBreakpoints.txt

python ${scriptRepo}/argweaver/plotArgweaverConvergence.py ${numberOfSamples} \
                  ${sequenceLength} \
                  argweaverOutput/arg-sample.stats \
                  argweaverConvergencePlot.png

#unzip argweaverOutput/arg-sample.2000.smc.gz -d argweaverOutput/lastArgFolder
gzip -dk argweaverOutput/arg-sample.2000.smc.gz
conda activate py2
# get last ARGweaver ARG
/fs/cbcb-lab/ekmolloy/group/software/argweaver-fork/argweaver/bin/smc2arg argweaverOutput/arg-sample.2000.smc argweaverOutput/arg-sample.2000.arg

conda activate arg
python /fs/cbcb-lab/ekmolloy/vwray/args/gitRepoScripts/what-is-an-arg-paper/convertArgweaverOutput.py argweaverOutput/arg-sample.2000.arg argweaverLastARG.trees

# Compare trees from all 3 methods
python ${scriptRepo}/compareEstimatedToTrueTrees_argweaverLastARGandConsensusTrees.py ${numberOfSamples} \
                  ${sequenceLength} \
                  msprimeNewick.newick \
                  msprimeBreakpoints.txt \
                  relateNewick.txt \
                  relateBreakpoints.txt \
                  rentOutput.txt.trees \
                  rentNewick.txt \
                  argweaverNewick.txt \
                  argweaverBreakpoints.txt \
                  argweaverLastARG.trees \
                  argweaverLastARGNewick.txt \
                  treeComparison.png \
                  msprimeTrees.trees \
                  averageRFDistance.csv

# Compare TMRCA from all 3 methods
python ${scriptRepo}/compareTMRCA.py ${numberOfSamples} \
                  ${sequenceLength} \
                  msprimeNewick.newick \
                  msprimeBreakpoints.txt \
                  relateNewick.txt \
                  relateBreakpoints.txt \
                  rentOutput.txt.trees \
                  rentNewick.txt \
                  argweaverTMRCA.txt \
                  treeTMRCA.png \
                  10000


# Compare TMRCA error from all 3 methods
python ${scriptRepo}/compareTMRCAError.py ${numberOfSamples} \
                  ${sequenceLength} \
                  msprimeNewick.newick \
                  msprimeBreakpoints.txt \
                  relateNewick.txt \
                  relateBreakpoints.txt \
                  rentOutput.txt.trees \
                  rentNewick.txt \
                  argweaverTMRCA.txt \
                  treeTMRCAerror.png \
                  10000 \
                  averageTMRCAError.csv
