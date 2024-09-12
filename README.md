A project to benchmark methods for inferring Ancestral Recombination Graphs (ARGs).

## Prerequisites

It is recommended to use a compute cluster for running the scripts in this project. ARGweaver, Relate, and RENT+ must be installed. Conda is used in the scripts to manage dependencies and environment. Conda expects an environment `arg` with Python 3 installed and other dependencies required by `tskit` and `msprime`. Conda also expects an environment `py2` with Python 2 installed, which is required by ARGweaver. The scripts in this repo switch between the two Conda environments.

## Running the scripts

Use the script [run_arg_methods_fixed.sh](scripts/run/run_arg_methods_fixed.sh) to generate msprime simulated data, run ARGweaver, Relate, and RENT+, and plot accuracy results with Python scripts. Modify the parameters at the top to change, sequence length, recombination rate, etc. In addition, a base directory containing an `output` directory should be specified. These scripts have been run using the base directory `/fs/cbcb-lab/ekmolloy/vwray/args` on the CBCB cluster. These scripts create a new timestamped directory inside the output directory during each run.

Use the script [run_arg_methods_deboraData_ARGweaverModified](scripts/run/run_arg_methods_deboraData_ARGweaverModified) for recreating the benchmarking done in Debora Brandt's evaluation of ARG methods paper.

Both of the above scripts call the other Python scripts in this repo for generating msprime data, converting file types, computing accuracy, and plotting results.
