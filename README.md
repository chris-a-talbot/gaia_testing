## Pipeline for a single replicate

### Overview

`run_slim.slim` -> `coalescent_check.py` -> `generate_samples.py` -> `subset_tree.py`

-> `generate_samples.py --subsets` -> `validate_pre_gaia.py` -> `run_gaia.R` -> `analyze_gaia.py`

### Steps

1. Run `run_slim.slim` to generate the tree sequence.
2. Run `coalescent_check.py` to ensure SLiM simulation ran to coalescence. Returns `0` if successful and `1` otherwise.
3. Run `generate_samples.py` to generate sample location data from tree sequence.
4. Run `subset_tree.py` to subset the tree sequence to various sampling schemes.
5. Run `generate_samples.py` to generate sample location data from subset tree sequences.
6. Run `validate_pre_gaia.py` to validate the sample location data. Returns `0` if there are no issues and `1` otherwise.
7. Run `run_gaia.R` to run gaia on the subset tree sequences and sample location data.
8. Run `analyze_gaia.py` to analyze the accuracy of geographic inference from gaia.

### Example
To run the pipeline for a single replicate named 'tree-S0.3-R1', run the following commands.

Replacing `<replicate ID>` with `1`, `<dispersal>` with `0.3`, and `<CWD>` with the current working directory:

`slim -d "PWD=<CWD>" -d "S=<dispersal>" -d "REP=<replicate ID>" -d "RUNTIME=40000" ./scripts/run_slim.slim`

`python ./scripts/coalescent_check.py "tree-S0.3-R1"`

`python ./scripts/generate_samples.py "tree-S0.3-R1"`

`python ./scripts/subset_tree.py "tree-S0.3-R1"`

`python ./scripts/generate_samples.py "tree-S0.3-R1" --subsets`

`python ./scripts/validate_pre_gaia.py "tree-S0.3-R1"`

`Rscript ./scripts/run_gaia.R "tree-S0.3-R1"`

`python ./scripts/analyze_gaia.py "tree-S0.3-R1"`