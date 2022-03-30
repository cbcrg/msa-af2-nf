# msa-af2-nf
Data, documentation, analysis and nextflow pipeline for the computation of structure-based MSAs with AlphaFold2 models

## Credits
This work has been carried out in [Notredame Lab](https://github.com/cbcrg) at the [Centre for Genomic Regulation - CRG](https://www.crg.eu/)

The authors who contributed to the analysis and manuscript are:

* Athanasios Baltzis
* Leila Mansouri
* Suzanne Jin
* Bjorn Langer
* Ionas Erb
* Cedric Notredame

## Notebooks
This repository contains a series of [Jupyter Notebooks](https://github.com/athbaltzis/msa-af2-nf/tree/main/notebook) that contain the steps for replicating the analysis, tables and figures in the manuscript using [R](https://www.r-project.org/).

## Pipeline and containers
The pipeline for predicting the AF2 models and producing the MSAs is built using [Nextflow](https://www.nextflow.io/). It comes with a singularity container (the recipe is available [here](https://github.com/athbaltzis/msa-af2-nf/blob/main/containers/AF2.def)) for running AF2 and a docker container (available on DockerHub [here](https://hub.docker.com/r/athbaltzis/pred)).

## Usage
- Download the genetic databases required for [AlphaFold2](https://github.com/deepmind/alphafold) using the provided [script](https://github.com/deepmind/alphafold/blob/main/scripts/download_all_data.sh).
- Download and format the database used for PSI-Coffee blast search (by default [Uniref50](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz)).
- Make sure you have singularity installed in your system.
- Install the Nextflow runtime by running the following command:

	$ curl -fsSL get.nextflow.io | bash

- You can launch the pipeline execution by entering the command shown below:
	
	$ nextflow run athbaltzis/msa-af2-nf

By default the pipeline is executed against the provided example [dataset](https://github.com/athbaltzis/msa-af2-nf/tree/main/data). You can modify the input data as well as the other available parameteres listed below:

#### `--input_fasta` 
Input sequences (FASTA)
#### `--list` 
Input lists of sequences
#### `--template` 
Input template lists
#### `--pdbs`
Input experimentally determined PDB structures
#### `--db`
Input path to Database for PSI-Coffee
#### `--predict`
Predict structures with AF2 [true or false(default)]
#### `--AF2`
Path to AF2 predicted models (if --predict false)
#### `--pdb_for_dssp`
Input PDB structures for secondary structure assignment
