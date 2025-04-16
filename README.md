# Bioengineering S2025 Capstone Group J - Computational Model for Scoring sgRNAs for CRISPRi/a with Experimental Validation in *C. elegans*

## Overview

A computational model has been developed for scoring sgRNAs for CRISPRi/a with experimental validation in *C. elegans*. The model is first fine tuned on prior literature, where ~13 features are derived from the sgRNA sequence and target gene, relating to some performance score. Then, experimental data is to be generated to fine tune the model towards this data. Not enough data was generated within the scope of the experiment to fully fine tune the model for *C. elegans*; however, the pereliminary model framework is given.

## Scoring Model Workflow

Scoring model is split into the following files:

- Directories
  - data: containing all prior training and testing data
  - models: containing all fitted models
- Jupyter Notebook Files
  - structuring_data.ipynb: initial notebook for looking at ways to generate features from sequences
  - ML_models.ipynb: initial notebook for looking at different model structures
  - validiting_exp_data.ipynb: exploratory analysis on experimental data.
- Python files
  - scoring_model_interface.py: reads in a target gene and generates all possible sgRNAs based on the PAM (NGG) and 1000 bps in either direction of the TSS of the target gene - returns all of these to be used in later functions
  - structuring_data.py: functions used in structuring_data - takes in all sgRNAs and sequences and creates all features
  - scoring_model.py: contains the scoring model itself and applies scores based on dataframes structured using functions outlined in structuring_data.py
  - main.py: has the option of being run directly or called from commmand line. Takes in a target gene (as well as many other optional arguments specifying directory of genome and number of bases looked at in either direction) and generates a ranked list of sgRNAs corresponding to the fitted model's score for each sgRNA detected within a specified window of the specified gene.
- Misc.
  - cas-offinder.exe: executable to run cas-offinder for off-target scoring of sgRNAs

IMPORTANT NOTES:
- Genomes are required to be downloaded separately from this repository. *C. elegans* genome used for this model is N2. Model directory is specified in main.py using the "--dir" argument. Default directory is "../genomes/c_elegans_n2". For easiest installation of genome FASTA files, contact Jackson Rutecki at jackson.rutecki@gmail.com.

## Usage of the model

Usage of the model should be as simple as running main.py. As is given in the main method, there are two ways to run the file - in command line, or running directly in a compiler. Additionally, there are two optional arguments:
- --longest_dist: the longest distnance looked at (in bps) from the TSS of the targeted gene
- --dir: directory of genome

Name of gene is specified as a required argument of the script.