# script to verify that all of the sgRNAs designed are indeed in the c. elegans genome
# without any mismatches

import pandas as pd
import numpy as np
import subprocess

df = pd.read_csv("scoring_model/data/test_sgRNA_sequences.csv")
sequences = df.loc[:, "sequence"]

max_len = sequences.apply(len).max()

# pad strings with Ns until they are all the same length
def normalize_str_len(str, max_len):
    while len(str) < max_len:
        str = "N" + str
    return(str)

sequences = sequences.apply(normalize_str_len, args = (max_len,))

# preparing input strings for cas-offinder
pam = "NGG"
full_pam = ("N" * max_len) + pam

# construct input file string
co_input = f"../genomes/c_elegans_n2\n{full_pam} 0 0"
for i, seq in enumerate(sequences, start = 1):
    co_input = co_input + f"\n{seq}NNN 0 Seq{i}"
    
# write input string to input.txt for cas-offinder
with open("scoring_model/data/input_elegans.txt", "w") as file:
    file.write(co_input)

# running subprocess for cas-offinder
subprocess.run(f".\scoring_model/cas-offinder scoring_model/data/input_elegans.txt G scoring_model/data/output_elegans.txt")

# read the input and output .csv files again, and add sequence IDs onto the output
input_csv = pd.read_csv("scoring_model/data/input_elegans.txt", sep = " ", header = None, skiprows = 2)
input_csv.columns = ["sequence", "MMs", "seq_ID"]

output_csv = pd.read_csv("scoring_model/data/output_elegans.txt", sep = "\t", header = None)
output_csv.columns = ["sequence", "chromID", "loc", "MM_sequence", "strand", "MMs"]

# output the merged dataframe with sequence IDs on output as .csv
merged_csv = pd.merge(input_csv[["sequence", "seq_ID"]], output_csv, on = "sequence", how = "right")
merged_csv.to_csv("scoring_model/data/merged_csv.csv", index = False)