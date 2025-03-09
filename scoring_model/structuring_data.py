import pandas as pd
import subprocess
import numpy as np
import ViennaRNA as vr

# each base pair
BPS = ["A", "C", "T", "G"]

# combinations of each base pair
COMBO_BPS = []
for i, b1 in enumerate(BPS):
    for j, b2 in enumerate(BPS):
        COMBO_BPS.append(b1 + b2)

# each part of the folding structure
PARTS = [".", "(", ")"]
    
# add GC content to the data frame
def adding_GC_content(df):
    for base in ("A", "C", "T", "G"):
        df[("frac_" + base)] = df.apply(lambda row: row["sequence"].count(base) / row["length"], axis = 1)
    
    df["GC_content"] = df["frac_C"] + df["frac_G"]
    df["GA_content"] = df["frac_G"] + df["frac_A"]
    df["CA_content"] = df["frac_C"] + df["frac_A"]
    
    return df

# encode each monomer with one hot encoding
def __encode_monomer__(seq, max_len):
    one_hot_vector = np.zeros((max_len, len(BPS)))
    
    for i, bp in enumerate(seq):
        if bp in BPS:
            one_hot_vector[i, BPS.index(bp)] = 1
    return one_hot_vector

# add one hot encoding of each monomer to the dataframe
def one_hot_encode_monomers(df, max_len):
    encoded_monomers = df["sequence"].apply(lambda x: __encode_monomer__(x, max_len))
    temp_cols = {}

    for i in range(max_len):
        for j, base in enumerate(BPS):
            temp_cols[f"pos_{i+1}_{base}"] = encoded_monomers.apply(lambda x: x[i, j])
            
    df = pd.concat([df, pd.DataFrame(temp_cols)], axis = 1)
    
    return df

# encode each dimer with one hot encoding
def __encode_dimer__(seq, max_len):
    one_hot_vector = np.zeros((max_len, len(COMBO_BPS)))
    
    # iterating over two bases at a time, in order
    for i in range(len(seq) - 1):
        dimer = seq[i] + seq[i + 1]
        if dimer in COMBO_BPS:
            one_hot_vector[i, COMBO_BPS.index(dimer)] = 1
    return one_hot_vector

# add one hot encoding for each dimer to the dataframe
def one_hot_encode_dimers(df, max_len):
    encoded = df["sequence"].apply(lambda x: __encode_dimer__(x, max_len))
    temp_cols = {}

    for i in range(max_len - 1): # taking 1 off as we add 2 later on, so adding 2 would get to 26
        for j, bases in enumerate(COMBO_BPS):
            temp_cols[f"pos_{i+1}{i+2}_{bases}"] = encoded.apply(lambda x: x[i, j])
            
    df = pd.concat([df, pd.DataFrame(temp_cols)], axis = 1)
    
    return df

# normalizes sequence strings to be the same length as the max length sgRNA, as needed for COF
def __normalize_str_len__(str, max_len):
    while len(str) < max_len:
        str = "N" + str
    return(str)

def __construct_cof_input__(df, max_len, dir):
    # normalize all sgRNAs to have the same length as the max length sgRNA (by adding Ns in front)
    df["sequence"] = df["sequence"].apply(__normalize_str_len__, args = (max_len,))

    # preparing input strings for cas-offinder
    pam = "NGG"
    full_pam = ("N" * max_len) + pam

    co_input = f"{dir}\n{full_pam} 2 1"
    for i, seq in enumerate(df["sequence"], start = 1):
        co_input = co_input + f"\n{seq}NNN 3 Seq{i}"

    # write input string to input.txt for cas-offinder
    with open("temp_dir/input.txt", "w") as file:
        file.write(co_input)
    
# function to determine what indices of a given string are lowercase
def __lowercase_indices__(string):
    indices = []
    
    for i, c in enumerate(string):
        if c.islower():
            # adding 1 to the index so it is the position in sgRNA (based at 1, not 0)
            indices.append(i + 1)

    return indices      

# function to update a temp df with the running count of the number of MMs at a given position for each sequence
# - row: the row from the OT output
# - row_index: the index of the row from the OT output
# - df: the temp dataframe storing all info about positions and MMs at each position
# - max_len: constant, maximum length sgRNA found in the sequence
def __populate_temp_df__(row, row_index, df, max_len):
    num_Ns = row["WT"].count("N") - 3 # subtract 3 because of PAM
    lowercases = pd.Series(__lowercase_indices__(row["OT"])) - num_Ns # calculate indices of lowercases, subtract number of Ns
    
    for j in lowercases:
        # some MMs can be in the PAM - we don't care about these
        if j <= max_len:
            # add 1 to the running count if a mismatch in a position is found
            df.at[row_index, f"num_MMs_pos{j}"] += 1

# using CasOFFinder for off-target analysis
def run_cof(df, max_len, dir):
    # construct input for and run casoffinder
    __construct_cof_input__(df.copy(), max_len, dir)
    subprocess.run(f"scoring_model/cas-offinder temp_dir/input.txt G temp_dir/output.txt")
    
    # read the output from casoffinder and map them to the dataframe
    co_output = pd.read_csv("temp_dir/output.txt", sep = "\t", comment = "#", header = None)
    co_output.columns = ["WT", "chromosome", "posn", "OT", "dir", "mismatches"]

    temp_df = {}

    # loop through each all numbers of mismatches and count how many each sequence has for each number
    for MM in (0, 1, 2, 3):
        temp_df[f"{MM}_MMs"] = df["sequence"].apply(lambda seq: len(co_output.loc[(co_output["WT"].str.contains(seq)) &
                                                                             (co_output["mismatches"] == MM)]))
    df = pd.concat([df, pd.DataFrame(temp_df)], axis = 1)

    # initializing a temp dataframe with 0s for mismatches between position 1 and position max_len for all sequences
    temp_df = pd.DataFrame()
    for i in range(max_len):
        temp_df[f"num_MMs_pos{(i + 1)}"] = np.zeros(len(df))
        
    # iterate through the actual dataframe, looking at each sequence and counting the number of MMs at each position using the above defined functions
    for i, seq in enumerate(df["sequence"]):
        co_output[co_output["WT"].str.contains(seq)].apply(lambda row: __populate_temp_df__(row, i, temp_df, max_len), axis = 1)
        
    df = pd.concat([df.reset_index(drop = True), pd.DataFrame(temp_df)], axis = 1)
    
    return df

# determine the MFE of the folding structure of a sequence
def __structure_mfe__(seq):
    fc = vr.fold_compound(seq)
    return fc.mfe()

def __encode_struct__(struct, max_len):
    one_hot_vector = np.zeros((max_len, len(PARTS)))
    
    for i, part in enumerate(struct):
        if part in PARTS:
            one_hot_vector[i, PARTS.index(part)] = 1
    return one_hot_vector

# add structure and MFE for each sequence
def sgRNA_struct(df, max_len):
    # determine MFE and folding structure
    structs_mfes = df["sequence"].apply(lambda x: pd.Series(__structure_mfe__(x)))
    df = pd.concat([df, pd.DataFrame({"struct": structs_mfes[0], "mfe": structs_mfes[1]})], axis = 1)
    
    # one hot encodde the secondary structure for each position
    encoded_structs = df["struct"].apply(lambda x: __encode_struct__(x, max_len))
    temp_cols = {}

    for i in range(max_len):
        for j, part in enumerate(PARTS):
            temp_cols[f"pos_{i+1}_{part}"] = encoded_structs.apply(lambda x: x[i, j])
            
    df = pd.concat([df, pd.DataFrame(temp_cols)], axis = 1)
    
    return df

# determine the longest homopolymer in a given sequence
def __longest_homopolymer_in_seq__(seq):
    max_count = {"longest_" + bp: 0 for bp in BPS}
    
    for bp in BPS:
        current_count = 0
        for base in seq:
            if bp == base:
                current_count += 1
                if current_count > max_count["longest_" + bp]:
                    max_count["longest_" + bp] = current_count
            else:
                current_count = 0
    return max_count

# for the dataframe, determine the longest homopolymer in each sequence
def longest_homopolymer(df): 
    temp_df = df["sequence"].apply(__longest_homopolymer_in_seq__)
    df = pd.concat([df, temp_df.apply(pd.Series)], axis = 1)
    
    return df

def structure_dataframe(df, dir):
    df["length"] = df["sequence"].apply(len)
    max_len = df["length"].agg(max)
    
    df = adding_GC_content(df)
    df = one_hot_encode_monomers(df, max_len)
    df = one_hot_encode_dimers(df, max_len)
    df = run_cof(df, max_len, dir)
    df = sgRNA_struct(df, max_len)
    df = longest_homopolymer(df)
    
    return df