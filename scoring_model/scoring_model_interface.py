# interface for usage of the scoring model itself

import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np

SHORTEST_SGRNA = 18
LONGEST_SGRNA = 25

# import the c. elegans genome from multiple FASTA files into one single sequence
def get_seq(genome_dir):
    combined_seq = Seq("")

    # combine all fasta files into one singular sequence
    for filename in os.listdir(genome_dir):
        if filename.endswith("fa"):
            file_path = os.path.join(genome_dir, filename)
            
            for record in SeqIO.parse(file_path, "fasta"):
                combined_seq += record.seq
                
    return combined_seq

# get primary, secondary TSS and strand of the targeted gene
def get_TSSs(gene_name):
    chen_df = pd.read_excel("scoring_model/data/chen_c.elegans_TSS.xlsx", sheet_name = "all TSSs")

    TSSs = chen_df[chen_df["gene name"] == gene_name].reset_index(drop = True)

    # get the locations of primary and secondary TSS for the gene
    prim_TSS = TSSs.loc[TSSs["number of aligned reads in the best-fit Gaussian distribution"].idxmax(), :]["start position"]
    sec_TSS = TSSs.loc[TSSs["number of aligned reads in the best-fit Gaussian distribution"].nlargest(2).index[1], :]["start position"]

    # get the strand of the gene
    strand = np.unique(TSSs["strand"])
    
    if len(strand) > 1:
        raise ValueError("strands are non-unique, go back and check .csv input")
    
    return(prim_TSS, sec_TSS, strand[0])

def __get_seqs__(seq, strand, gene_strand, prim_TSS, sec_TSS, max_distance):
    output_list = []
    
    # 1 for same strand, -1 for opposite
    same_strand = (strand == gene_strand) * 2 - 1
    
    parse_range = range(SHORTEST_SGRNA + 3, len(seq))
    
    for i in parse_range:
        if (seq[i] == "G" and seq[i - strand] == "G"):
            # i have absolutely 0 idea why this is needed, but its the only way i could get it to work
            if strand == 1:
                shifter = -2
            else:
                shifter = -1
            current_seq = seq[np.max((i - LONGEST_SGRNA + shifter, 0)):i + shifter]
            
            for j in range(LONGEST_SGRNA - SHORTEST_SGRNA + 1):
                if current_seq[j] == "G":
                    pam_loc = prim_TSS + ((max_distance - i) * -strand) - np.abs(strand - 1)
                    obs_seq = "".join(current_seq[j:])
                    
                    # this if gate is to get rid of some edge cases where sequences were far too short at the ends of the observed region
                    if (len(obs_seq) >= SHORTEST_SGRNA) and (len(obs_seq) <= LONGEST_SGRNA):
                        output_list.append({
                            "sequence": obs_seq,
                            "strand": strand * -1, # inverse, because strand targeted is opposite of strand sequence is designed for
                            "PAMloc": pam_loc,
                            # including 5p and 3p distances as same for now, may change this in the future
                            "prim_TSS_dist5p": same_strand * strand * (pam_loc - prim_TSS), # PLEASE double check this
                            "prim_TSS_dist3p": same_strand * strand * (pam_loc - prim_TSS),
                            "sec_TSS_dist5p": same_strand * strand * (pam_loc - sec_TSS),
                            "sec_TSS_dist3p": same_strand * strand * (pam_loc - sec_TSS),
                        })
                    
    return output_list
                

# given the location of the primary and secondary TSS, determine all possible sgRNA locations within a certain distance
def get_all_sgRNA_sequences(prim_TSS, sec_TSS, max_distance, genome, gene_strand, ):
    seq_front = genome[(prim_TSS - max_distance - 1 - LONGEST_SGRNA):(prim_TSS + max_distance - 1 + LONGEST_SGRNA)]
    seq_rev = seq_front.reverse_complement()
    
    forward_list = __get_seqs__(seq_front, 1, gene_strand, prim_TSS, sec_TSS, max_distance + LONGEST_SGRNA)
    reverse_list = __get_seqs__(seq_rev, -1, gene_strand, prim_TSS, sec_TSS, max_distance + LONGEST_SGRNA)
    
    # return pd.DataFrame(forward_list)
    return pd.concat([pd.DataFrame(forward_list), pd.DataFrame(reverse_list)])

# better way to write this?
def conv_strand(str) :
    CONV_STRAND = {"+": 1, "-": -1}
    return CONV_STRAND[str]

def main():
    genome_dir = "../genomes/c_elegans_n2" # directory for c. elegans genome
    gene_name = "Y53C10A.12.1" # gene name for hsf-1
    
    genome = get_seq(genome_dir)
    prim_TSS, sec_TSS, strand = get_TSSs(gene_name)
    
    potential_sgRNAs = get_all_sgRNA_sequences(prim_TSS, sec_TSS, 1000, genome, conv_strand(strand))
    potential_sgRNAs.to_csv("scoring_model/data/all_possible_sgRNAs.csv", index = False) 
    # strand targeted != strand that the sgRNA is designed for, by the naming convention used here

if __name__ == "__main__":
    main()