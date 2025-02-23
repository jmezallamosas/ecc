# interface for usage of the scoring model itself

import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np

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

def __get_seqs__(seq, strand, gene_strand, SHORTEST_SGRNA, LONGEST_SGRNA, prim_TSS, sec_TSS, max_distance):
    output_list = []
    
    # 1 for same strand, -1 for opposite
    same_strand = (strand == gene_strand) * 2 - 1
    
    parse_range = range(SHORTEST_SGRNA + 3, len(seq))
    if strand == 1:
        parse_range = reversed(parse_range)
    
    for i in parse_range:
        if (seq[i] == "G" and seq[i - strand] == "G"):
            current_seq = seq[np.max((i - LONGEST_SGRNA - 2, 0)):i - 1]
            
            for j in range(LONGEST_SGRNA - SHORTEST_SGRNA):
                if current_seq[j] == "G":
                    pam_loc = prim_TSS + ((max_distance - i) * -strand) - np.abs(strand - 1)
                    output_list.append({
                        "seq": "".join(current_seq[j:]),
                        "strand": strand,
                        "PAMloc": pam_loc,
                        "PAMdist_prim": same_strand * strand * (pam_loc - prim_TSS), # PLEASE double check this
                        "PAMdist_sec": same_strand * strand * (pam_loc - sec_TSS),
                    })
                    
    return output_list
                

# given the location of the primary and secondary TSS, determine all possible sgRNA locations within a certain distance
def get_all_sgRNA_sequences(prim_TSS, sec_TSS, max_distance, genome, gene_strand, SHORTEST_SGRNA = 16, LONGEST_SGRNA = 26):
    seq_front = genome[(prim_TSS - max_distance - 1):(prim_TSS + max_distance - 1)]
    seq_rev = seq_front.reverse_complement()
    
    forward_list = __get_seqs__(seq_front, 1, gene_strand, SHORTEST_SGRNA, LONGEST_SGRNA, prim_TSS, sec_TSS, max_distance)
    reverse_list = __get_seqs__(seq_rev, -1, gene_strand, SHORTEST_SGRNA, LONGEST_SGRNA, prim_TSS, sec_TSS, max_distance)
    
    return pd.concat([pd.DataFrame(forward_list), pd.DataFrame(reverse_list)])

CONV_STRAND = {"+": 1, "-": -1}

def main():
    genome_dir = "../genomes/c_elegans_n2"
    gene_name = "Y53C10A.12.1" # gene name for hsf-1    
    
    genome = get_seq(genome_dir)
    prim_TSS, sec_TSS, strand = get_TSSs(gene_name)
    
    potential_sgRNAs = get_all_sgRNA_sequences(prim_TSS, sec_TSS, 1000, genome, CONV_STRAND[strand])
    potential_sgRNAs.to_csv("testing.csv", index = False)

if __name__ == "__main__":
    main()