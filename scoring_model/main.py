import argparse
from scoring_model_interface import get_seq, get_TSSs, get_all_sgRNA_sequences, conv_strand
from structuring_data import structure_dataframe
from scoring_model import score_sgRNAs
import pandas as pd

def main(targeted_gene, longest_dist = 1000, dir = "../genomes/c_elegans_n2"):    
    genome = get_seq(dir)
    # prim_TSS, sec_TSS, strand = get_TSSs(targeted_gene)
    prim_TSS = 44281043
    sec_TSS = 44281043
    strand = "-"
                                         
    potential_sgRNAs = get_all_sgRNA_sequences(prim_TSS, sec_TSS, longest_dist, genome, conv_strand(strand))
    
    structured_df = structure_dataframe(potential_sgRNAs.copy(), dir)
    
    score_df = score_sgRNAs(structured_df, strand)
    outcome_df = pd.concat([potential_sgRNAs.reset_index(drop = True), pd.DataFrame({"score": score_df})], axis = 1)
    outcome_df = outcome_df.sort_values(by = "score", ascending = False)
    
    outcome_df.to_csv("temp_dir/outcome.csv", index = False)
    
def parse_args():
    parser = argparse.ArgumentParser(description = "testing")
    
    # required arguments
    parser.add_argument("targeted_gene", type = str, help = "Gene to be targeted")
    
    # optional arguments
    parser.add_argument("--longest_dist", type = int, default = 1000, help = "Longest distance from TSS")
    parser.add_argument("--dir", type  = str, default = "../genomes/c_elegans_n2")
    
    return parser.parse_args()
    
if __name__ == "__main__":
    # main("Y53C10A.12.1")
    
    args = parse_args()
    main(args.targeted_gene, args.longest_dist, args.dir)