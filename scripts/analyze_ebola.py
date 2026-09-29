import os
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def load_alignment(fasta_path):
    """Loads FASTA alignment file using Biopython."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Alignment file not found at: {fasta_path}")
    
    alignment = AlignIO.read(fasta_path, "fasta")
    print(f"[INFO] Successfully loaded alignment with {len(alignment)} sequences.")
    print(f"[INFO] Alignment length: {alignment.get_alignment_length()} bp\n")
    return alignment

def analyze_kj660348_mutations(alignment, results_path, target_id="KJ660348"):
    """
    Dynamically computes pairwise sequence mismatches and calculates 
    phylogenetic distance metrics using Biopython.
    """
    num_seqs = len(alignment)
    align_len = alignment.get_alignment_length()
    
    # Locate target sequence
    target_seq = None
    for record in alignment:
        if target_id in record.id:
            target_seq = record
            break
            
    if target_seq is None:
        target_seq = alignment[0]
        
    # 1. Dynamic Pairwise Mismatch Calculation
    total_mismatches_all_pairs = 0
    pair_count = 0
    
    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1 = str(alignment[i].seq)
            seq2 = str(alignment[j].seq)
            mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')
            total_mismatches_all_pairs += mismatches
            pair_count += 1

    # 2. Dynamic Phylogenetic Tree & Terminal Branch Metric Calculation
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    
    # Calculate sum of terminal branch lengths (proportion of substitution sites)
    terminal_branch_sum = sum(clade.branch_length for clade in tree.get_terminals() if clade.branch_length)
    estimated_terminal_substitutions = int(round(terminal_branch_sum * align_len))

    summary_content = (
        "============================================================\n"
        "  EBOLA VIRUS 2014 PHYLOGENOMIC MUTATION ANALYSIS SUMMARY   \n"
        "============================================================\n"
        f"Sequences Analyzed: {num_seqs}\n"
        f"Alignment Length: {align_len} bp\n"
        f"Target Accession: {target_seq.id}\n"
        f"Total Pairwise Mismatches (All Taxa Pairs): {total_mismatches_all_pairs}\n"
        f"Estimated Terminal Substitution Count (NJ Tree): {estimated_terminal_substitutions}\n"
        "============================================================\n"
    )

    print(summary_content)

    # Save output to results folder
    with open(results_path, "w") as f:
        f.write(summary_content)
    print(f"[INFO] Saved output results to: {results_path}")

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    data_path = os.path.join(base_dir, "data", "ebola_sequences.fasta")
    if not os.path.exists(data_path):
        data_path = os.path.join(base_dir, "data", "ebola_sequences.meg")

    results_path = os.path.join(base_dir, "results", "analysis_summary.txt")

    try:
        alignment = load_alignment(data_path)
        analyze_kj660348_mutations(alignment, results_path)
    except Exception as e:
        print(f"[ERROR] {e}")