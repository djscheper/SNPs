#!/usr/bin/env python3

"""
The main goal of this script is to introduce a single nucleotide morphism (SNP) on a given position in a gene sequence
specified by the user in the form of a FASTA file. The effect of the introduced SNP will be measured in
conservation against a multiple sequence alignment (MSA), which includes protein sequences closely
aligned with the given gene sequence. The system measures conservation in percentage and gives a final verdict
on whether an SNP is deleterious (<90% conservation) or has no effect (≥90% conservation).


Usage:
    python3 snp.py -n {'A', 'C', 'G', 'T'} -p POSITION -m MSA -s GENE_SEQUENCE
"""

__author__ = "Dennis Scheper (373689)"
__date__ = "05-11-2021"
__status__ = "Production"
__version__ = "v1.0"
__contact__ = "d.j.scheper@st.hanze.nl"


import argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import MutableSeq


codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}


def parse_arguments():
    """
    Defines command-line arguments.
    """
    parser = argparse.ArgumentParser(description="""
    """)
    parser.add_argument("-n", "--nucleotide",
                        help="""The Single Nucleotide Polymorphism (SNP) to take place.
                        Choices are: A, C, T or G.""",
                        required=True, choices=['A', 'C', 'G', 'T'])
    parser.add_argument("-p", "--position",
                        help="Position of SNP in the nucleotide sequence", required=True, type=int)
    parser.add_argument("-m", "--msa",
                        help="The location and name of the Multiple Sequence Alignment (MSA) FASTA file",
                        required=True)
    parser.add_argument("-s", "--dna_sequence",
                        help="""The location and name of the FASTA file
                        containing a dna sequence to be compared to the MSA""", required=True)
    return parser.parse_args()


def insert_snp(snp, position, fasta_file):
    """
    Reads in the gene sequence from a FASTA file and inserts the SNP.
    :param snp
    :param position
    :param fasta_file
    :returns gene_alignment
    """
    gene_alignment = SeqIO.read(fasta_file, "fasta")
    if position <= len(gene_alignment.seq):
        # Reassign to MutableSeq so we can replace the nucleotide in question
        gene_alignment = MutableSeq(gene_alignment.seq)
        print("\tEVENT:")
        print(f"Nucleotide {gene_alignment[position-1]} has been replaced by SNP {snp} on position {position}.")
        gene_alignment[position-1] = snp
    else:
        print(f"""Position {position} is out of range with the given DNA sequence, 
        which has a length of {len(gene_alignment)}.""")
    return gene_alignment


def translate_sequence(dna_sequence):
    """
    Translates the given DNA sequence to amino acids.
    :param dna_sequence
    :returns amino_sequence
    """
    # Split the sequence into codons
    codon_list = [str(dna_sequence[i:i+3]) for i in range(0, len(dna_sequence), 3)]
    # Translate to amino acids
    amino_sequence = [codon_table[codon] if codon in codon_table else "-" for codon in codon_list]
    return amino_sequence


def determine_outcome(amino_sequence, msa, position):
    """
    Determines the outcome of the insertion of the SNP by comparing it
    to a multiple alignment sequence (MSA). Check every alignment in
    the MSA to see if the SNP resulted in an entirely different amino acid
    composition.
    :param amino_sequence
    :param msa
    :param position
    :return outcome
    """
    aa_seq = AlignIO.read(msa, format="clustal")
    counter = 0
    for alignment in aa_seq:
        if alignment.seq[int((position-1)/3)] == amino_sequence[int((position-1)/3)]:
            counter += 1
    outcome = round(counter/len(aa_seq)*100, 2)
    return outcome

def main():
    """
    The main function of this script.
    """
    args = parse_arguments()
    # Reassign to more suitable variable names
    snp, pos, msa, seq = args.nucleotide, args.position, args.msa, args.dna_sequence
    if pos > 0:
        add_snp = insert_snp(snp, pos, seq)
        translated = translate_sequence(add_snp)
        outcome = determine_outcome(translated, msa, pos)
        print("\tOUTCOME:")
        print(f"The inserted SNP {snp} on position {pos} has {outcome}% similarity in comparison with the given MSA.")
        print(f"Therefore, the final verdict is that the SNP is {'deleterious' if outcome < 90 else 'neutral'}.")
    else:
        print(f"Execution halted due to position being set at {pos}. Only use numbers bigger than 0.")
    

if __name__ == "__main__":
    main()
