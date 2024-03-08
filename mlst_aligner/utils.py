"""utils.py"""
import os
from typing import List, Tuple, Union
from pysam import FastaFile


def read_fasta(file_path: str) -> Union[FastaFile, None]:
    """Read sequences from a FASTA file.

    Validates the provided file path and ensures it has a .fasta extension before reading. 
    If the file is valid, a FastaFile object is created and returned.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        FastaFile: An object representing the FASTA file and its contents.

    Raises:
        FileNotFoundError: If the file does not exist at the provided path.
        ValueError: If the file path does not point to a file with a .fasta or .fa extension.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    if not (file_path.lower().endswith('.fasta') or file_path.lower().endswith('.fa')):
        raise ValueError("File extension must be .fasta or .fa")

    sequences_object = FastaFile(file_path)
    return sequences_object

def fetch_references(file_path: str) -> List[Tuple[str, str]]:
    """
    Utilizes the read_fasta function to parse a FASTA file, extracting gene names and sequences.
    
    Parameters:
    - file_path: str, path to the FASTA file.
    
    Returns:
    List[Tuple[str, str]]: A list of tuples, where each tuple contains a gene name and its sequence.
    """
    try:
        fasta_obj = read_fasta(file_path)
        if fasta_obj is None:
            return []

        gene_sequences = []
        for reference in fasta_obj.references:
            sequence = fasta_obj.fetch(reference)
            gene_sequences.append((reference, sequence))
        
        return gene_sequences

    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading FASTA file: {e}")
        return []

def weighted_average(scores: List[Tuple[int, int]]) -> float:
    """
    Computes the weighted average for a list of values and their weights.

    Args:
        scores (List[Tuple[int, int]]): A list of tuples, where each tuple contains a value and its corresponding weight.

    Returns:
        float: The weighted average of the values.
    """
    if not scores:
        return 0

    total_product_sum = sum(value * weight for value, weight in scores)
    total_weight_sum = sum(weight for _, weight in scores)

    return total_product_sum / total_weight_sum if total_weight_sum else 0


def subset_fasta(original_fasta: str, subset_count: 1000, output_fasta: str):
    """
    Subsets a FASTA file and saves the subset to a new file.

    Parameters:
    - original_fasta: str, path to the original FASTA file.
    - subset_count: int, number of reads to include in the subset.
    - output_fasta: str, path where the subset FASTA file will be saved.

    The function reads the original FASTA file, extracts the specified number of reads,
    and writes them to the output FASTA file.
    """
    with FastaFile(original_fasta) as fasta:
        with open(output_fasta, 'w') as outfile:
            for i in range(min(subset_count, len(fasta.references))):
                sequence = fasta.fetch(fasta.references[i])
                outfile.write(f">{fasta.references[i]}\n{sequence}\n")

    print(f"Subset FASTA file saved to {output_fasta} with {subset_count} reads.")
