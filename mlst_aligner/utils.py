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


def save_results(results, file_path):
    # Function to save the computed results to a file.
    pass
