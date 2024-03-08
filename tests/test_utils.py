"""utils.py"""

import pytest
import os
from pysam import FastaFile
from mlst_aligner.utils import read_fasta  # Adjust import path as necessary


def create_temp_fasta_file(tmp_path, content=">seq1\nATCG"):
    """
    Helper function to create a temporary FASTA file.
    """
    file_path = tmp_path / "temp.fasta"
    with open(file_path, "w") as f:
        f.write(content)
    return file_path

def test_read_fasta_valid_file(tmp_path):
    """
    Test that read_fasta returns a FastaFile object for a valid FASTA file.
    """
    file_path = create_temp_fasta_file(tmp_path)
    result = read_fasta(str(file_path))
    assert isinstance(result, FastaFile), "The result should be a FastaFile object"

def test_read_fasta_file_not_found():
    """
    Test that read_fasta raises FileNotFoundError for a non-existing file.
    """
    with pytest.raises(FileNotFoundError):
        read_fasta("non_existing_file.fasta")

def test_read_fasta_invalid_extension(tmp_path):
    """
    Test that read_fasta raises ValueError for a file with invalid extension.
    """
    file_path = tmp_path / "invalid.txt"
    file_path.touch()  # Create an empty file with invalid extension
    with pytest.raises(ValueError):
        read_fasta(str(file_path))