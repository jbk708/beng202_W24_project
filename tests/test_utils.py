"""utils.py"""

import pytest
import os
from pysam import FastaFile
from mlst_aligner.utils import read_fasta, weighted_average, subset_fasta
from unittest.mock import patch, mock_open


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


@pytest.mark.parametrize(
    "scores, expected_average",
    [
        # Typical case
        ([(10, 2), (20, 3), (30, 5)], (10 * 2 + 20 * 3 + 30 * 5) / (2 + 3 + 5)),
        # Equal weights (simple average)
        ([(1, 1), (2, 1), (3, 1)], (1 + 2 + 3) / 3),
        # Varying weights
        ([(100, 1), (200, 2), (300, 3)], (100 * 1 + 200 * 2 + 300 * 3) / (1 + 2 + 3)),
        # Empty list
        ([], 0),
        # This test case is commented out because normally weights should be positive.
        # ([(1, -1), (2, 1)], 0),
    ])
def test_weighted_average(scores, expected_average):
    assert weighted_average(scores) == pytest.approx(
        expected_average), "The weighted average does not match the expected value."


def test_subset_fasta_creates_output_file(tmp_path):
    """
    Test that the subset_fasta function creates an output file at the specified location.
    """
    original_fasta = "path/to/original/test_data.fasta"
    output_fasta = tmp_path / "output.fasta"
    subset_count = 5

    # Mock open to simulate reading and writing files
    with patch("builtins.open", mock_open(read_data=">seq1\nATGC\n>seq2\nATGG\n")) as mocked_file:
        # Assuming FastaFile is a custom class you're using to read FASTA files
        with patch("mlst_aligner.utils.FastaFile") as MockFastaFile:
            mock_fasta_file = MockFastaFile.return_value
            mock_fasta_file.references = ['seq1', 'seq2']
            mock_fasta_file.fetch.side_effect = ['ATGC', 'ATGG']

            subset_fasta(original_fasta, subset_count, str(output_fasta))
            mocked_file.assert_called_with(str(output_fasta),
                                           'w')  # Check if the output file was attempted to be opened for writing
