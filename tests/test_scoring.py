"""test_scoring.py"""
import pytest
from unittest.mock import patch, MagicMock
from mlst_aligner.scoring import merge_scores, GeneScore

@pytest.mark.parametrize("input_dicts, expected_output", [
    # Case 1: Merging dictionaries with unique positions
    ([{1: 10, 2: 20}, {3: 30, 4: 40}], {1: [10], 2: [20], 3: [30], 4: [40]}),
    # Case 2: Merging dictionaries with overlapping positions
    ([{1: 10, 2: 20}, {1: 15, 2: 25, 3: 35}], {1: [10, 15], 2: [20, 25], 3: [35]}),
    # Case 3: Handling an empty list of dictionaries
    ([], {}),
    # Case 4: Single dictionary input
    ([{1: 10, 2: 20}], {1: [10], 2: [20]}),
])

def test_merge_scores(input_dicts, expected_output):
    assert merge_scores(input_dicts) == expected_output
