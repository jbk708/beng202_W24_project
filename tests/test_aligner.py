"""test_aligner.py"""
import os
import pytest
from mlst_aligner.aligner import positional_alignment

@pytest.mark.parametrize("match_reward, mismatch_penalty, indel_penalty, s, t, expected_score, expected_aligned_s, expected_aligned_t", [
    (3, 3, 1, "AGC", "ATC", 4, "A-GC", "AT-C"),
    (1, 1, 1, "TAACG", "ACGTG", 3, "ACG", "ACG"),
    (3, 2, 1, "CAGAGATGGCCG", "ACG", 6, "CG", "CG"),
    (2, 3, 1, "CTT", "AGCATAAAGCATT", 5, "C-TT", "CATT"),
])
def test_positional_alignment(match_reward, mismatch_penalty, indel_penalty, s, t, expected_score, expected_aligned_s, expected_aligned_t):
    # Run the alignment
    score, aligned_s, aligned_t, scores_at_positions = positional_alignment(match_reward, mismatch_penalty, indel_penalty, s, t)

    # Assert that the returned score and alignments are as expected
    assert score == expected_score
    assert aligned_s == expected_aligned_s
    assert aligned_t == expected_aligned_t
    print(scores_at_positions)