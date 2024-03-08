"""test_aligner.py"""
import os
import pytest
from mlst_aligner.aligner import positional_alignment


@pytest.mark.parametrize(
    "match_reward, mismatch_penalty, indel_penalty, s, t, expected_score, expected_aligned_s, expected_aligned_t, expected_position_score", [
        (3, 3, 1, "AGC", "ATC", 4, "A-GC", "AT-C", {1: (1, 3), 3: (4, 3)}),
        (1, 1, 1, "TAACG", "ACGTG", 3, "ACG", "ACG", {2: (1, 3), 3: (3, 3), 4: (2, 3), 5: (1, 3)}),
        (3, 2, 1, "CAGAGATGGCCG", "ACG", 6, "CG", "CG", {2: (2, 2), 3: (6, 2)}),
        (2, 3, 1, "CTT", "AGCATAAAGCATT", 5, "C-TT", "CATT", {5: (2, 4), 6: (1, 4), 12: (2, 4), 13: (5, 4)}),
    ])
def test_positional_alignment(match_reward, mismatch_penalty, indel_penalty, s, t, expected_score, expected_aligned_s,
                              expected_aligned_t, expected_position_score):
    # Run the alignment
    score, aligned_s, aligned_t, scores_at_positions = positional_alignment(match_reward, mismatch_penalty, indel_penalty, s,
                                                                            t)

    # Assert that the returned score and alignments are as expected
    assert score == expected_score
    assert aligned_s == expected_aligned_s
    assert aligned_t == expected_aligned_t
    assert scores_at_positions == expected_position_score
