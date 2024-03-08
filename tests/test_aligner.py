"""test_aligner.py"""
import os
import pytest
from mlst_aligner.aligner import positional_alignment


@pytest.mark.parametrize(
    "match_reward, mismatch_penalty, indel_penalty, s, t, expected_score, expected_aligned_s, expected_aligned_t, expected_position_score", [
        (3, 3, 1, "AGC", "ATC", 4, "A-GC", "AT-C", {
            1: 1,
            2: 0,
            3: 4
        }),
        (1, 1, 1, "TAACG", "ACGTG", 3, "ACG", "ACG", {
            1: 0,
            2: 1,
            3: 3,
            4: 2,
            5: 1
        }),
        (3, 2, 1, "CAGAGATGGCCG", "ACG", 6, "CG", "CG", {
            1: 0,
            2: 2,
            3: 6
        }),
        (2, 3, 1, "CTT", "AGCATAAAGCATT", 5, "C-TT", "CATT", {
            1: 0,
            2: 0,
            3: 0,
            4: 0,
            5: 2,
            6: 1,
            7: 0,
            8: 0,
            9: 0,
            10: 0,
            11: 0,
            12: 2,
            13: 5
        }),
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
