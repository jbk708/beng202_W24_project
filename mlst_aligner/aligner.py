"""aligner.py"""
from typing import Tuple, Dict


def positional_alignment(match_reward: 2, mismatch_penalty: -2, indel_penalty: 1, s: str,
                         t: str) -> Tuple[int, str, str, Dict[int, int]]:
    """
    Perform local sequence alignment between two strings using dynamic programming.

    This function computes the local alignment between string s (source) and t (target) based on the given scoring parameters. It returns the highest alignment score, the corresponding local alignment for s and t, and a dictionary containing the scores at each position in the target string where an actual alignment (match/mismatch) occurred.

    Args:
        match_reward (int): The score to reward when characters match.
        mismatch_penalty (int): The penalty (negative score) to assign for character mismatches.
        indel_penalty (int): The penalty (negative score) to assign for insertions and deletions.
        s (str): The source string to align.
        t (str): The target string to align with the source string.

    Returns:
        max_score (int): The highest score achieved in the alignment.
        aligned_s (str): The aligned version of the source string with gaps ('-') as necessary.
        aligned_t (str): The aligned version of the target string with gaps ('-') as necessary.
        scores_at_positions (Dict[int, int]): A dictionary mapping each position in the target string where an actual alignment occurred to the score achieved at that position.

    Examples:
        >>> local_alignment_w_scores(1, -1, -1, 'GATTACA', 'GCATGCU')
        (3, 'GATTA', 'GAT--', {1: 0, 2: 1, 3: 2, 5: 3})
    """
    dp = {0: {0: 0}}
    backtrack = {0: {0: 0}}
    scores_at_positions = {}

    for i in range(1, len(s) + 1):
        dp[i] = {0: 0}
        backtrack[i] = {0: 0}
    for j in range(1, len(t) + 1):
        dp[0][j] = 0
        backtrack[0][j] = 0

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, len(s) + 1):
        for j in range(1, len(t) + 1):
            match = match_reward if s[i - 1] == t[j - 1] else -mismatch_penalty
            scores = {0: 0, 1: dp[i - 1][j] - indel_penalty, 2: dp[i][j - 1] - indel_penalty, 3: dp[i - 1][j - 1] + match}
            dp[i][j] = max(scores.values())
            backtrack[i][j] = max(scores, key=scores.get)

            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_pos = (i, j)

            scores_at_positions[j] = dp[i][j]

    i, j = max_pos
    aligned_s, aligned_t = '', ''
    while i > 0 and j > 0 and dp[i][j] > 0:
        if backtrack[i][j] == 3:
            aligned_s = s[i - 1] + aligned_s
            aligned_t = t[j - 1] + aligned_t
            i -= 1
            j -= 1
        elif backtrack[i][j] == 1:
            aligned_s = s[i - 1] + aligned_s
            aligned_t = '-' + aligned_t
            i -= 1
        elif backtrack[i][j] == 2:
            aligned_s = '-' + aligned_s
            aligned_t = t[j - 1] + aligned_t
            j -= 1
        else:
            break

    return max_score, aligned_s, aligned_t, scores_at_positions
