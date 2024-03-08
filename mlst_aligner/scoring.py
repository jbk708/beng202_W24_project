"""mlst.py"""
from typing import List

from mlst_aligner.aligner import positional_alignment
from mlst_aligner.utils import read_fasta, weighted_average


def merge_scores(scores_at_positions_list: List[Dict[int, int]]) -> Dict[int, List[int]]:
    """
    Merges scoring dictionaries from multiple alignments into a single dictionary,
    where each key is a position, and each value is a list of scores for that position.

    Args:
        scores_at_positions_list (List[Dict[int, int]]): A list of dictionaries, where each dictionary
                                                          represents scores at positions for a single alignment.

    Returns:
        Dict[int, List[int]]: A dictionary where each key is a position, and each value is a list of scores for that position.
    """
    merged_scores = {}

    for scores_at_positions in scores_at_positions_list:
        for position, score in scores_at_positions.items():
            if position not in merged_scores:
                merged_scores[position] = [score]
            else:
                merged_scores[position].append(score)

    return merged_scores


class GeneScore:
    """docstring to come soonTM
    """

    def __init__(self, read_fp: str, reference: str, **kwargs):
        """GeneScore Initialization"""
        self.reads = read_fasta(read_fp)
        self.reference = reference
        self.scoring_parameters = (kwargs.get("match", 2), kwargs.get("mismatch", -2), kwargs.get("indel", -1))
        self.scores = None

    def get_scores(self):
        """
        Fetches each read from the FASTA file, performs local sequence alignment against a reference sequence,
        and merges the scores at each position into a single dictionary.

        This method iterates over each read in the FASTA file specified by the read file path provided during
        object initialization. It performs local sequence alignment of each read against the reference sequence
        using the scoring parameters. The scores at each position from these alignments are then aggregated across 
        all reads to generate ascoring dictionary for the aligned segments of the reference. The resulting merged 
        scores dictionary combines all scores on the position to which they were aligned.

        Attributes:
            scores (Dict[int, List[int]]): A dictionary where each key is a position in the reference sequence,
                                        and each value is a list of scores for that position from all reads'
                                        alignments. This dictionary provides a comprehensive overview of how
                                        each position in the reference sequence aligns with the reads.
        """
        score_dicts = []
        for read_name in self.reads.references:
            read_sequence = self.reads.fetch(read_name)
            alignment_score, aligned_read, aligned_reference, scores_at_positions = positional_alignment(
                *self.scoring_parameters, s=read_sequence, t=self.reference)
            score_dicts.append((scores_at_positions))
        self.scores = merge_scores(score_dicts)
