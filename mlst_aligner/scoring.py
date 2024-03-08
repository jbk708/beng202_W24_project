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

    def __init__(self, read_fp: str, reference: str, **kwargs):
        self.reads = read_fasta(read_fp)
        self.reference = reference
        self.scoring_parameters = (kwargs.get("match", 2), kwargs.get("mismatch", -2), kwargs.get("indel", -1))
        self.scores = None
        
    def get_scores(self):
        score_dicts = []
        for read_name in self.reads.references:
            read_sequence = self.reads.fetch(read_name)
            alignment_score, aligned_read, aligned_reference, scores_at_positions = positional_alignment(
                *self.scoring_parameters, s=read_sequence, t=self.reference
            )
            score_dicts.append((scores_at_positions))  
        self.scores = merge_scores(score_dicts)