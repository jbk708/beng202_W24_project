"""mlst.py"""
from mlst_aligner.utils import fetch_references
from mlst_aligner.scoring import GeneScore


class ScoreMLST(GeneScore):
    """
    A class for scoring MLST (Multi-Locus Sequence Typing) based on gene sequences against a set of reference sequences.

    Attributes:
        reads_fp (str): File path to the FASTA file containing read sequences.
        references_fp (str): File path to the FASTA file containing reference sequences.
        references (List[Tuple[str, str]]): A list of tuples, where each tuple contains a gene name and its sequence,
                                             extracted from the references FASTA file.
    
    Inherits:
        GeneScore: Inherits from the GeneScore class to utilize its scoring mechanisms.
    
    Args:
        reads_fp (str): File path to the reads FASTA file.
        references_fp (str): File path to the references FASTA file.
        **kwargs: Arbitrary keyword arguments passed to the GeneScore initializer.
    """

    def __init__(self, reads_fp: str, references_fp: str, **kwargs):
        """
        Initializes ScoreMLST
        """
        super().__init__(reads_fp, "", **kwargs)
        self.reads_fp = reads_fp
        self.references = fetch_references(references_fp)

    def score_mlst(self):
        """
        Scores each gene in the references against the reads, calculating a total score for each gene.

        Iterates through each reference gene, sets it as the current reference in the superclass, and calculates
        the total score for that gene using the superclass's scoring mechanism.

        Returns:
            List[Tuple[str, int]]: A list of tuples, where each tuple contains a gene name and its corresponding
                                   total score. The scores are computed based on alignments with the reads.
        """
        gene_scores = []
        for gene_name, sequence in self.references:
            self.reference = sequence
            gene_score = self.get_t_score()
            gene_scores.append((gene_name, gene_score))

        return gene_scores
