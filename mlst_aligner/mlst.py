"""mlst.py"""
from mlst_aligner.utils import fetch_references
from mlst_aligner.scoring import GeneScore


class ScoreMLST(GeneScore):

    def __init__(self, reads_fp: str, references_fp: str, **kwargs):
        super().__init__(reads_fp, "", **kwargs)
        self.reads_fp = reads_fp
        self.references = fetch_references(references_fp)

    def score_mlst(self):
        gene_scores = []
        for gene_name, sequence in self.references:
            self.reference = sequence
            gene_score = self.get_t_score()
            gene_scores.append((gene_name, gene_score))

        return gene_scores
