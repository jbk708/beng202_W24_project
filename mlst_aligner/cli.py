import click
import time
from mlst_aligner.scoring import GeneScore
from mlst_aligner.utils import subset_fasta
from mlst_aligner.mlst import ScoreMLST


@click.group()
def cli():
    """MLST Aligner CLI."""
    pass


@click.command()
@click.argument('read_fp')
@click.argument('reference')
@click.option('--match', default=2, help='Match score.')
@click.option('--mismatch', default=-2, help='Mismatch penalty.')
@click.option('--indel', default=-1, help='Indel penalty.')
def score(read_fp, reference, match, mismatch, indel):
    """
    Compute and print the gene scores based on alignments.
    """
    start_time = time.time()
    gene_score = GeneScore(read_fp, reference, match=match, mismatch=mismatch, indel=indel)
    final_score = gene_score.get_t_score()
    click.echo(f"Final Score: {final_score}")
    end_time = time.time()
    print(f"Completed in {end_time - start_time:.2f} seconds.")

@click.command()
@click.argument('original_fasta_fp', type=click.Path(exists=True))
@click.option('--subset_count', default=1000, help='Number of reads to include in the subset.', type=int)
@click.argument('output_fasta_fp', type=click.Path())
def subset(original_fasta_fp, subset_count, output_fasta_fp):
    """
    Subsets a FASTA file and saves the subset to a new file.
    
    ORIGINAL_FASTA_FP: Path to the original FASTA file.
    
    OUTPUT_FASTA_FP: Path where the subset FASTA file will be saved.
    """
    subset_fasta(original_fasta=original_fasta_fp, subset_count=subset_count, output_fasta=output_fasta_fp)

@click.command()
@click.argument('reads_fp', type=click.Path(exists=True))
@click.argument('mlst_fp', type=click.Path(exists=True))
@click.option('--match', default=2, help='Match score.')
@click.option('--mismatch', default=-2, help='Mismatch penalty.')
@click.option('--indel', default=-1, help='Indel penalty.')
def score_mlst(reads_fp, mlst_fp, match, mismatch, indel):
    """
    Compute and print the MLST scores for multiple genes based on alignments.
    """
    start_time = time.time()
    mlst_scorer = ScoreMLST(reads_fp=reads_fp, references_fp=mlst_fp, match=match, mismatch=mismatch, indel=indel)
    gene_scores = mlst_scorer.score_mlst()
    for gene_name, gene_score in gene_scores:
        click.echo(f"Gene: {gene_name}, Score: {gene_score}")
    end_time = time.time()
    print(f"Completed in {end_time - start_time:.2f} seconds.")
cli.add_command(score)
cli.add_command(subset)
cli.add_command(score_mlst)
if __name__ == '__main__':
    
    cli()
