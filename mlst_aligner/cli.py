import click
from mlst_aligner.scoring import GeneScore
from mlst_aligner.utils import subset_fasta

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
    gene_score = GeneScore(read_fp, reference, match=match, mismatch=mismatch, indel=indel)
    gene_score.get_scores()
    final_score = gene_score.get_t_score()
    click.echo(f"Final Score: {final_score}")

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
cli.add_command(score)
cli.add_command(subset)

if __name__ == '__main__':
    cli()
