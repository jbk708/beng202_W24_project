"""cli.py"""
import click
from mlst_aligner.scoring import GeneScore

@click.group()
def cli():
    """GeneScore CLI"""
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

cli.add_command(score)

if __name__ == '__main__':
    cli()