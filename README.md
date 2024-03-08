# Positional Aligner for MLST Based Strain Classification
This package provides tools for Multi-Locus Sequence Typing (MLST) based strain classification using positional alignment, an ~~abomination~~ riff on local alignment. It includes functionalities for aligning genetic sequences to reference strains and computing gene scores, facilitating the identification and classification of bacterial strains.

## Requirements

Before installing and running this package, ensure you have the following installed:

- [Python](https://www.python.org/) 3.10 or higher
- [Poetry](https://python-poetry.org/) for dependency management

This package also depends on several third-party libraries, including:

- Click for creating the CLI
- PySAM for handling FASTA files
- pytest, yapf, and pylint for testing, formatting, and linting

All dependencies will be installed automatically when using Poetry to install the package.


## Installation:

1. Clone git repo into directory of choice
2. Open main repo directory
3. run `poetry install`
4. Run your module of choice with `poetry run {module}`

## Using the CLI for GeneScore

The CLI provides a simple and interactive way to compute gene scores based on local sequence alignment. Here's how to use it:

### Basic Command
```sh
poetry run mlst_aligner score [OPTIONS] READ_FP REFERENCE
```
`READ_FP` is the file path to your reads in FASTA format, and `REFERENCE` is the reference sequence string.

Options:
- `match`: The score to reward when characters match. Default is 2.
- `mismatch`: The penalty (negative score) to assign for character mismatches. Default is -2.
- `indel`: The penalty (negative score) to assign for insertions and deletions. Default is -1.

### Example Usage
To compute gene scores with custom scoring parameters:

```
poetry run mlst_aligner score path/to/reads.fasta "ACTG" --match 3 --mismatch -3 --indel -2
```

This command computes and prints the final gene score based on the provided reads file path, reference sequence, and scoring parameters.

## Code Testing, Formatting, and Linting Standards

For code testing, run from the root folder:

```sh
poetry run pytest
```

For code formatting, simply run [yapf](https://github.com/google/yapf) from the root folder:

```sh
poetry run yapf --in-place --recursive ./mlst_aligner ./tests
```

For linting, simply run [pylint](https://pylint.pycqa.org/en/latest/) from the root folder:
```sh
poetry run pylint ./mlst_aligner ./tests