# Metagenomics-and-Read-Mapping

## About

Scripts for DNA Read-mappipng and Metagenomic Analysis

## Prerequisites

Before running the script, ensure you have the following installed:
- Python 3.x
- Required Python packages:

```sh
pip install biopython
```

## Usage

To run the `aligner-final.py` script, open your command line interface and execute the following command:

```sh
& C:/Python313/python.exe "aligner-final.py"
```

Or you can open the file in VSCode and run it using the 'Run Python file' button

To specify which reference file and which reads file to run from, modify refAddress and readsAddress from main, or line 83 and 84.

```python
if __name__ == "__main__":
    refAddress = '.\\reference_genome.fasta'
    readsAddress = '.\\with_error_paired_reads.fasta'
```

The output should be a predictions.csv within the same directory as the code

To run the `meta.py` script, open your command line interface and execute:
    ```sh
    python meta.py
    ```
Or press the run button from VS Code

Change the string value of these two variables in meta.py:
```python
reference_genomes_directory = "./references"
reads_file = "./reads.fasta"
```
In this case, reference genomes are all in the references directory, and the reads file is in the same
directory as the meta.py.
