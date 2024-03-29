import csv
import os

import pandas as pd
import numpy as np

from os.path import join

# from ncbi.datasets.package import dataset


DATA_DIR = 'data'
CORRESPONDENCES_DIR = join(DATA_DIR, 'correspondences')
SEQUENCES_DIR = join(DATA_DIR, 'sequences')
DOWNLOADED_DATA_DIR = join(SEQUENCES_DIR, 'downloads')
TMP_DIR = join(DATA_DIR, 'tmp')


def read_protein_correspondences(filename):
    """
    Read correspondences data from a CSV file.
    """
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        # skip header
        try:
            next(reader)
        except StopIteration:
            return []
        return list(reader)


def save_protein_correspondences(data, filename):
    """
    Save correspondences data as a CSV file.
    """
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['specie1_protein_id', 'specie2_protein_id', 'identity'])
        writer.writerows(data)


def get_genome_pair_filename(genome1, genome2):
    """
    Returns the name of the correspondence file for a pair of genomes.
    """
    filename = '{}-{}.csv'.format(*sorted([genome1, genome2]))
    return join(CORRESPONDENCES_DIR, filename)


def get_dataset_filename(accession):
    return join(DOWNLOADED_DATA_DIR, accession + '.faa.zip')


def get_sequences_filename(accession):
    return join(SEQUENCES_DIR, accession + '.faa')


def get_diamond_db_filename(accession):
    return join(SEQUENCES_DIR, accession + '.dmnd')


def get_mmseqs_db_filename(accession):
    return join(DATA_DIR, 'mmseqs_dbs', accession)


def get_mmseqs_alignment_db_filename(accession1, accession2):
    return join(TMP_DIR, 'alignments', accession1 + '_' + accession2)


def get_mmseqs_results_filename(accession1, accession2):
    return join(TMP_DIR, 'alignments', accession1 + '_' + accession2 + '.tsv')


def read_correspondences():
    """
    Read all correspondence data from the file system.
    """
    correspondence_arrs = {}
    for filename in os.listdir(CORRESPONDENCES_DIR):
        if filename.endswith('.csv'):
            genome1, genome2 = filename.rstrip('.csv').split('-')
            key = (genome1, genome2)
            correspondences = read_protein_correspondences(join(CORRESPONDENCES_DIR, filename))
            correspondence_arrs[key] = correspondences
    return correspondence_arrs


def read_correspondences_df():
    """
    Read all correspondence data from the file system as a Pandas DataFrame.
    """
    correspondences_dfs = []
    for filename in os.listdir(CORRESPONDENCES_DIR):
        if filename.endswith('.csv'):
            species1, species2 = filename.rstrip('.csv').split('-')
            df = pd.read_csv(join(CORRESPONDENCES_DIR, filename))
            row_count = df.shape[0]
            df['species1'] = np.repeat(species1, row_count)
            df['species2'] = np.repeat(species2, row_count)
            correspondences_dfs.append(df)
    merged = pd.concat(correspondences_dfs)
    return merged.reset_index()
