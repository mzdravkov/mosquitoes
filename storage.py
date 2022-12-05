import csv
from os.path import join


DATA_DIR = 'data'
CORRESPONDANCES_DIR = join(DATA_DIR, 'correspondances')
SEQUENCES_DIR = join(DATA_DIR, 'sequences')


def read_protein_correspondances(filename):
    """
    Read correspondances data from a CSV file.
    """
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        return list(reader)


def save_protein_correspondances(data, filename):
    """
    Save correspondances data as a CSV file.
    """
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['specie1_protein_id', 'specie2_protein_id', 'identity'])
        writer.writerows(data)


def get_genome_pair_filename(genome1, genome2):
    """
    Returns the name of the correspondance file for a pair of genomes.
    """
    filename = '{}-{}.csv'.format(*sorted([genome1, genome2]))
    return join(CORRESPONDANCES_DIR, filename)


def get_dataset_filename(accession):
    return join(SEQUENCES_DIR, accession + '.faa.zip')