import logging
import os

from proteins import download_dataset
from storage import CORRESPONDENCES_DIR, SEQUENCES_DIR


def __get_genome_combinations(genomes):
    combinations = []
    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):
            combinations.append((genomes[i], genomes[j]))
    return combinations


def get_genome_pairs_for_processing(genomes):
    """
    Returns a list of genome pairs that are not yet processed.
    """
    processed_genome_pairs = set()
    for filename in os.listdir(CORRESPONDENCES_DIR):
        if filename.endswith('.csv'):
            genome1, genome2 = filename.split('-')
            processed_genome_pairs.add((genome1, genome2))
    # sort each pair to match the filename format
    genome_combinations = (tuple(sorted(pair)) for pair in __get_genome_combinations(genomes))
    return [pair for pair in genome_combinations if pair not in processed_genome_pairs]


def download_data_for_genomes(genomes):
    """
    Downloads data for the given genomes.
    """
    downloaded_data = set()
    for filename in os.listdir(SEQUENCES_DIR):
        if filename.endswith('.faa.zip'):
            downloaded_data.add(filename)
    for genome in genomes:
        if genome + '.faa.zip' not in downloaded_data:
            logging.info('Downloading data for {}'.format(genome))
            download_dataset(genome)