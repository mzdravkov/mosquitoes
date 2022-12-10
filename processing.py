import logging
import os

from multiprocessing import Pool

from proteins import download_dataset
from proteins import get_protein_correspondence_table
from storage import CORRESPONDENCES_DIR, SEQUENCES_DIR, get_genome_pair_filename, save_protein_correspondences


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
            
            
def __get_protein_correspondence_table_wrapper(args):
    specie1, specie2 = args
    correspondences, avg_identity = get_protein_correspondence_table(specie1, specie2)
    return (specie1, specie2), correspondences, avg_identity


def get_correspondence_tables_in_parallel(genome_pairs):
    """
    Returns a list of correspondence tables for the given genome pairs.
    """
    with Pool() as p:
        results = p.map(__get_protein_correspondence_table_wrapper, genome_pairs)
    
    for (specie1, specie2), correspondences, avg_identity in results:
        filename = get_genome_pair_filename(specie1, specie2)
        save_protein_correspondences(correspondences, filename)
        print('Average identity score for {}-{}: {}'.format(specie1, specie2, avg_identity))
        logging.info('Average identity score for {}-{}: {}'.format(specie1, specie2, avg_identity))