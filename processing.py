import logging
import os

from multiprocessing import Pool
from genomes import get_accessions

from proteins import download_dataset, extract_protein_data, has_protein_data
from proteins import get_protein_correspondence_table
from storage import CORRESPONDENCES_DIR, SEQUENCES_DIR, get_genome_pair_filename, save_protein_correspondences
from taxonomy import get_species_in_taxon_subtree
from utils import print_table


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
            
            
def create_correspondance_table(args):
    specie1, specie2 = args
    correspondences, avg_identity = get_protein_correspondence_table(specie1, specie2)
    filename = get_genome_pair_filename(specie1, specie2)
    save_protein_correspondences(correspondences, filename)
    print('Average identity score for {}-{}: {}'.format(specie1, specie2, avg_identity), flush=True)
    logging.info('Average identity score for {}-{}: {}'.format(specie1, specie2, avg_identity))


def create_correspondence_tables_in_parallel(genome_pairs):
    """
    Creates a set of correspondence tables for the given
    genome pairs and saves them to the file system.
    """
    with Pool() as p:
        p.map(create_correspondance_table, genome_pairs)
        

def align(args):
    specie_taxons = list(get_species_in_taxon_subtree(args.taxon))
    print('Got {} taxons in the subtree under taxon {}'.format(len(specie_taxons), args.taxon))

    accessions = get_accessions(specie_taxons)

    print("Downloading data for {} genomes...".format(len(accessions)))
    download_data_for_genomes(accessions.values())

    status_table = []
    prot_assemblies = []
    # extract the datasets to get the protein fasta files
    # and filter out the accessions that don't have protein data
    for taxon_id, taxon_name, taxon_rank in specie_taxons:
        accession = accessions.get(taxon_id)
        if accession:
            extract_protein_data(accession)
            if has_protein_data(accession):
                prot_assemblies.append(accession)
                status_table.append((taxon_id, taxon_name, taxon_rank, accession, 'OK'))
            else:
                logging.warning('No protein data for {}'.format(accession))
                status_table.append((taxon_id, taxon_name, taxon_rank, accession, 'N/A'))
        else:
            status_table.append((taxon_id, taxon_name, taxon_rank, 'N/A', 'N/A'))

    status_table.sort(key=lambda x: x[1])

    print_table(status_table, columns=('Taxon', 'Name', 'Rank', 'Accession', 'Protein'))

    # accessions = filtered_accessions
    print('Got {} accessions with protein data in their assembly'.format(len(prot_assemblies)))
    logging.info('Got {} accessions with protein data in their assembly'.format(len(prot_assemblies)))

    # proteins = [
    #     # 'GCF_000005575.2', # Anopheles gambiae
    #     # 'GCF_016920715.1', # Anopheles arabiensis
    #     # 'GCF_016801865.1' # Culex pipiens pallens
    #     # 'GCF_943734845.2',
    #     # 'GCF_013141755.1'
    #     ]
    # genome_pairs = get_genome_pairs_for_processing(proteins)

    genome_pairs = get_genome_pairs_for_processing(prot_assemblies)

    logging.info('Got {} genome pairs for processing'.format(len(genome_pairs)))
    print('Got {} genome pairs for processing'.format(len(genome_pairs)))

    create_correspondence_tables_in_parallel(genome_pairs)
