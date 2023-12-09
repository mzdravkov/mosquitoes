import logging
import os

from genomes import get_accessions

from proteins import download_dataset, extract_sequence_data, has_sequence_data
from proteins import get_homologs_table
from storage import CORRESPONDENCES_DIR, SEQUENCES_DIR, get_genome_pair_filename, save_protein_correspondences
from taxonomy import get_species_in_taxon_subtree
from utils import print_table


def __get_genome_combinations(genomes):
    combinations = []
    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):
            if genomes[i] != genomes[j]:
                combinations.append((genomes[i], genomes[j]))
    return combinations


def get_genome_pairs_for_processing(genomes, overwrite=False):
    """
    Returns a list of genome pairs that are not yet processed
    (except if overwrite is True, in which case all pairs will be returned).
    """
    # sort each pair to match the filename format
    genome_combinations = (tuple(sorted(pair)) for pair in __get_genome_combinations(genomes))

    if overwrite:
        return list(genome_combinations)

    # if not overwriting, only return the genome pairs that are not yet processed
    processed_genome_pairs = set()
    for filename in os.listdir(CORRESPONDENCES_DIR):
        if filename.endswith('.csv'):
            genome1, genome2 = filename.rstrip('.csv').split('-')
            processed_genome_pairs.add((genome1, genome2))
    return {pair for pair in genome_combinations if pair not in processed_genome_pairs}


def download_sequences_data(accessions, genome=False):
    """
    Downloads data for the given genomes.
    """
    downloaded_data = set()
    for filename in os.listdir(SEQUENCES_DIR):
        if filename.endswith('.faa.zip'):
            downloaded_data.add(filename)
    for accession in accessions:
        if accession + '.faa.zip' not in downloaded_data:
            logging.info('Downloading data for {}'.format(accession))
            download_dataset(accession, genome)
            
            
def create_homologs_table(specie1, specie2, sensitivity, genome):
    homologs, avg_identity = get_homologs_table(specie1, specie2, sensitivity, genome)
    filename = get_genome_pair_filename(specie1, specie2)
    save_protein_correspondences(homologs, filename)
    print('Average identity for {}-{}: {}'.format(specie1, specie2, avg_identity), flush=True)
    logging.info('Average identity for {}-{}: {}'.format(specie1, specie2, avg_identity))


def align(args):
    specie_taxons = list(get_species_in_taxon_subtree(args.taxon))
    print('Got {} taxons in the subtree under taxon {}'.format(len(specie_taxons), args.taxon))

    accessions = get_accessions(specie_taxons)

    print("Downloading data for {} genomes...".format(len(accessions)))
    sequence_type = 'gene' if args.genome else 'protein'
    download_sequences_data(accessions.values(), sequence_type)

    status_table = []
    assemblies = []
    # extract the datasets to get the sequence fasta files
    # and filter out the accessions that don't have sequence data
    for taxon_id, taxon_name, taxon_rank in specie_taxons:
        accession = accessions.get(taxon_id)
        if accession:
            extract_sequence_data(accession, sequence_type)
            if has_sequence_data(accession):
                assemblies.append(accession)
                status_table.append((taxon_id, taxon_name, taxon_rank, accession, 'OK'))
            else:
                logging.warning('No sequence data for {}'.format(accession))
                status_table.append((taxon_id, taxon_name, taxon_rank, accession, 'N/A'))
        else:
            status_table.append((taxon_id, taxon_name, taxon_rank, 'N/A', 'N/A'))
            
    if args.additional_genomes:
        additional_genomes = args.additional_genomes.split(',')
        print('Got {} additional genomes'.format(len(additional_genomes)))
        for accession in additional_genomes:
            extract_sequence_data(accession, sequence_type)
            if has_sequence_data(accession):
                assemblies.append(accession)
                status_table.append(('', '', '', accession, 'OK'))
            else:
                logging.warning('No sequence data for {}'.format(accession))
                status_table.append(('', '', '', accession, 'N/A'))

    status_table.sort(key=lambda x: x[1])

    print_table(status_table, columns=('Taxon', 'Name', 'Rank', 'Accession', 'Protein'))

    print('Got {} accessions with sequence data in their assembly'.format(len(assemblies)))
    logging.info('Got {} accessions with sequence data in their assembly'.format(len(assemblies)))

    genome_pairs = get_genome_pairs_for_processing(assemblies, args.overwrite)

    logging.info('Got {} genome pairs for processing'.format(len(genome_pairs)))
    print('Got {} genome pairs for processing'.format(len(genome_pairs)))

    for specie1, specie2 in genome_pairs:
        create_homologs_table(specie1, specie2, args.sensitivity, args.genome)
