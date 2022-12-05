from ncbi.datasets.openapi.model.v1_assembly_dataset_request import V1AssemblyDatasetRequest

import argparse
import logging

from genomes import get_accessions
from processing import download_data_for_genomes
from processing import get_genome_pairs_for_processing
from proteins import get_protein_correspondance_table
from storage import get_genome_pair_filename, save_protein_correspondances
from taxonomy import get_species_in_taxon_subtree


logging.basicConfig(filename='logs.log', level=logging.DEBUG)


def print_table(table):
    for row in table:
        print('\t'.join([str(x) for x in row]))


argparser = argparse.ArgumentParser()
argparser.add_argument('taxon', help='The taxon id of the subtree start node')
argparser.add_argument("-w", "--overwrite", help='Overwrite existing data.', action="store_true")
args = argparser.parse_args()

# specie_taxons = get_species_in_taxon_subtree('7157')
specie_taxons = list(get_species_in_taxon_subtree('7157'))[:5]
print_table(specie_taxons)
print('Got {} taxons'.format(len(specie_taxons)))

accessions = get_accessions(specie_taxons)

download_data_for_genomes(accessions)


proteins = [
    'GCF_000005575.2', # Anopheles gambiae
    'GCF_016920715.1', # Anopheles arabiensis
    'GCF_016801865.1' # Culex pipiens pallens
    ]


# genome_pairs = get_genome_pairs_for_processing(accessions)
genome_pairs = get_genome_pairs_for_processing(proteins)

print(genome_pairs)

for specie1, specie2 in genome_pairs:
    logging.info('Starting to build a correspondance table between {}-{}'.format(specie1, specie2))
    correspondances, avg_identity = get_protein_correspondance_table(specie1, specie2)
    filename = get_genome_pair_filename(specie1, specie2)
    save_protein_correspondances(correspondances, filename)
    print('Average identity score for {}-{}: {}'.format(specie1, specie2, avg_identity))
    logging.info('Average identity score for {}-{}: {}'.format(specie1, specie2, avg_identity))

# protein_correspondances, avg_identity = get_protein_correspondance_table(proteins[0], proteins[2])

# for pc in protein_correspondances:
#     print(*pc)

# save_protein_correspondances(protein_correspondances, proteins[0] + '-' + proteins[2] + '.csv')
    
# print('Average identity score: {}'.format(avg_identity))
