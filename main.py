import argparse
import logging

from ncbi.datasets.openapi.model.v1_assembly_dataset_request import V1AssemblyDatasetRequest

from genomes import get_accessions
from processing import download_data_for_genomes, get_correspondence_tables_in_parallel
from processing import get_genome_pairs_for_processing
from proteins import extract_protein_data, get_protein_correspondence_table, has_protein_data
from storage import get_genome_pair_filename, save_protein_correspondences
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
specie_taxons = list(get_species_in_taxon_subtree('7157'))
print_table(specie_taxons)
print('Got {} taxons'.format(len(specie_taxons)))

accessions = get_accessions(specie_taxons)

download_data_for_genomes(accessions)

filtered_accessions = []
# extract the datasets to get the protein fasta files
# and filter out the accessions that don't have protein data
for accession in accessions:
    extract_protein_data(accession)
    if has_protein_data(accession):
        filtered_accessions.append(accession)
    else:
        logging.warning('No protein data for {}'.format(accession))
accessions = filtered_accessions

logging.info('Got {} accessions with protein data.'.format(len(accessions)))

proteins = [
    'GCF_000005575.2', # Anopheles gambiae
    'GCF_016920715.1', # Anopheles arabiensis
    # 'GCF_016801865.1' # Culex pipiens pallens
    'GCF_943734665.1'
    ]


# genome_pairs = get_genome_pairs_for_processing(accessions)

genome_pairs = get_genome_pairs_for_processing(proteins)
print(genome_pairs)

get_correspondence_tables_in_parallel(genome_pairs)


# for specie1, specie2 in genome_pairs:
#     logging.info('Starting to build a correspondence table between {}-{}'.format(specie1, specie2))
#     correspondences, avg_identity = get_protein_correspondence_table(specie1, specie2)


# scores = []
# for i in range(10):
#     correspondences, avg_identity = get_protein_correspondence_table(proteins[1], proteins[0], subset=100)
#     scores.append(avg_identity)
#     print('Average identity score for {}-{}: {}'.format(proteins[1], proteins[0], avg_identity))
# print("TOTAL average = ", sum(scores) / len(scores))


# protein_correspondences, avg_identity = get_protein_correspondence_table(proteins[0], proteins[2])

# for pc in protein_correspondences:
#     print(*pc)

# save_protein_correspondences(protein_correspondences, proteins[0] + '-' + proteins[2] + '.csv')

# print('Average identity score: {}'.format(avg_identity))
