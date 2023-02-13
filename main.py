import argparse
import logging
import pprint

from ncbi.datasets.openapi.model.v1_assembly_dataset_request import V1AssemblyDatasetRequest
from analysis import analyse_protein
from analysis import get_top_proteins_and_validate
from analysis import top_by_relevance

from processing import align
from storage import read_correspondences


logging.basicConfig(filename='logs.log', level=logging.DEBUG)


parser = argparse.ArgumentParser(description="Attempts to identify genes relevant " \
+ "for some phenotype by comparing the proteomes of related species.")

# Create a subcommand parser
subparsers = parser.add_subparsers(title="align", dest="subcommand")

# Create a parser for the "align" subcommand
parser_align = subparsers.add_parser("align", help="Aligns the complete proteomes of all pairs of species in the subtree under the given taxon id.")

# Add arguments to the "align" subcommand parser
parser_align.add_argument('taxon', help='The taxon id of the subtree parent node')
parser_align.add_argument("-w", "--overwrite", help='Overwrite existing data.', action="store_true")
parser_align.add_argument("-a", "--additional_genomes", help='Additional genome accessions to include (as comma-separated list)')


# Create a parser for the "analyse" subcommand
parser_align = subparsers.add_parser("analyse", help="Analyse the already generated correspendences.")

# Add arguments to the "analyse" subcommand parser
parser_align.add_argument('-t', '--top', help='Find the top N most relevant proteins.', default=10)
parser_align.add_argument('-p', '--protein', help='Analyse specefic protein by accession.')
parser_align.add_argument('-v', '--validate', help='Run analysis on N-1 species and test on the remaining one', action='store_true')

if __name__ == '__main__':
    args = parser.parse_args()

    if args.subcommand == 'align':
        align(args)
    elif args.subcommand == 'analyse':
        if args.protein:
            analyse_protein(args.protein)
        else:
            correspondences = read_correspondences()
            if args.validate:
                get_top_proteins_and_validate(correspondences, args)
            else:
                top_by_relevance(correspondences, int(args.top))
    else:
        parser.print_help()