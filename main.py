import argparse
import logging

from ncbi.datasets.openapi.model.v1_assembly_dataset_request import V1AssemblyDatasetRequest
from analysis import analyse_protein, top_by_relevance

from processing import align


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


# Create a parser for the "analyse" subcommand
parser_align = subparsers.add_parser("analyse", help="Analyse the already generated correspendences.")

# Add arguments to the "analyse" subcommand parser
parser_align.add_argument('-t', '--top', help='Find the top N most relevant proteins.', default=10)
parser_align.add_argument('-p', '--protein', help='Analyse specefic protein by accession.')

if __name__ == '__main__':
    args = parser.parse_args()

    if args.subcommand == 'align':
        align(args)
    elif args.subcommand == 'analyse':
        if args.protein:
            analyse_protein(args.protein)
        else:
            top_by_relevance(int(args.top))
    else:
        parser.print_help()