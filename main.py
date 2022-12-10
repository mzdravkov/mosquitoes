import argparse
import logging

from ncbi.datasets.openapi.model.v1_assembly_dataset_request import V1AssemblyDatasetRequest

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


if __name__ == '__main__':
    args = parser.parse_args()

    if args.subcommand == 'align':
        align(args)
    else:
        parser.print_help()