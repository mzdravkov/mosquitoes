import argparse
import logging

from analysis import analyse_protein
from analysis import moving_window_validation
from analysis import get_top_proteins_and_validate
from analysis import top_by_relevance

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
parser_align.add_argument("-a", "--additional_genomes", help='Additional genome accessions to include (as comma-separated list)')


# Create a parser for the "analyse" subcommand
parser_align = subparsers.add_parser("analyse", help="Analyse the already generated correspendences.")

# Add arguments to the "analyse" subcommand parser
parser_align.add_argument('-t', '--top', help='Find the top N most relevant proteins.', type=int, default=10)
parser_align.add_argument('-p', '--protein', help='Analyse specefic protein by accession.')
parser_align.add_argument('-v', '--validate', help='Run analysis on N-1 species and test on the remaining one', action='store_true')
parser_align.add_argument('-w', '--window', help='Run validation on a window with size W across the top N proteins', action='store_true')
parser_align.add_argument('-W', '--window_size', help='Size of the moving window for the validation', type=int, default=5)
parser_align.add_argument('-s', '--window_step', help='Window step', type=int, default=1)

if __name__ == '__main__':
    args = parser.parse_args()

    if args.subcommand == 'align':
        align(args)
    elif args.subcommand == 'analyse':
        if args.protein:
            analyse_protein(args.protein)
        else:
            if args.validate:
                if args.window:
                    moving_window_validation(args.top, args.window_size, args.window_step)
                else:
                    get_top_proteins_and_validate(args.top)
            else:
                top_by_relevance(args.top)
    else:
        parser.print_help()