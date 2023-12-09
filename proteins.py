import csv
import logging
import os
import shutil
import subprocess
import sys

from os.path import join

from Bio.SeqIO.FastaIO import SimpleFastaParser
from ncbi.datasets import GenomeApi
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
# from ncbi.datasets.package import dataset

from ncbi.datasets.openapi.model.v1_annotation_for_assembly_type import V1AnnotationForAssemblyType

from storage import SEQUENCES_DIR, get_diamond_db_filename, get_sequences_filename
from storage import get_mmseqs_db_filename
from storage import get_mmseqs_alignment_db_filename
from storage import get_mmseqs_results_filename
from storage import get_dataset_filename
from storage import TMP_DIR


def diamond_align(specie1, specie2, reverse=False, sensitivity='sensitive'):
    specie1_file = get_sequences_filename(specie1)
    specie2_file = get_sequences_filename(specie2)

    # prepare diamond-formatted DB file for specie2 if missing
    specie2_db_file = get_diamond_db_filename(specie2)
    if not os.path.isfile(specie2_db_file):
        makedb = ['diamond', 'makedb', '--in', specie2_file, '--db', specie2_db_file]
        subprocess.run(makedb, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # execute diamond as a subprocess
    cmd = """
    diamond blastp --db {db}
                   -q {query}
                   --max-target-seqs 1
                   --{sensitivity}
                   -f 6 qseqid sseqid pident

    """.format(db=specie2_db_file, query=specie1_file, sensitivity=sensitivity)
    proc = subprocess.run(cmd.split(), capture_output=True, text=True)

    # parse the output of diamond and get a correspondence_table
    correspondences = {}
    for line in proc.stdout.splitlines():
        query, hit, identity = line.rstrip().split('\t')
        key = (query, hit)
        if reverse:
            key = (hit, query)
        correspondences[key] = float(identity)

    return correspondences


def make_mmseqs_db(input, output):
    if not os.path.isfile(output):
        makedb = ['mmseqs', 'createdb', input, output]
        subprocess.run(makedb, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def mmseqs_align(species1, species2, reverse=False, sensitivity='sensitive'):
    print('Aligning {} to {}'.format(species1, species2))
    species1_file = get_sequences_filename(species1)
    species2_file = get_sequences_filename(species2)

    # prepare mmseqs-formatted DB files for species1 if missing
    species1_db = get_mmseqs_db_filename(species1)
    print(species1_file)
    print(species1_db)
    make_mmseqs_db(species1_file, species1_db)

    # prepare mmseqs-formatted DB files for species2 if missing
    species2_db = get_mmseqs_db_filename(species2)
    make_mmseqs_db(species2_file, species2_db)

    alignment_db = get_mmseqs_alignment_db_filename(species1, species2)

    # execute mmseqs as a subprocess
    cmd = """
    mmseqs search --search-type 3
                  {query}
                  {target}
                  {alignment}
                  {tmp}
    """.format(query=species1_db, target=species2_db, alignment=alignment_db, tmp=TMP_DIR)
    subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    result_file = get_mmseqs_results_filename(species1, species2)

    cmd = """
    mmseqs createtsv {query} {target} {alignment} {result}
    """.format(query=species1_db, target=species2_db, alignment=alignment_db, result=result_file)
    subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # parse the output of mmseqs and get a homologs table
    homologs = {}
    with open(result_file) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            query, hit, _ , identity = row[:4]
            key = (query, hit)
            if reverse:
                key = (hit, query)
            homologs[key] = float(identity)

    return homologs


def add_missing_homologs(homologs, specie1, specie2):
    """
    Adds sequences that are not already present in the homologs table.
    """
    all_specie1_seqs = {row[0] for row in homologs}
    with open(get_sequences_filename(specie1)) as handle:
        for header, _ in SimpleFastaParser(handle):
            id = header.split(' ')[0]
            if id not in all_specie1_seqs:
                homologs.append((id, None, 0))

    # all_specie2_seqs = {row[1] for row in homologs}
    # with open(get_sequences_filename(specie2)) as handle:
    #     for header, _ in SimpleFastaParser(handle):
    #         id = header.split(' ')[0]
    #         if id not in all_specie2_seqs:
    #             homologs.append((None, id, 0))


def get_homologs_table(specie1, specie2, sensitivity, genome):
    """
    Takes the acession ids of two species and returns a
    homolog table and an average identity percentage.
    The homolog table contains records of the following type:
    (specie1_seq_id, specie2_seq_id, identity)
    """
    aligner = mmseqs_align if genome else diamond_align
    # Search specie1 sequences in specie2's genome/proteome
    homologs_forw = aligner(specie1, specie2, sensitivity=sensitivity)
    # Search specie2 sequences in specie1's genome/proteome
    # homologs_back = aligner(specie2, specie1, reverse=True)

    homologs = []
    total_sum = 0
    # # Merge the two results into a single table
    # all_pairs = set(homologs_forw.keys()).union(homologs_back.keys())
    # for pair in all_pairs:
    #     score = homologs_forw.get(pair, homologs_back.get(pair))
    #     if pair in homologs_forw and pair in homologs_back:
    #         # Should be the same in both places, but just in case
    #         score = (homologs_forw[pair] + homologs_back[pair])/2
    #     homologs.append((pair[0], pair[1], score))
    #     total_sum += score

    for pair in homologs_forw:
        identity = homologs_forw[pair]
        homologs.append((pair[0], pair[1], identity))
        total_sum += identity

    # if there are sequences in specie1 that don't have a match in specie2 or vice versa,
    # add them to the table with a score of 0.
    add_missing_homologs(homologs, specie1, specie2)

    avg_identity = total_sum/len(homologs)

    return homologs, avg_identity


def download_dataset(accession, sequence_type="protein"):
    """
    Downloads the protein fasta file for the genome assembly with the given accession number.

    Arguments:
    accession -- species' accession number

    Keyword arguments:
    sequence_type -- protein/gene (default: protein)
    """
    if sequence_type not in ('protein', 'gene'):
        raise ValueError(
            'sequence_type must be "protein" or "gene". Got {}'.format(sequence_type))

    zipfile_name = get_dataset_filename(accession)

    # Check if zipfile_name already exists, skipping the download if it does.
    try:
        with open(zipfile_name, 'r') as f:
            logging.info('File %s already exists. Skipping.', zipfile_name)
            return
    except FileNotFoundError:
        pass

    # download the data package through the api-client.
    with DatasetsApiClient() as api_client:
        genome_api = GenomeApi(api_client)
        try:
            annotation_type = 'CDS_FASTA' if sequence_type == 'gene' else 'PROT_FASTA'
            exclude_genome = sequence_type == 'protein'
            gene_dataset_download = genome_api.download_assembly_package(
                [accession],
                include_annotation_type=[V1AnnotationForAssemblyType(annotation_type)],
                exclude_sequence=exclude_genome,
                filename=zipfile_name,
                _preload_content=False,
            )

            with open(zipfile_name, "wb") as f:
                f.write(gene_dataset_download.data)
        except DatasetsApiException as e:
            sys.exit(f"Exception when calling GeneApi: {e}\n")


def extract_sequence_data(accession, sequence_type):
    """
    Checks if the archive with data for a genome is already
    downloaded and contains sequence data of the given type.
    If it does, it will extract the archive and return the
    path to the sequences file.

    Arguments:
    accession     -- species' accession number
    sequence_type -- type of sequence data. Should be "protein" or "gene".
    """
    if sequence_type == 'protein':
        file_type = 'PROTEIN_FASTA'
    else:
        # file_type = 'GENOMIC_NUCLEOTIDE_FASTA'
        file_type = 'CDS_NUCLEOTIDE_FASTA'

    sequences_target_filepath = get_sequences_filename(accession)

    # If the data is already extracted, do nothing
    if os.path.isdir(sequences_target_filepath):
        return

    # If the data is not downloaded, do nothing
    zipfile = get_dataset_filename(accession)
    if not os.path.isfile(zipfile):
        return

    # Don't extract if there is no sequence data
    package = dataset.GeneDataset(zipfile)
    sequence_files = package.get_file_names_by_type(file_type)
    print(sequence_files)
    if len(sequence_files) == 0:
        return

    tmp_dir = join(SEQUENCES_DIR, 'tmp')

    shutil.unpack_archive(zipfile, tmp_dir, format='zip')

    # move the extracted data to the sequence directory
    # src = join(tmp_dir, 'ncbi_dataset', 'data', accession)
    # TODO: check if we may have more than one file
    local_path = next(file for file in sequence_files if file.startswith(accession))
    src = join(tmp_dir, 'ncbi_dataset', 'data', local_path)
    dest = join(SEQUENCES_DIR, accession + '.faa')
    shutil.move(src, dest)

    # clean up the directory left from extracting the protein data
    shutil.rmtree(tmp_dir)


def has_sequence_data(accession):
    """
    Checks if the sequence fasta file for the given accession is available.
    """
    fasta = get_sequences_filename(accession)
    return os.path.isfile(fasta)


def get_protein_sequence(protein, species):
    """
    Takes a protein accession number and genome accession number
    and returns the sequence for the protein as found in the
    local sequences storage.
    """
    fasta_file = get_sequences_filename(species)
    with open(fasta_file) as handle:
        for header, sequence in SimpleFastaParser(handle):
            if protein in header:
                return sequence
    return None